## =============================================================================
## 34_fit_mortality.R
##
## FVS-CONUS: Individual tree mortality driver for mort_gompit.stan
## (complementary log-log / gompit link with log(YEARS) exposure offset,
## nested EPA L1 / L2 / L3 + species random intercept hierarchy, harmonized
## with DG / HG / HCB).
##
## The cloglog link derives from a continuous-time proportional hazards model:
##   p_die = 1 - exp(-exp(eta))
## where eta includes a log(YEARS) offset so that the linear predictor on the
## log-hazard scale works with variable remeasurement intervals.
##
## Response:
##   alive = 1 if STATUS2 == 1 (survived), 0 if died from natural causes
##   (AGENTCD < 80 to exclude harvest mortality).
##
## Inputs:
##   calibration/data/conus_remeasurement_pairs.rds (metric, per 30 Step 5b)
##   calibration/traits/species_traits.rds          (optional, for trait prior)
##
## Usage:
##   Rscript scripts/34_fit_mortality.R --site cspi
##   Rscript scripts/34_fit_mortality.R --site bgi --traits 0
##
## Author: A. Weiskittel, 2026-04-16
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(optparse)
})

## Shared trait columns (must match 31 / 32 / 33 / 36 for consistency)
DEFAULT_TRAIT_COLS <- c(
  "wood_specific_gravity",
  "shade_tolerance_num",
  "softwood",
  "leaf_longevity_months",
  "sla_mm2_mg",
  "seed_mass_mg",
  "max_ht_m",
  "max_dbh_cm",
  "vulnerability_score",
  "sensitivity"
)

## ---- 1. Data prep ----------------------------------------------------------
prepare_mort_data <- function(dat, site_var = "cspi", min_sp_n = 5000L) {
  message("Preparing mortality data. site_var = ", site_var,
          "  min_sp_n = ", min_sp_n)
  setDT(dat)

  ## Mortality flag: alive if STATUS2 == 1; dead from natural causes if
  ## STATUS2 == 2 and AGENTCD2 < 80 (exclude harvest). Trees harvested
  ## (AGENTCD2 >= 80) or with missing status are dropped.
  dat <- dat[
    STATUS1 == 1 &
    !is.na(STATUS2) &
    (STATUS2 == 1 | (STATUS2 == 2 & !is.na(AGENTCD2) & AGENTCD2 < 80))
  ]
  dat[, alive := as.integer(STATUS2 == 1)]

  dat <- dat[
    !is.na(DBH1) & DBH1 >= 1.0 &
    !is.na(CR1)  & CR1 > 0 & CR1 <= 1.0 &
    !is.na(BA1)  & BA1 > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    YEARS >= 1 & YEARS <= 20 &
    !is.na(get(site_var)) & get(site_var) > 0
  ]

  ## Rare species pooling
  sp_n <- dat[, .N, by = SPCD]
  rare <- sp_n[N < min_sp_n, SPCD]
  dat[SPCD %in% rare, SPCD := -1L]

  ## Derived covariates matching mort_gompit.stan
  cspi_shift <- max(1.0 - min(dat[[site_var]], na.rm = TRUE), 0.01)
  dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]
  dat[, log_years     := log(YEARS)]
  dat[, dbh2          := DBH1^2]
  dat[, bal_over_ba   := fifelse(BA1 > 0, BAL1 / BA1, 0)]
  dat[, sqrt_ba_rd    := sqrt(BA1 * rd_add)]

  ## Factor indices
  sp_levels <- sort(unique(dat$SPCD))
  L1_levels <- sort(unique(dat$EPA_L1_CODE))
  L2_levels <- sort(unique(dat$EPA_L2_CODE))
  L3_levels <- sort(unique(dat$EPA_L3_CODE))

  dat[, sp_idx := match(SPCD, sp_levels)]
  dat[, L1_idx := match(EPA_L1_CODE, L1_levels)]
  dat[, L2_idx := match(EPA_L2_CODE, L2_levels)]
  dat[, L3_idx := match(EPA_L3_CODE, L3_levels)]

  mort_rate <- 1 - mean(dat$alive)
  message("  Obs:        ", format(nrow(dat), big.mark = ","))
  message("  Mortality:  ", round(100 * mort_rate, 2), "%")
  message("  Species:    ", length(sp_levels),
          "  (pooled rare = ", length(rare), ")")
  message("  EPA L1/L2/L3: ",
          length(L1_levels), " / ",
          length(L2_levels), " / ",
          length(L3_levels))

  list(data = dat,
       sp = sp_levels, L1 = L1_levels, L2 = L2_levels, L3 = L3_levels,
       cspi_shift = cspi_shift, site_var = site_var)
}

## ---- 2. Trait design matrix ------------------------------------------------
build_trait_matrix <- function(prep, traits_path,
                               trait_cols = DEFAULT_TRAIT_COLS,
                               enable = TRUE, verbose = TRUE) {
  if (!enable) {
    if (verbose) message("Trait prior disabled.")
    return(list(W = matrix(0, nrow = length(prep$sp), ncol = 0),
                trait_cols = character(0), P_trait = 0L))
  }
  if (!file.exists(traits_path)) {
    message("Trait file not found at ", traits_path,
            ".  Falling back to non-trait RE.")
    return(list(W = matrix(0, nrow = length(prep$sp), ncol = 0),
                trait_cols = character(0), P_trait = 0L))
  }
  tr <- readRDS(traits_path)
  setDT(tr)
  keep <- intersect(trait_cols, names(tr))
  if (length(keep) == 0L) {
    message("No usable trait columns found; disabling trait prior.")
    return(list(W = matrix(0, nrow = length(prep$sp), ncol = 0),
                trait_cols = character(0), P_trait = 0L))
  }
  tr_sp <- tr[match(prep$sp, SPCD), c("SPCD", keep), with = FALSE]
  tr_sp[, SPCD := prep$sp]
  informative <- character(0)
  for (cc in keep) {
    v <- tr_sp[[cc]]
    med <- median(v, na.rm = TRUE)
    if (is.na(med)) med <- 0
    v[is.na(v)] <- med
    s <- sd(v)
    if (!is.finite(s) || s < 1e-12) {
      if (verbose) message("  dropping zero-variance trait: ", cc)
      tr_sp[, (cc) := NULL]
      next
    }
    tr_sp[, (cc) := (v - mean(v)) / s]
    informative <- c(informative, cc)
  }
  if (length(informative) == 0L) {
    message("No informative trait columns after imputation; disabling trait prior.")
    return(list(W = matrix(0, nrow = length(prep$sp), ncol = 0),
                trait_cols = character(0), P_trait = 0L))
  }
  W <- as.matrix(tr_sp[, informative, with = FALSE])
  if (verbose) {
    message("Trait design W:  ", nrow(W), " species x ", ncol(W), " traits")
    message("  traits kept: ", paste(informative, collapse = ", "))
  }
  list(W = W, trait_cols = informative, P_trait = ncol(W))
}

## ---- 3. Stan data packer ---------------------------------------------------
make_stan_data_mort <- function(prep, tmat) {
  d <- prep$data
  W <- tmat$W
  if (ncol(W) == 0L) W <- matrix(0, nrow = length(prep$sp), ncol = 1)
  list(
    N_obs         = nrow(d),
    N_sp          = length(prep$sp),
    N_L1          = length(prep$L1),
    N_L2          = length(prep$L2),
    N_L3          = length(prep$L3),
    P_trait       = tmat$P_trait,
    alive         = d$alive,
    log_years     = d$log_years,
    dbh           = d$DBH1,
    dbh2          = d$dbh2,
    bal_over_ba   = d$bal_over_ba,
    cr            = d$CR1,
    sqrt_ba_rd    = d$sqrt_ba_rd,
    ln_cspi_shift = d$ln_cspi_shift,
    sp_idx        = d$sp_idx,
    L1_idx        = d$L1_idx,
    L2_idx        = d$L2_idx,
    L3_idx        = d$L3_idx,
    W             = W
  )
}

## ---- 4. Driver -------------------------------------------------------------
if (!interactive() && Sys.getenv("FVS_CONUS_SKIP_DRIVER") != "1") {

  opts <- OptionParser(option_list = list(
    make_option("--data",       type = "character",
                default = "calibration/data/conus_remeasurement_pairs.rds"),
    make_option("--stan_file",  type = "character",
                default = "scripts/stan/mort_gompit.stan"),
    make_option("--site",       type = "character", default = "cspi"),
    make_option("--traits_rds", type = "character",
                default = "calibration/traits/species_traits.rds"),
    make_option("--traits",     type = "integer",   default = 1L,
                help = "1 = enable trait-informed species prior, 0 = disable"),
    make_option("--out",        type = "character",
                default = "calibration/output/conus/mortality"),
    make_option("--n",          type = "integer",   default = 0L),
    make_option("--chains",     type = "integer",   default = 4L),
    make_option("--warmup",     type = "integer",   default = 1000L),
    make_option("--sampling",   type = "integer",   default = 1000L),
    make_option("--seed",       type = "integer",   default = 42L),
    make_option("--min_sp_n",   type = "integer",   default = 5000L),
    make_option("--save_draws", action = "store_true", default = FALSE)
  )) |> parse_args()

  stopifnot(file.exists(opts$stan_file))
  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

  cat("==============================================================\n")
  cat("FVS-CONUS Mortality (gompit) fit (driver)\n")
  cat("Stan file:    ", opts$stan_file, "\n")
  cat("Site variable:", opts$site, "\n")
  cat("Trait prior:  ", if (opts$traits == 1) "enabled" else "disabled", "\n")
  cat("Data:         ", opts$data, "\n")
  cat("Chains:       ", opts$chains, "  warmup:", opts$warmup,
      "  sampling:", opts$sampling, "\n")
  cat("==============================================================\n\n")

  dat_in <- readRDS(opts$data)
  setDT(dat_in)
  if (opts$n > 0 && opts$n < nrow(dat_in)) {
    set.seed(opts$seed)
    dat_in <- dat_in[sample(.N, opts$n)]
    message("Sub-sampled to ", format(nrow(dat_in), big.mark = ","), " rows")
  }

  prep  <- prepare_mort_data(dat_in, site_var = opts$site,
                             min_sp_n = opts$min_sp_n)
  tmat  <- build_trait_matrix(prep, opts$traits_rds,
                              enable = (opts$traits == 1L))
  sdata <- make_stan_data_mort(prep, tmat)

  mod <- cmdstan_model(opts$stan_file)

  t0  <- Sys.time()
  fit <- mod$sample(
    data            = sdata,
    chains          = opts$chains,
    parallel_chains = opts$chains,
    iter_warmup     = opts$warmup,
    iter_sampling   = opts$sampling,
    seed            = opts$seed,
    adapt_delta     = 0.95,      # tighter for cloglog
    max_treedepth   = 13,
    refresh         = 200
  )
  t1 <- Sys.time()
  cat("\nWall time:", round(as.numeric(t1 - t0, units = "mins"), 1), "min\n")

  fixef_pat <- "^(m[0-9]+|gamma\\[)"
  summ_fx   <- fit$summary() |>
    filter(grepl(fixef_pat, variable) |
           variable %in% c("sigma_sp","sigma_L1","sigma_L2","sigma_L3"))

  tag <- sprintf("mort_gompit_%s_traits%d", opts$site, opts$traits)
  write_csv(summ_fx, file.path(opts$out, sprintf("%s_summary.csv", tag)))
  saveRDS(list(
    site_var   = opts$site,
    traits     = opts$traits == 1L,
    trait_cols = tmat$trait_cols,
    stan_file  = opts$stan_file,
    prep_meta  = list(sp = prep$sp, L1 = prep$L1, L2 = prep$L2, L3 = prep$L3,
                      cspi_shift = prep$cspi_shift),
    summary    = summ_fx,
    n_obs      = nrow(prep$data),
    wall_min   = as.numeric(t1 - t0, units = "mins")
  ), file.path(opts$out, sprintf("%s_meta.rds", tag)))

  if (opts$save_draws) {
    fit$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))
  }

  cat("\nDone. Outputs under:", opts$out, "\n")
}
