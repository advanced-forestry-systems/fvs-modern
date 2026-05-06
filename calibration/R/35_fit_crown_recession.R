## =============================================================================
## 35_fit_crown_recession.R
##
## FVS-CONUS: Crown recession (ΔHCB) driver for cr_recession_hh.stan
## Following the incremental method of Hann & Hanus (2004, Can. J. For. Res.
## 34: 1993-2003), adapted for CONUS-wide application with nested EPA
## L1/L2/L3 + species random intercepts and optional trait-informed prior.
##
## The response is the ratio r = ΔHCB / (CL_S + ΔH) on (0, 1), where:
##   ΔHCB  = HCB2 - HCB1 (crown recession during the interval)
##   CL_S  = HT1 - HCB1  (crown length at start)
##   ΔH    = HT2 - HT1   (height growth during the interval)
##   CL_S + ΔH = HT2 - HCB1 (maximum possible crown recession)
##
## This ratio is modeled with a Beta likelihood using the Hann & Hanus
## logistic denominator form: r = 1/(1+exp(eta)) = inv_logit(-eta).
## The resulting ΔHCB is always bounded [0, CL_S + ΔH].
##
## During FVS simulation, crown ratio is updated as:
##   HCB2 = HCB1 + ΔHCB
##   CR2  = 1 - HCB2 / HT2
##
## Data requirements: trees with CR1, CR2, HT1, HT2 all populated so that
## HCB can be derived as HT * (1 - CR). More restrictive than the static
## HCB model (script 33) because CR2 is also required.
##
## Inputs:
##   calibration/data/conus_remeasurement_pairs.rds (metric, per 30 Step 5b)
##   calibration/traits/species_traits.rds          (optional, for trait prior)
##
## Usage:
##   Rscript scripts/35_fit_crown_recession.R --site cspi
##   Rscript scripts/35_fit_crown_recession.R --site bgi --traits 0
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

## Shared trait columns (must match 31 / 32 / 33 / 34 / 36 for consistency)
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
prepare_cr_data <- function(dat, site_var = "cspi", min_sp_n = 5000L) {
  message("Preparing crown recession data. site_var = ", site_var,
          "  min_sp_n = ", min_sp_n)
  setDT(dat)

  ## Derive HCB from HT and CR (HCB not stored directly in FIA data)
  dat[, HCB1 := HT1 * (1 - CR1)]
  dat[, HCB2 := HT2 * (1 - CR2)]

  ## Require both HCB and HT valid at start AND end
  dat <- dat[
    STATUS1 == 1 & STATUS2 == 1 &
    !is.na(HT1)  & HT1  >= 1.37 &
    !is.na(HT2)  & HT2  >= 1.37 &
    !is.na(CR1)  & CR1 > 0 & CR1 < 1.0 &
    !is.na(CR2)  & CR2 > 0 & CR2 < 1.0 &
    !is.na(HCB1) & HCB1 > 0 &
    !is.na(HCB2) & HCB2 > 0 &
    !is.na(DBH1) & DBH1 >= 1.0 &
    !is.na(BA1)  & BA1  > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    !is.na(get(site_var)) & get(site_var) > 0 &
    YEARS >= 1 & YEARS <= 20
  ]

  ## Compute ΔHCB and the bounding quantities
  dat[, delta_hcb := HCB2 - HCB1]
  dat[, cl_s      := HT1 - HCB1]        # crown length at start
  dat[, delta_h   := HT2 - HT1]         # height growth
  dat[, max_dhcb  := cl_s + delta_h]     # = HT2 - HCB1, maximum possible

  ## Filter: crown recession must be non-negative and within bounds
  ## Allow small negative ΔHCB (measurement error) but clamp to near zero
  ## Also require positive max_dhcb (otherwise ratio is undefined)
  dat <- dat[max_dhcb > 0.1]             # at least 10 cm of potential recession
  dat[delta_hcb < 0, delta_hcb := 0]     # clamp epicormic / measurement error

  ## Ratio r = ΔHCB / (CL_S + ΔH), clamped to (0.001, 0.999) for Beta
  dat[, cr_ratio := delta_hcb / max_dhcb]
  dat[, cr_ratio := pmin(pmax(cr_ratio, 0.001), 0.999)]

  ## Drop trees with no recession at all (r ~ 0) if they dominate
  ## Actually keep them; the Beta handles near-zero values with phi control.
  ## But do filter height shrinkage (HT2 < HT1 by more than measurement error)
  dat <- dat[delta_h >= -0.5]            # allow small measurement error in HT

  ## Rare species pooling
  sp_n <- dat[, .N, by = SPCD]
  rare <- sp_n[N < min_sp_n, SPCD]
  dat[SPCD %in% rare, SPCD := -1L]

  ## Derived covariates matching cr_recession_hh.stan
  cspi_shift <- max(1.0 - min(dat[[site_var]], na.rm = TRUE), 0.01)
  dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]
  dat[, ln_cr         := log(CR1)]
  dat[, sqrt_ba       := sqrt(BA1)]
  dat[, ln_bal_ba     := log(BAL1 / BA1 + 1)]
  dat[, cr_over_rd    := CR1 / rd_add]

  ## Factor indices
  sp_levels <- sort(unique(dat$SPCD))
  L1_levels <- sort(unique(dat$EPA_L1_CODE))
  L2_levels <- sort(unique(dat$EPA_L2_CODE))
  L3_levels <- sort(unique(dat$EPA_L3_CODE))

  dat[, sp_idx := match(SPCD, sp_levels)]
  dat[, L1_idx := match(EPA_L1_CODE, L1_levels)]
  dat[, L2_idx := match(EPA_L2_CODE, L2_levels)]
  dat[, L3_idx := match(EPA_L3_CODE, L3_levels)]

  message("  Obs:          ", format(nrow(dat), big.mark = ","))
  message("  Species:      ", length(sp_levels),
          "  (pooled rare = ", length(rare), ")")
  message("  EPA L1/L2/L3: ",
          length(L1_levels), " / ",
          length(L2_levels), " / ",
          length(L3_levels))
  message("  ΔHCB range:   [", round(min(dat$delta_hcb), 3), ", ",
          round(max(dat$delta_hcb), 3), "] m")
  message("  ratio range:  [", round(min(dat$cr_ratio), 4), ", ",
          round(max(dat$cr_ratio), 4), "]")
  message("  Mean ΔHCB:    ", round(mean(dat$delta_hcb), 3), " m")
  message("  Median CR1:   ", round(median(dat$CR1), 3))
  message("  Mean YEARS:   ", round(mean(dat$YEARS), 1))

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
make_stan_data_cr <- function(prep, tmat) {
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
    ratio         = d$cr_ratio,
    ln_cr         = d$ln_cr,
    cr            = d$CR1,
    sqrt_ba       = d$sqrt_ba,
    ln_bal_ba     = d$ln_bal_ba,
    cr_over_rd    = d$cr_over_rd,
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
                default = "scripts/stan/cr_recession_hh.stan"),
    make_option("--site",       type = "character", default = "cspi"),
    make_option("--traits_rds", type = "character",
                default = "calibration/traits/species_traits.rds"),
    make_option("--traits",     type = "integer",   default = 1L,
                help = "1 = enable trait-informed species prior, 0 = disable"),
    make_option("--out",        type = "character",
                default = "calibration/output/conus/crown_recession"),
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
  cat("FVS-CONUS Crown Recession (Hann & Hanus 2004) fit (driver)\n")
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

  prep  <- prepare_cr_data(dat_in, site_var = opts$site,
                           min_sp_n = opts$min_sp_n)
  tmat  <- build_trait_matrix(prep, opts$traits_rds,
                              enable = (opts$traits == 1L))
  sdata <- make_stan_data_cr(prep, tmat)

  mod <- cmdstan_model(opts$stan_file)

  t0  <- Sys.time()
  fit <- mod$sample(
    data            = sdata,
    chains          = opts$chains,
    parallel_chains = opts$chains,
    iter_warmup     = opts$warmup,
    iter_sampling   = opts$sampling,
    seed            = opts$seed,
    adapt_delta     = 0.9,
    max_treedepth   = 10,
    refresh         = 200
  )
  t1 <- Sys.time()
  cat("\nWall time:", round(as.numeric(t1 - t0, units = "mins"), 1), "min\n")

  fixef_pat <- "^(r[0-9]+|gamma\\[|trait_effect\\[)"
  ## Targeted summary: ask cmdstan for only the variables we report.
  all_vars  <- fit$metadata()$model_params
  keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                        all_vars %in% c("phi","sigma_sp","sigma_L1",
                                        "sigma_L2","sigma_L3")]
  summ_fx   <- fit$summary(variables = keep_vars)

  print(summ_fx |> select(variable, mean, median, sd, rhat, ess_bulk), n = 50)

  tag <- sprintf("cr_recession_%s_traits%d", opts$site, opts$traits)
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
