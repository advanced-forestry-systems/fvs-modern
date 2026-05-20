## =============================================================================
## 36_fit_ht_dbh.R
##
## FVS-CONUS: unified driver for the static height-diameter (H-D) equation.
## Three candidate forms are fit in parallel (selectable via --form):
##
##   chapman  -> scripts/stan/ht_dbh_chapman.stan  (Chapman-Richards on A)
##   wykoff   -> scripts/stan/ht_dbh_wykoff.stan   (Wykoff 1986 log-HT form)
##   schnute  -> scripts/stan/ht_dbh_schnute.stan  (Schnute 2-parameter)
##
## Shared covariate block (Aaron, 2026-04-15): BAL, sqrt_BA, ln(CSPI_shift),
## BA x RD, BAL x RD. All forms use nested EPA L1 / L2 / L3 + species
## random intercepts on the growth-scaling parameter, with an optional
## trait-informed species prior:
##
##      z_sp ~ Normal( W . gamma, sigma_sp )
##
## where W is a species x trait design matrix built from
## calibration/traits/species_traits.rds (script 42). Set --traits 0 to
## disable the trait prior (reduces to a plain hierarchical RE).
##
## Usage:
##   Rscript scripts/36_fit_ht_dbh.R --form chapman
##   Rscript scripts/36_fit_ht_dbh.R --form wykoff  --site cspi
##   Rscript scripts/36_fit_ht_dbh.R --form schnute --traits 0
##
## Author: A. Weiskittel, 2026-04-15
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(optparse)
})

KNOWN_FORMS <- c("chapman", "wykoff", "schnute")

## Traits kept by default (must exist in species_traits.rds). Additional
## columns are welcome; just add them here and they will be centered and
## scaled before entering W. Missing values are median-imputed per column.
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
prepare_htdbh_data <- function(dat, site_var = "cspi", min_sp_n = 2000L) {
  message("Preparing HT-DBH data. site_var = ", site_var,
          "  min_sp_n = ", min_sp_n)
  setDT(dat)

  ## Use the time-1 HT / DBH observations; each remeasurement row
  ## contributes one cross-sectional observation. (We can optionally
  ## include time-2 later if sample is thin.)
  dat <- dat[
    STATUS1 == 1 &
    !is.na(HT1) & HT1 >= 1.37 & HT1 < 85 &
    !is.na(DBH1) & DBH1 >= 1.0 & DBH1 < 250 &
    !is.na(BA1)  & BA1  > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    !is.na(get(site_var)) & get(site_var) > 0
  ]

  ## Rare species pooling keeps the species RE stable without flooding it
  ## with singletons.
  sp_n <- dat[, .N, by = SPCD]
  rare <- sp_n[N < min_sp_n, SPCD]
  dat[SPCD %in% rare, SPCD := -1L]

  cspi_shift <- max(1.0 - min(dat[[site_var]], na.rm = TRUE), 0.01)
  dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]
  dat[, sqrt_ba       := sqrt(BA1)]
  dat[, ba_x_rd       := BA1 * rd_add]
  dat[, bal_x_rd      := BAL1 * rd_add]

  sp_levels <- sort(unique(dat$SPCD))
  L1_levels <- sort(unique(dat$EPA_L1_CODE))
  L2_levels <- sort(unique(dat$EPA_L2_CODE))
  L3_levels <- sort(unique(dat$EPA_L3_CODE))

  dat[, sp_idx := match(SPCD, sp_levels)]
  dat[, L1_idx := match(EPA_L1_CODE, L1_levels)]
  dat[, L2_idx := match(EPA_L2_CODE, L2_levels)]
  dat[, L3_idx := match(EPA_L3_CODE, L3_levels)]

  message("  Obs:        ", format(nrow(dat), big.mark = ","))
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
##
## Returns a list:
##   $W          numeric matrix [N_sp, P_trait] in order of prep$sp
##   $trait_cols kept columns after NA-screening
##   $P_trait    ncol(W)
##
## If enable = FALSE (or traits file missing, or no usable columns), returns
## an empty matrix and P_trait = 0 so the Stan code falls back to the plain
## non-centered RE parameterization.
build_trait_matrix <- function(prep, traits_path, trait_cols = DEFAULT_TRAIT_COLS,
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

  ## Align to modeled species order
  tr_sp <- tr[match(prep$sp, SPCD),
              c("SPCD", keep), with = FALSE]
  tr_sp[, SPCD := prep$sp]

  ## Per-column median impute, then center and scale.
  ## Drop zero-variance columns (all NAs or constant across modeled species)
  ## so they do not create unidentifiable gamma parameters.
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

## ---- 3. Stan data packers --------------------------------------------------
shared_htdbh_pack <- function(prep, tmat) {
  d <- prep$data
  W <- tmat$W
  if (ncol(W) == 0L) W <- matrix(0, nrow = length(prep$sp), ncol = 1)  # dummy for Stan
  list(
    N_obs         = nrow(d),
    N_sp          = length(prep$sp),
    N_L1          = length(prep$L1),
    N_L2          = length(prep$L2),
    N_L3          = length(prep$L3),
    P_trait       = tmat$P_trait,
    ht_obs        = d$HT1,
    dbh           = d$DBH1,
    bal           = d$BAL1,
    sqrt_ba       = d$sqrt_ba,
    ln_cspi_shift = d$ln_cspi_shift,
    ba_x_rd       = d$ba_x_rd,
    bal_x_rd      = d$bal_x_rd,
    sp_idx        = d$sp_idx,
    L1_idx        = d$L1_idx,
    L2_idx        = d$L2_idx,
    L3_idx        = d$L3_idx,
    W             = W
  )
}

make_stan_data_chapman <- function(prep, tmat) shared_htdbh_pack(prep, tmat)
make_stan_data_wykoff  <- function(prep, tmat) shared_htdbh_pack(prep, tmat)

make_stan_data_schnute <- function(prep, tmat, dbh_ref = 25.0) {
  base <- shared_htdbh_pack(prep, tmat)
  base$dbh_ref <- dbh_ref
  base
}

stan_data_packer_ht <- function(form) {
  switch(form,
         chapman = make_stan_data_chapman,
         wykoff  = make_stan_data_wykoff,
         schnute = make_stan_data_schnute)
}

## ---- 4. Driver -------------------------------------------------------------
if (!interactive() && Sys.getenv("FVS_CONUS_SKIP_DRIVER") != "1") {

  opts <- OptionParser(option_list = list(
    make_option("--form",       type = "character", default = "chapman",
                help = paste("H-D form:", paste(KNOWN_FORMS, collapse = " | "))),
    make_option("--data",       type = "character",
                default = "calibration/data/conus_remeasurement_pairs.rds"),
    make_option("--stan_dir",   type = "character", default = "scripts/stan"),
    make_option("--site",       type = "character", default = "cspi"),
    make_option("--traits_rds", type = "character",
                default = "calibration/traits/species_traits.rds"),
    make_option("--traits",     type = "integer",   default = 1L,
                help = "1 = enable trait-informed species prior, 0 = disable"),
    make_option("--out",        type = "character",
                default = "calibration/output/conus/ht_dbh"),
    make_option("--n",          type = "integer",   default = 0L),
    make_option("--chains",     type = "integer",   default = 4L),
    make_option("--warmup",     type = "integer",   default = 1000L),
    make_option("--sampling",   type = "integer",   default = 1000L),
    make_option("--seed",       type = "integer",   default = 42L),
    make_option("--min_sp_n",   type = "integer",   default = 2000L),
    make_option("--dbh_ref",    type = "double",    default = 25.0,
                help = "Schnute anchor DBH (cm); ignored for other forms"),
    make_option("--save_draws", action = "store_true", default = FALSE)
  )) |> parse_args()

  stopifnot(opts$form %in% KNOWN_FORMS)

  stan_file <- file.path(opts$stan_dir,
                         switch(opts$form,
                                chapman = "ht_dbh_chapman.stan",
                                wykoff  = "ht_dbh_wykoff.stan",
                                schnute = "ht_dbh_schnute.stan"))
  stopifnot(file.exists(stan_file))

  ## Schnute has the flexible extra shape parameter and needs slightly
  ## tighter adaptation; the log-linear Wykoff is cheapest.
  form_adapt_delta   <- if (opts$form == "schnute") 0.95 else 0.9
  form_max_treedepth <- if (opts$form == "schnute") 13   else 12

  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

  cat("==============================================================\n")
  cat("FVS-CONUS HT-DBH Multi-Model Fit (driver)\n")
  cat("Form:         ", opts$form, "\n")
  cat("Stan file:    ", stan_file, "\n")
  cat("Site variable:", opts$site, "\n")
  cat("Trait prior:  ", if (opts$traits == 1) "enabled" else "disabled", "\n")
  cat("Data:         ", opts$data, "\n")
  cat("==============================================================\n\n")

  dat_in <- readRDS(opts$data)
  setDT(dat_in)
  if (opts$n > 0 && opts$n < nrow(dat_in)) {
    set.seed(opts$seed)
    dat_in <- dat_in[sample(.N, opts$n)]
    message("Sub-sampled to ", format(nrow(dat_in), big.mark = ","), " rows")
  }

  prep <- prepare_htdbh_data(dat_in, site_var = opts$site,
                             min_sp_n = opts$min_sp_n)
  tmat <- build_trait_matrix(prep, opts$traits_rds,
                             enable = (opts$traits == 1L))
  packer <- stan_data_packer_ht(opts$form)
  sdata  <- if (opts$form == "schnute") packer(prep, tmat, dbh_ref = opts$dbh_ref)
            else                         packer(prep, tmat)

  mod <- cmdstan_model(stan_file)

  t0 <- Sys.time()
  fit <- mod$sample(
    data            = sdata,
    chains          = opts$chains,
    parallel_chains = opts$chains,
    iter_warmup     = opts$warmup,
    iter_sampling   = opts$sampling,
    seed            = opts$seed,
    adapt_delta     = form_adapt_delta,
    max_treedepth   = form_max_treedepth,
    refresh         = 200
  )
  t1 <- Sys.time()
  cat("\nWall time:", round(as.numeric(t1 - t0, units = "mins"), 1), "min\n")

  fixef_pat <- "^(a0|a_bal|a_ba|a_cspi|a_bard|a_blrd|b0|b1|b_rate|c_shape|a_rate|b_shape|s0|s1|gamma\\[|trait_effect\\[)"
  ## Targeted summary: skip slow full fit$summary().
  all_vars  <- fit$metadata()$model_params
  keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                        all_vars %in% c("sigma_sp","sigma_L1",
                                        "sigma_L2","sigma_L3")]
  summ_fx   <- fit$summary(variables = keep_vars)

  tag <- sprintf("htdbh_%s_%s_traits%d", opts$form, opts$site, opts$traits)
  write_csv(summ_fx, file.path(opts$out, sprintf("%s_summary.csv", tag)))
  saveRDS(list(
    form        = opts$form,
    site_var    = opts$site,
    traits      = opts$traits == 1L,
    trait_cols  = tmat$trait_cols,
    stan_file   = stan_file,
    prep_meta   = list(sp = prep$sp, L1 = prep$L1, L2 = prep$L2, L3 = prep$L3,
                       cspi_shift = prep$cspi_shift),
    summary     = summ_fx,
    n_obs       = nrow(prep$data),
    wall_min    = as.numeric(t1 - t0, units = "mins")
  ), file.path(opts$out, sprintf("%s_meta.rds", tag)))

  if (opts$save_draws) {
    fit$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))
  }

  cat("\nDone. Outputs under:", opts$out, "\n")
}
