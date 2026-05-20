## =============================================================================
## 31_fit_dg_multimodel.R
##
## FVS-CONUS: Unified driver for the four candidate diameter growth model
## forms (ORGANON, Weiskittel 2016 SS_REAL, Kuehne 2022 FVS_ACD, Greg
## Johnson compound ratio). All four share the same nested EPA L1 / L2 /
## L3 + species random intercept hierarchy and the same annual growth
## likelihood with Cao (2000) variance scaling so LOO and held out plot
## MAB differences attach to the mean function alone.
##
## Usage:
##   Rscript scripts/31_fit_dg_multimodel.R --form weiskittel
##   Rscript scripts/31_fit_dg_multimodel.R --form kuehne    --site cspi
##   Rscript scripts/31_fit_dg_multimodel.R --form compound  --n 100000
##   Rscript scripts/31_fit_dg_multimodel.R --form organon   --chains 4
##
## Forms accepted: organon | weiskittel | kuehne | compound
##
## Upstream assumption:
##   conus_remeasurement_pairs.rds has already been converted to metric
##   in script 30 Step 5b. DBH in cm, HT in m, BA / BAL in m^2/ha, RD is
##   additive SDI / SDImax_brms, BAL_SW1 and BAL_HW1 populated, EPA_L1_CODE
##   and EPA_L2_CODE parsed from NA_L3CODE.
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

## ---- Configuration ---------------------------------------------------------

KNOWN_FORMS <- c("organon", "weiskittel", "kuehne", "compound", "organon_gj")

## Traits kept by default for the species-level prior. Columns absent from
## species_traits.rds are dropped; remaining columns are median-imputed,
## centered, and scaled before being multiplied by gamma inside Stan.
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

## Note: option parsing and driver setup live inside the driver guard at the
## bottom so this file can be safely sourced for its helper functions without
## consuming caller script argv.

## ---- 1. Data prep ----------------------------------------------------------
##
## Shared preparation. Units are assumed metric (script 30 Step 5b).
## Produces a list with $data (data.table), $sp (factor levels),
## $L1, $L2, $L3 (factor levels).
## ---------------------------------------------------------------------------
prepare_dg_data <- function(dat, site_var = "cspi", min_sp_n = 5000L) {
  message("Preparing DG data. site_var = ", site_var,
          "  min_sp_n = ", min_sp_n)
  setDT(dat)

  ## Observed annual DBH growth
  dat[, dg_obs_a := (DBH2 - DBH1) / YEARS]

  ## Basic plausibility filters (tighter filters live in 39_qc)
  dat <- dat[
    STATUS1 == 1 & STATUS2 == 1 &
    !is.na(DBH1) & DBH1 > 0 &
    !is.na(DBH2) & DBH2 > 0 &
    dg_obs_a > -0.5 & dg_obs_a < 3.0 &            # cm/yr
    !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(BAL1) & !is.na(BAL_SW1) & !is.na(BAL_HW1) &
    !is.na(HT1) & HT1 >= 1.37 &                    # compound form needs HT
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    YEARS >= 1 & YEARS <= 20 &
    !is.na(get(site_var)) & get(site_var) > 0
  ]

  ## Rare species pooling. Below threshold pool into a generic "rare"
  ## group so the species RE does not over-shrink.
  sp_n <- dat[, .N, by = SPCD]
  rare <- sp_n[N < min_sp_n, SPCD]
  dat[SPCD %in% rare, SPCD := -1L]

  ## Derived covariates used across the four forms
  cspi_shift <- 1.0 - min(dat[[site_var]], na.rm = TRUE)
  cspi_shift <- max(cspi_shift, 0.01)
  dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]
  dat[, ln_dbh        := log(DBH1)]
  dat[, ln_cr_adj     := log((CR1 + 0.2) / 1.2)]
  dat[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
  dat[, sqrt_ba       := sqrt(BA1)]
  dat[, sqrt_years    := sqrt(YEARS)]

  ## Factor indices
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
  message("  EPA L1:     ", length(L1_levels))
  message("  EPA L2:     ", length(L2_levels))
  message("  EPA L3:     ", length(L3_levels))

  list(data = dat,
       sp = sp_levels, L1 = L1_levels, L2 = L2_levels, L3 = L3_levels,
       cspi_shift = cspi_shift, site_var = site_var)
}

## ---- 2. Trait design matrix ------------------------------------------------
##
## Returns list(W, trait_cols, P_trait). If enable = FALSE, traits file
## missing, or no usable columns, returns an empty 0-column matrix and
## P_trait = 0 so the Stan code falls back to the plain non-centered RE.
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

  ## Per-column median impute, then center and scale.
  ## Track which columns survive: any column that ends up zero-variance

  ## (all NAs imputed to the same value, or constant across modeled species)
  ## is dropped so it does not create an unidentifiable gamma parameter.
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

## ---- 3. Stan data packers ---------------------------------------------------
## One packer per form; each returns a named list consumed by the
## corresponding Stan file. tmat carries the species x trait design;
## when tmat$P_trait == 0, a dummy 1-column zero matrix is sent so Stan's
## matrix[N_sp, P_trait > 0 ? P_trait : 1] W stays well defined.
## ---------------------------------------------------------------------------
shared_idx <- function(d, prep, tmat) {
  W <- tmat$W
  if (ncol(W) == 0L) W <- matrix(0, nrow = length(prep$sp), ncol = 1)
  list(
    N_obs       = nrow(d),
    N_sp        = length(prep$sp),
    N_L1        = length(prep$L1),
    N_L2        = length(prep$L2),
    N_L3        = length(prep$L3),
    P_trait     = tmat$P_trait,
    dg_obs_a    = d$dg_obs_a,
    sqrt_years  = d$sqrt_years,
    sp_idx      = d$sp_idx,
    L1_idx      = d$L1_idx,
    L2_idx      = d$L2_idx,
    L3_idx      = d$L3_idx,
    W           = W
  )
}

make_stan_data_weiskittel <- function(prep, tmat) {
  d <- prep$data
  c(shared_idx(d, prep, tmat), list(
    ln_dbh        = d$ln_dbh,
    dbh           = d$DBH1,
    ln_cr_adj     = d$ln_cr_adj,
    bal_sw        = d$BAL_SW1,
    bal_hw        = d$BAL_HW1,
    rd            = d$rd_add,
    ln_cspi_shift = d$ln_cspi_shift,
    ba_x_rd       = d$ba_x_rd,
    bal_x_rd      = d$bal_x_rd
  ))
}

make_stan_data_kuehne <- function(prep, tmat) {
  d <- prep$data
  c(shared_idx(d, prep, tmat), list(
    ln_dbh        = d$ln_dbh,
    dbh           = d$DBH1,
    ln_cr_adj     = d$ln_cr_adj,
    ln_bal_sw_adj = d$ln_bal_sw_adj,
    bal_hw        = d$BAL_HW1,
    ln_csi        = d$ln_cspi_shift,
    ba_x_rd       = d$ba_x_rd,
    bal_x_rd      = d$bal_x_rd
  ))
}

make_stan_data_compound <- function(prep, tmat) {
  d <- prep$data
  c(shared_idx(d, prep, tmat), list(
    dbh           = d$DBH1,
    ht            = pmax(d$HT1, 1.37),
    cr            = d$CR1,
    bal_sw        = d$BAL_SW1,
    bal_hw        = d$BAL_HW1,
    rd            = d$rd_add,
    ln_cspi_shift = d$ln_cspi_shift,
    ba_x_rd       = d$ba_x_rd,
    bal_x_rd      = d$bal_x_rd
  ))
}

make_stan_data_organon <- function(prep, tmat) {
  d <- prep$data
  c(shared_idx(d, prep, tmat), list(
    dbh           = d$DBH1,
    ln_cr_adj     = d$ln_cr_adj,
    ln_cspi_shift = d$ln_cspi_shift,
    bal           = d$BAL1,
    sqrt_ba       = d$sqrt_ba,
    ba_x_rd       = d$ba_x_rd,
    bal_x_rd      = d$bal_x_rd
  ))
}

make_stan_data_organon_gj <- function(prep, tmat) {
  ## Greg Johnson simplified ORGANON. No CSPI, no BA/BAL x RD interactions.
  d <- prep$data
  c(shared_idx(d, prep, tmat), list(
    dbh       = d$DBH1,
    ln_cr_adj = d$ln_cr_adj,
    bal       = d$BAL1,
    sqrt_ba   = d$sqrt_ba
  ))
}

stan_data_packer <- function(form) {
  switch(form,
         organon    = make_stan_data_organon,
         weiskittel = make_stan_data_weiskittel,
         kuehne     = make_stan_data_kuehne,
         compound   = make_stan_data_compound,
         organon_gj = make_stan_data_organon_gj)
}

## ---- 3. Convergence / diagnostic helpers -----------------------------------
check_convergence <- function(fit) {
  summ <- fit$summary(variables = c("sigma", "sigma_sp",
                                    "sigma_L1", "sigma_L2", "sigma_L3"))
  summ$rhat_ok  <- summ$rhat < 1.05
  summ$ess_ok   <- summ$ess_bulk > 100
  summ
}

## ---- 4. Driver --------------------------------------------------------------
## Guard so that pre-flight and other callers can source this file purely to
## grab prepare_dg_data() and stan_data_packer() without triggering a full
## production fit. Sentinel env var lets sourcing scripts opt out.
if (!interactive() && Sys.getenv("FVS_CONUS_SKIP_DRIVER") != "1") {

  opts <- OptionParser(option_list = list(
    make_option("--form",       type = "character", default = "weiskittel",
                help = paste("Model form:",
                             paste(KNOWN_FORMS, collapse = " | "))),
    make_option("--data",       type = "character",
                default = "calibration/data/conus_remeasurement_pairs.rds"),
    make_option("--stan_dir",   type = "character",
                default = "scripts/stan"),
    make_option("--stan_file",  type = "character", default = "",
                help = "Optional override; if set, used instead of --form lookup"),
    make_option("--site",       type = "character", default = "cspi"),
    make_option("--traits_rds", type = "character",
                default = "calibration/traits/species_traits.rds"),
    make_option("--traits",     type = "integer",   default = 1L,
                help = "1 = enable trait-informed species prior, 0 = disable"),
    make_option("--out",        type = "character",
                default = "calibration/output/conus/dg"),
    make_option("--n",          type = "integer",   default = 0L,
                help = "Optional subsample size; 0 = use all"),
    make_option("--chains",     type = "integer",   default = 4L),
    make_option("--warmup",     type = "integer",   default = 1000L),
    make_option("--sampling",   type = "integer",   default = 1000L),
    make_option("--seed",       type = "integer",   default = 42L),
    make_option("--min_sp_n",   type = "integer",   default = 5000L),
    make_option("--save_draws", action = "store_true", default = FALSE),
    make_option("--max_treedepth", type = "integer", default = 10L),
    make_option("--adapt_delta",   type = "double",  default = 0.9)
  )) |> parse_args()

  stopifnot(opts$form %in% KNOWN_FORMS)

  stan_file <- switch(opts$form,
    organon    = file.path(opts$stan_dir, "dg_organon.stan"),
    weiskittel = file.path(opts$stan_dir, "dg_weiskittel2016.stan"),
    kuehne     = file.path(opts$stan_dir, "dg_kuehne2022.stan"),
    compound   = file.path(opts$stan_dir, "dg_compound_gj.stan"),
    organon_gj = file.path(opts$stan_dir, "dg_organon_gj.stan"))
  if (nzchar(opts$stan_file)) stan_file <- opts$stan_file
  stopifnot(file.exists(stan_file))

  ## Form-specific sampler tuning. After the 2026-04-16 linearization of
  ## the compound form (k1/k2 -> c_num/c_den), all forms use the same
  ## defaults. Keep the conditional structure for future per-form tuning.
  form_adapt_delta   <- opts$adapt_delta
  form_max_treedepth <- opts$max_treedepth

  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

  cat("==============================================================\n")
  cat("FVS-CONUS DG Multi-Model Fit (driver)\n")
  cat("Form:         ", opts$form, "\n")
  cat("Stan file:    ", stan_file, "\n")
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

  prep   <- prepare_dg_data(dat_in, site_var = opts$site,
                            min_sp_n = opts$min_sp_n)
  tmat   <- build_trait_matrix(prep, opts$traits_rds,
                               enable = (opts$traits == 1L))
  packer <- stan_data_packer(opts$form)
  sdata  <- packer(prep, tmat)

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

  ## Summary on fixed effects only (plus variance components)
  ## Targeted: skip full fit$summary() (30+ min stall on 5+ GB fits).
  fixef_pat <- "^(a[0-9]+|b[0-9]+|c[0-9]+|g[0-9]+|s[0-9]+|K[0-9]?h?|k[12]|gamma\\[|trait_effect\\[)"
  all_vars  <- fit$metadata()$model_params
  keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                        all_vars %in% c("sigma","sigma_sp","sigma_L1",
                                        "sigma_L2","sigma_L3")]
  summ_fx   <- fit$summary(variables = keep_vars)

  tag <- sprintf("dg_%s_%s_traits%d", opts$form, opts$site, opts$traits)
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
