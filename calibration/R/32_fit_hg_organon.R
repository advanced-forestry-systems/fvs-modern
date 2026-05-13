## =============================================================================
## 32_fit_hg_organon.R
##
## FVS-CONUS: Production driver for the ORGANON-analogue height growth (HG)
## model, log-normal likelihood with fixed K1h / K2h / K4h shape constants.
##
## This driver mirrors 31_fit_dg_multimodel.R (trait-aware species prior,
## nested EPA L1/L2/L3 + species crossed random intercepts, Cao/Weiskittel
## variance scaling, same CLI surface) but is specialised for the HG Stan
## contract in calibration/stan/hg_organon_fixedK.stan (v4).
##
## Linear predictor (log annual height growth, log m/yr):
##   ln(dHT_a) = a0
##             + trait_effect[sp] + z_sp[sp]
##             + z_L1[l1] + z_L2[l2] + z_L3[l3]
##             + a1 ln(HT + K1h) + a2 HT^K2h
##             + a3 ln((CR + 0.2)/1.2)
##             + a4 ln(CSPI_shift)
##             + a5 CCFL / ln(HT + K4h)
##             + a6 sqrt(BA) + a7 BA*RD + a8 CCFL*RD
##
## Likelihood: log_hg_obs_a ~ Normal(eta, sigma / sqrt(YEARS))
##
## Usage:
##   Rscript calibration/R/32_fit_hg_organon.R \
##     --data calibration/data/conus_remeasurement_pairs_metric.rds \
##     --stan_file calibration/stan/hg_organon_fixedK.stan \
##     --site cspi --traits_rds calibration/traits/species_traits.rds \
##     --traits 1 --out calibration/output/conus/hg \
##     --n 0 --chains 2 --warmup 500 --sampling 1000 --seed 42 \
##     --min_sp_n 5000 --save_draws --fix_K --K1h 1.0 --K2h 0.5 --K4h 2.7
##
## Author: A. Weiskittel, 2026-04-22 (restored after 3e00a81 reset)
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(optparse)
})

## ---- Configuration ---------------------------------------------------------

## Default trait columns kept for the species-level prior. Missing columns
## in species_traits.rds are silently dropped; zero-variance ones are
## dropped inside build_trait_matrix().
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
prepare_hg_data <- function(dat, site_var = "cspi", min_sp_n = 5000L) {
  message("Preparing HG data. site_var = ", site_var,
          "  min_sp_n = ", min_sp_n)
  setDT(dat)

  ## Observed annual HT growth (m/yr)
  dat[, hg_obs_a := (HT2 - HT1) / YEARS]

  ## Basic plausibility filters. Require strictly positive increment for
  ## the log-normal likelihood.
  dat <- dat[
    !is.na(HT1) & HT1 >= 1.37 &
    !is.na(HT2) & HT2 > HT1 &
    hg_obs_a > 0.001 & hg_obs_a < 2.0 &                  # m/yr
    !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(CCFL1) & CCFL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    YEARS >= 1 & YEARS <= 20 &
    !is.na(get(site_var)) & get(site_var) > 0
  ]

  ## Rare species pooling into pseudo-SPCD = -1L
  sp_n <- dat[, .N, by = SPCD]
  rare <- sp_n[N < min_sp_n, SPCD]
  dat[SPCD %in% rare, SPCD := -1L]

  ## Shifted log of site productivity so values <= 0 do not break log()
  cspi_shift <- 1.0 - min(dat[[site_var]], na.rm = TRUE)
  cspi_shift <- max(cspi_shift, 0.01)
  dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]

  ## Derived covariates. ba_x_rd and ccfl_x_rd interactions match the
  ## Stan contract in hg_organon_fixedK.stan v4.
  dat[, log_hg_obs_a := log(hg_obs_a)]
  dat[, ln_cr_adj    := log((CR1 + 0.2) / 1.2)]
  dat[, sqrt_ba      := sqrt(BA1)]
  dat[, ba_x_rd      := BA1   * rd_add]
  dat[, ccfl_x_rd    := CCFL1 * rd_add]
  dat[, sqrt_years   := sqrt(YEARS)]

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
    v   <- tr_sp[[cc]]
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
make_stan_data_hg <- function(prep, tmat, K1h, K2h, K4h) {
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

    log_hg_obs_a  = d$log_hg_obs_a,
    sqrt_years    = d$sqrt_years,

    ht            = d$HT1,
    ln_cr_adj     = d$ln_cr_adj,
    ln_cspi_shift = d$ln_cspi_shift,
    ccfl          = d$CCFL1,
    sqrt_ba       = d$sqrt_ba,
    ba_x_rd       = d$ba_x_rd,
    ccfl_x_rd     = d$ccfl_x_rd,

    K1h           = K1h,
    K2h           = K2h,
    K4h           = K4h,

    sp_idx        = d$sp_idx,
    L1_idx        = d$L1_idx,
    L2_idx        = d$L2_idx,
    L3_idx        = d$L3_idx,

    W             = W
  )
}

## ---- 4. Diagnostics --------------------------------------------------------
check_convergence <- function(fit) {
  summ <- fit$summary(variables = c("sigma", "sigma_sp",
                                    "sigma_L1", "sigma_L2", "sigma_L3"))
  summ$rhat_ok <- summ$rhat < 1.05
  summ$ess_ok  <- summ$ess_bulk > 100
  summ
}

## ---- 5. Driver -------------------------------------------------------------
if (!interactive() && Sys.getenv("FVS_CONUS_SKIP_DRIVER") != "1") {

  opts <- OptionParser(option_list = list(
    make_option("--data",       type = "character",
                default = "calibration/data/conus_remeasurement_pairs_metric.rds"),
    make_option("--stan_file",  type = "character",
                default = "calibration/stan/hg_organon_fixedK.stan"),
    make_option("--site",       type = "character", default = "cspi"),
    make_option("--traits_rds", type = "character",
                default = "calibration/traits/species_traits.rds"),
    make_option("--traits",     type = "integer",   default = 1L,
                help = "1 = enable trait-informed species prior, 0 = disable"),
    make_option("--out",        type = "character",
                default = "calibration/output/conus/hg"),
    make_option("--n",          type = "integer",   default = 0L,
                help = "Optional subsample size; 0 = use all"),
    make_option("--chains",     type = "integer",   default = 2L),
    make_option("--warmup",     type = "integer",   default = 500L),
    make_option("--sampling",   type = "integer",   default = 1000L),
    make_option("--seed",       type = "integer",   default = 42L),
    make_option("--min_sp_n",   type = "integer",   default = 5000L),
    make_option("--save_draws", action = "store_true", default = FALSE),
    make_option("--fix_K",      action = "store_true", default = FALSE,
                help = "Flag required by fixedK Stan file (accepted for parity)"),
    make_option("--K1h",        type = "double",    default = 1.0),
    make_option("--K2h",        type = "double",    default = 0.5),
    make_option("--K4h",        type = "double",    default = 2.7),
    make_option("--adapt_delta", type = "double",   default = 0.9),
    make_option("--max_treedepth", type = "integer", default = 10L),
    make_option("--init_scale",   type = "double",    default = 0.1)
  )) |> parse_args()

  stopifnot(file.exists(opts$stan_file))
  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

  cat("==============================================================\n")
  cat("FVS-CONUS Height Growth Production Fit\n")
  cat("Stan file:     ", opts$stan_file, "\n")
  cat("Site variable: ", opts$site, "\n")
  cat("Trait prior:   ", if (opts$traits == 1L) "enabled" else "disabled", "\n")
  cat("K1h / K2h / K4h:", opts$K1h, "/", opts$K2h, "/", opts$K4h, "\n")
  cat("Data:          ", opts$data, "\n")
  cat("Chains:        ", opts$chains,
      " warmup:", opts$warmup,
      " sampling:", opts$sampling, "\n")
  cat("Output:        ", opts$out, "\n")
  cat("==============================================================\n\n")

  dat_in <- readRDS(opts$data)
  setDT(dat_in)

  if (opts$n > 0L && opts$n < nrow(dat_in)) {
    set.seed(opts$seed)
    dat_in <- dat_in[sample(.N, opts$n)]
    message("Sub-sampled to ", format(nrow(dat_in), big.mark = ","), " rows")
  }

  prep   <- prepare_hg_data(dat_in, site_var = opts$site,
                            min_sp_n = opts$min_sp_n)
  tmat   <- build_trait_matrix(prep, opts$traits_rds,
                               enable = (opts$traits == 1L))
  sdata  <- make_stan_data_hg(prep, tmat,
                              K1h = opts$K1h,
                              K2h = opts$K2h,
                              K4h = opts$K4h)

  mod <- cmdstan_model(opts$stan_file)

  t0 <- Sys.time()
  # Prior-centered init list per chain with small jitter.
  # This fixes 'Location parameter[1] is inf' at init caused by random
  # uniform(-scale, scale) draws landing in high-gradient regions.
  jitter_scale <- opts$init_scale
  jitter <- function(mu, sd) rnorm(1, mu, sd * jitter_scale)
  N_sp  <- sdata$N_sp
  N_L1  <- sdata$N_L1
  N_L2  <- sdata$N_L2
  N_L3  <- sdata$N_L3
  P_trait <- sdata$P_trait
  init_fn <- function(chain_id) {
    set.seed(opts$seed + chain_id)
    list(
      a0       = jitter(-2.5,   0.3),
      a1       = jitter( 0.5,   0.1),
      a2       = jitter(-0.03,  0.01),
      a3       = jitter( 0.8,   0.1),
      a4       = jitter( 0.4,   0.1),
      a5       = jitter(-0.005, 0.002),
      a6       = jitter(-0.03,  0.01),
      a7       = jitter( 0.0,   0.005),
      a8       = jitter( 0.0,   0.005),
      gamma    = rnorm(P_trait, 0, 0.1 * jitter_scale),
      z_sp_raw = rnorm(N_sp, 0, 0.1 * jitter_scale),
      z_L1_raw = rnorm(N_L1, 0, 0.1 * jitter_scale),
      z_L2_raw = rnorm(N_L2, 0, 0.1 * jitter_scale),
      z_L3_raw = rnorm(N_L3, 0, 0.1 * jitter_scale),
      sigma_sp = 0.1,
      sigma_L1 = 0.15,
      sigma_L2 = 0.1,
      sigma_L3 = 0.1,
      sigma    = 0.3
    )
  }
  inits <- lapply(seq_len(opts$chains), init_fn)

  fit <- mod$sample(
    data            = sdata,
    chains          = opts$chains,
    parallel_chains = opts$chains,
    iter_warmup     = opts$warmup,
    iter_sampling   = opts$sampling,
    seed            = opts$seed,
    adapt_delta     = opts$adapt_delta,
    max_treedepth   = opts$max_treedepth,
    init            = inits,
    refresh         = 100
  )
  t1 <- Sys.time()
  cat("\nWall time:", round(as.numeric(t1 - t0, units = "mins"), 1), "min\n")

  fixef_pat <- "^(a[0-9]+|gamma\\[|trait_effect\\[|sigma)"
  ## Targeted summary: skip the 30+ min full fit$summary() over thousands
  ## of random-effect entries on a 5+ GB fit.
  all_vars  <- fit$metadata()$model_params
  keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                        all_vars %in% c("sigma_sp","sigma_L1",
                                        "sigma_L2","sigma_L3")]
  summ_fx <- fit$summary(variables = keep_vars)

  tag <- sprintf("hg_organon_fixedK_%s_traits%d", opts$site, opts$traits)
  write_csv(summ_fx, file.path(opts$out, sprintf("%s_summary.csv", tag)))

  saveRDS(list(
    stan_file  = opts$stan_file,
    site_var   = opts$site,
    traits     = opts$traits == 1L,
    trait_cols = tmat$trait_cols,
    K1h        = opts$K1h, K2h = opts$K2h, K4h = opts$K4h,
    prep_meta  = list(sp = prep$sp, L1 = prep$L1, L2 = prep$L2, L3 = prep$L3,
                      cspi_shift = prep$cspi_shift),
    summary    = summ_fx,
    n_obs      = nrow(prep$data),
    wall_min   = as.numeric(t1 - t0, units = "mins")
  ), file.path(opts$out, sprintf("%s_meta.rds", tag)))

  if (opts$save_draws) {
    fit$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))
  }

  conv <- check_convergence(fit)
  cat("\nConvergence (variance components):\n")
  print(conv)

  cat("\nDone. Outputs under:", opts$out, "\n")
}
