##=============================================================================
## 33_fit_hcb_organon.R
## ORGANON form Height to Crown Base Model for CONUS wide FVS
##
## Fits the Hanus, Hann, and Marshall (2000) HCB_SWO logistic model with
## species x ecodivision Bayesian hierarchy via CmdStan.
##
## HCB is modeled as a static cross sectional equation fit to all FIA trees
## with measured crown ratio. Crown ratio is then derived as CR = 1 - HCB/HT.
##
## For projection, crown recession is computed as the difference:
##   HCBG = HCB_pred(t2 attributes) - HCB_pred(t1 attributes)
## with the constraint HCBG >= 0 (crown base does not move downward).
##
## ALL UNITS METRIC: HT in m, DBH in cm, BA in m2/ha
##
## Usage:
##   Rscript 33_fit_hcb_organon.R
##
## Author: Aaron Weiskittel, Greg Johnson
## Date: 2026-04-11
##=============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

## Source project utilities (skipped on Cardinal; other production drivers
## do not source 00_setup.R either. The fitted environment is expected to
## have R module + cmdstanr + required packages already installed).
# if (file.exists("calibration/R/00_setup.R")) source("calibration/R/00_setup.R")

##-----------------------------------------------------------------------------
## 1. Configuration
##-----------------------------------------------------------------------------

MIN_OBS_SPECIES <- 5000
STAN_MODEL <- "calibration/stan/organon_hcb_conus.stan"
OUT_DIR <- "calibration/output/conus"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

##-----------------------------------------------------------------------------
## 2. Data Preparation
##-----------------------------------------------------------------------------

prepare_hcb_data <- function(dat, site_var = "climate_si") {

  message("Preparing HCB data...")
  message("  Site productivity variable: ", site_var)

  ## Convert FIA imperial to metric
  dat <- dat %>%
    mutate(
      DBH1_cm    = DBH1 * 2.54,
      HT1_m      = HT1 * 0.3048,
      BA1_metric = BA1 * 0.2296,          # ft2/ac to m2/ha
      BAL1_metric = BAL1 * 0.2296,
      ## Derive HCB from crown ratio: HCB = HT * (1 - CR)
      HCB1_m     = HT1_m * (1 - CR1),
      ## Derived predictors
      ln_ba      = log(pmax(BA1_metric, 0.1)),
      dbh_ht_ratio = DBH1_cm / pmax(HT1_m, 0.1),
      site_prod  = pmax(.data[[site_var]], 0.1)
    )

  ## Compute Curtis relative density
  ## QMD derived from BA and TPA if not directly available
  dat <- dat %>%
    mutate(
      QMD1_metric = if ("QMD1" %in% names(.)) QMD1 * 2.54
                    else sqrt(BA1_metric / (TPA1 * 2.47105) * 40000 / pi),
      rd = BA1_metric / pmax(sqrt(QMD1_metric), 0.01)
    )

  ## Filter: valid observations
  dat <- dat %>%
    filter(
      DBH1_cm >= 2.54,                   # >= 1 inch
      HT1_m > 1.37,                      # above breast height
      CR1 > 0 & CR1 <= 1.0,
      HCB1_m >= 0 & HCB1_m < HT1_m,
      BA1_metric > 0,
      site_prod > 0,
      !is.na(CCFL1)
    )

  ## Species filtering
  sp_counts <- dat %>% count(SPCD) %>% filter(n >= MIN_OBS_SPECIES)
  dat <- dat %>% filter(SPCD %in% sp_counts$SPCD)

  ## Create integer indices
  sp_levels  <- sort(unique(dat$SPCD))
  eco_levels <- sort(unique(dat$ecodiv_code))

  dat <- dat %>%
    mutate(
      species_idx = match(SPCD, sp_levels),
      ecodiv_idx  = match(ecodiv_code, eco_levels)
    )

  message("  Trees: ", nrow(dat))
  message("  Species: ", length(sp_levels))
  message("  Ecodivisions: ", length(eco_levels))
  message("  Mean HCB: ", round(mean(dat$HCB1_m), 2), " m")
  message("  Mean CR: ", round(mean(dat$CR1), 3))

  list(
    data    = dat,
    species = sp_levels,
    ecodiv  = eco_levels
  )
}

##-----------------------------------------------------------------------------
## 3. Stan Data Assembly
##-----------------------------------------------------------------------------

make_stan_data_hcb <- function(prep) {
  dat <- prep$data

  list(
    N          = nrow(dat),
    N_species  = length(prep$species),
    N_ecodiv   = length(prep$ecodiv),
    hcb_obs    = dat$HCB1_m,
    ht         = dat$HT1_m,
    ccfl       = dat$CCFL1,
    ln_ba      = dat$ln_ba,
    dbh_ht_ratio = dat$dbh_ht_ratio,
    site_prod  = dat$site_prod,
    rd         = dat$rd,
    species_id = dat$species_idx,
    ecodiv_id  = dat$ecodiv_idx
  )
}

##-----------------------------------------------------------------------------
## 4. Model Fitting
##-----------------------------------------------------------------------------

fit_hcb <- function(prep, stan_file = STAN_MODEL,
                    chains = 4, iter_warmup = 1000,
                    iter_sampling = 1000, seed = 42) {

  message("Compiling Stan model: ", stan_file)
  mod <- cmdstan_model(stan_file)

  stan_data <- make_stan_data_hcb(prep)

  message("Running HMC sampling...")
  fit <- mod$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = min(chains, parallel::detectCores()),
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    seed            = seed,
    max_treedepth   = 10,
    adapt_delta     = 0.9,
    refresh         = 100
  )

  return(fit)
}

##-----------------------------------------------------------------------------
## 5. Convergence Diagnostics
##-----------------------------------------------------------------------------

check_convergence <- function(fit) {

  params <- c("mu_b0", "b1", "b2", "b3", "b4", "b5", "b6",
              "sigma_sp", "sigma_eco", "sigma")

  message("Checking convergence...")
  summ <- fit$summary(variables = params)
  print(summ)

  rhat_bad <- summ %>% filter(rhat > 1.05)
  if (nrow(rhat_bad) > 0) {
    warning("Parameters with Rhat > 1.05:")
    print(rhat_bad)
  } else {
    message("  All Rhat < 1.05: PASS")
  }

  ess_bad <- summ %>% filter(ess_bulk < 400)
  if (nrow(ess_bad) > 0) {
    warning("Parameters with ESS_bulk < 400:")
    print(ess_bad)
  } else {
    message("  All ESS_bulk > 400: PASS")
  }

  diag <- fit$diagnostic_summary()
  n_div <- sum(diag$num_divergent)
  message("  Divergent transitions: ", n_div)

  return(summ)
}

##-----------------------------------------------------------------------------
## 6. Model Evaluation
##-----------------------------------------------------------------------------

evaluate_hcb <- function(fit, prep) {

  message("Evaluating HCB model...")
  dat <- prep$data

  ## Extract posterior medians
  summ <- fit$summary()
  get_med <- function(var) summ %>% filter(variable == var) %>% pull(median)

  sp_int  <- fit$draws("b0_sp", format = "draws_matrix") %>% apply(2, median)
  eco_int <- fit$draws("b0_eco", format = "draws_matrix") %>% apply(2, median)

  ## Compute predicted HCB for each observation
  logit_cr <- get_med("mu_b0") +
    sp_int[dat$species_idx] +
    eco_int[dat$ecodiv_idx] +
    get_med("b1") * dat$HT1_m +
    get_med("b2") * dat$CCFL1 +
    get_med("b3") * dat$ln_ba +
    get_med("b4") * dat$dbh_ht_ratio +
    get_med("b5") * dat$site_prod +
    get_med("b6") * dat$rd

  hcb_pred <- dat$HT1_m / (1 + exp(logit_cr))
  cr_pred  <- 1 - hcb_pred / dat$HT1_m

  ## Residual statistics
  resid_hcb <- dat$HCB1_m - hcb_pred
  resid_cr  <- dat$CR1 - cr_pred

  stats <- tibble(
    metric    = c("HCB", "CR"),
    n         = c(length(resid_hcb), length(resid_cr)),
    bias      = c(mean(resid_hcb), mean(resid_cr)),
    rmse      = c(sqrt(mean(resid_hcb^2)), sqrt(mean(resid_cr^2))),
    pseudo_r2 = c(
      1 - sum(resid_hcb^2) / sum((dat$HCB1_m - mean(dat$HCB1_m))^2),
      1 - sum(resid_cr^2) / sum((dat$CR1 - mean(dat$CR1))^2)
    )
  )

  message("  Evaluation results:")
  print(stats)

  ## Species level comparison
  by_species <- dat %>%
    mutate(hcb_pred = hcb_pred, cr_pred = cr_pred) %>%
    group_by(SPCD) %>%
    summarise(
      n          = n(),
      obs_cr     = mean(CR1),
      pred_cr    = mean(cr_pred),
      obs_hcb    = mean(HCB1_m),
      pred_hcb   = mean(hcb_pred),
      .groups    = "drop"
    )

  message("  Species level CR comparison (top 20):")
  print(by_species %>% arrange(desc(n)) %>% head(20))

  list(
    stats      = stats,
    by_species = by_species,
    hcb_pred   = hcb_pred,
    cr_pred    = cr_pred
  )
}

##-----------------------------------------------------------------------------
## 7. Dynamic Crown Ratio Prediction
##-----------------------------------------------------------------------------

#' Predict crown ratio change over a projection period using the static
#' HCB model evaluated at beginning and end of period attributes.
#'
#' CR at time t is derived as: CR_t = 1 - HCB_pred(t) / HT(t)
#' Crown recession: HCBG = max(HCB_pred(t2) - HCB_pred(t1), 0)
#' Updated CR: CR2 = 1 - (HCB1 + HCBG) / HT2
#'
#' @param params Named list with mu_b0, b1:b6, sp_int, eco_int vectors
#' @param dat Data frame with t1 and t2 tree/stand attributes
#' @return Data frame with predicted HCB2, CR2, and HCBG

predict_dynamic_cr <- function(params, dat) {

  ## HCB at time 1 (from current attributes)
  logit_cr1 <- params$mu_b0 +
    params$sp_int[dat$species_idx] +
    params$eco_int[dat$ecodiv_idx] +
    params$b1 * dat$HT1_m +
    params$b2 * dat$CCFL1 +
    params$b3 * log(pmax(dat$BA1_metric, 0.1)) +
    params$b4 * (dat$DBH1_cm / pmax(dat$HT1_m, 0.1)) +
    params$b5 * dat$site_prod +
    params$b6 * dat$rd1

  hcb1_pred <- dat$HT1_m / (1 + exp(logit_cr1))

  ## HCB at time 2 (from projected attributes)
  logit_cr2 <- params$mu_b0 +
    params$sp_int[dat$species_idx] +
    params$eco_int[dat$ecodiv_idx] +
    params$b1 * dat$HT2_m +
    params$b2 * dat$CCFL2 +
    params$b3 * log(pmax(dat$BA2_metric, 0.1)) +
    params$b4 * (dat$DBH2_cm / pmax(dat$HT2_m, 0.1)) +
    params$b5 * dat$site_prod +
    params$b6 * dat$rd2

  hcb2_pred <- dat$HT2_m / (1 + exp(logit_cr2))

  ## Crown recession: constrained non negative (ORGANON convention)
  hcbg <- pmax(hcb2_pred - hcb1_pred, 0)

  ## Updated HCB and CR
  hcb2_final <- hcb1_pred + hcbg
  cr2_pred   <- 1 - hcb2_final / pmax(dat$HT2_m, 0.1)
  cr2_pred   <- pmax(pmin(cr2_pred, 1.0), 0.01)  # bound to (0.01, 1.0)

  tibble(
    hcb1_pred  = hcb1_pred,
    hcb2_pred  = hcb2_final,
    hcbg       = hcbg,
    cr2_pred   = cr2_pred
  )
}

##-----------------------------------------------------------------------------
## 8. Posterior Export
##-----------------------------------------------------------------------------

export_posteriors <- function(fit, prep, out_dir = OUT_DIR) {

  ## Fixed effects
  draws <- fit$draws(
    variables = c("mu_b0", "b1", "b2", "b3", "b4", "b5", "b6",
                  "sigma_sp", "sigma_eco", "sigma"),
    format = "draws_df"
  )
  write_csv(draws, file.path(out_dir, "hcb_organon_fixed_draws.csv"))

  ## Species intercepts
  sp_draws <- fit$draws("b0_sp", format = "draws_matrix")
  sp_medians <- apply(sp_draws, 2, median)
  sp_summary <- tibble(
    SPCD  = prep$species,
    b0_sp = sp_medians[seq_along(prep$species)]
  )
  write_csv(sp_summary, file.path(out_dir, "hcb_organon_species_intercepts.csv"))

  ## Ecodivision intercepts
  eco_draws <- fit$draws("b0_eco", format = "draws_matrix")
  eco_medians <- apply(eco_draws, 2, median)
  eco_summary <- tibble(
    ecodiv = prep$ecodiv,
    b0_eco = eco_medians[seq_along(prep$ecodiv)]
  )
  write_csv(eco_summary, file.path(out_dir, "hcb_organon_ecodiv_intercepts.csv"))

  message("Posteriors exported to: ", out_dir)
  list(fixed = draws, species = sp_summary, ecodiv = eco_summary)
}

##-----------------------------------------------------------------------------
## 9. Main Execution
##-----------------------------------------------------------------------------

if (sys.nframe() == 0) {

  args <- commandArgs(trailingOnly = TRUE)
  site_var <- if ("--bgi" %in% args) "bgi" else "climate_si"

  message("=== ORGANON HCB CONUS wide Model Fitting ===")
  message("Site variable: ", site_var)

  ## Load data
  data_file <- "calibration/data/conus_remeasurement_pairs.rds"
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file,
         "\nRun 30_build_conus_dataset.R first.")
  }
  raw_data <- readRDS(data_file)

  ## Prepare, fit, evaluate, export
  prep <- prepare_hcb_data(raw_data, site_var = site_var)
  fit  <- fit_hcb(prep)
  summ <- check_convergence(fit)
  eval <- evaluate_hcb(fit, prep)
  posts <- export_posteriors(fit, prep)

  message("=== Done ===")
}
