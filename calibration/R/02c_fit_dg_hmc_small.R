#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Wykoff DG Model (HMC on reduced data, v3)
# For variants where ADVI diverged (sigma >> 1).
# Strategy: MAP -> HMC on 5K obs -> constrained ADVI fallback
#
# Key changes from 02/02b:
#   1. Reduces max observations to 5000 for tractable HMC
#   2. Tries HMC first (feasible at 5K obs)
#   3. Constrained ADVI fallback with informative sigma prior
#   4. Uses wykoff_dg_v3.stan with sigma ~ normal(0.8, 0.3) prior
#
# Usage: Rscript 02c_fit_dg_hmc_small.R --variant cr
#

library(tidyverse)
library(data.table)
library(cmdstanr)
library(jsonlite)
library(posterior)
library(bayesplot)
library(tidybayes)
library(parallel)
library(progress)
library(logger)

# ============================================================================
# Parse Command Line Arguments
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

variant <- "ca"
species_codes <- NULL

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variant <- args[i + 1]
    }
    if (args[i] == "--species" & i < length(args)) {
      species_codes <- as.integer(strsplit(args[i + 1], ",")[[1]])
    }
  }
}

# ============================================================================
# Configuration
# ============================================================================

project_root <- Sys.getenv("FVS_PROJECT_ROOT", unset = "")
if (project_root == "") {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) getwd()
  )
}
calibration_dir <- file.path(project_root, "calibration")
processed_data_dir <- file.path(calibration_dir, "data", "processed")
output_dir <- file.path(calibration_dir, "output", "variants", variant)
stan_dir <- file.path(calibration_dir, "stan")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MAX_OBS <- as.integer(Sys.getenv("FVS_MAX_OBS", "5000"))

# HMC sampler config (env-tunable for convergence improvement)
HMC_WARMUP       <- as.integer(Sys.getenv("FVS_HMC_WARMUP",       "500"))
HMC_SAMPLING     <- as.integer(Sys.getenv("FVS_HMC_SAMPLING",     "500"))
HMC_ADAPT_DELTA  <- as.numeric(Sys.getenv("FVS_HMC_ADAPT_DELTA",  "0.90"))
HMC_MAX_TREEDEPTH<- as.integer(Sys.getenv("FVS_HMC_TREEDEPTH",    "10"))
HMC_CHAINS       <- as.integer(Sys.getenv("FVS_HMC_CHAINS",       "4"))
logger::log_info("HMC config: chains={HMC_CHAINS} warmup={HMC_WARMUP} sampling={HMC_SAMPLING} adapt_delta={HMC_ADAPT_DELTA} max_treedepth={HMC_MAX_TREEDEPTH}")

log_file <- file.path(calibration_dir, "logs", paste0("02c_fit_dg_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting DG HMC-small fitting for variant {variant}")
logger::log_info("Max observations: {MAX_OBS}")

# ============================================================================
# Load Data (same as 02b)
# ============================================================================

data_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")
if (!file.exists(data_file)) stop("Data file not found: ", data_file)

dg_data <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble()
logger::log_info("Loaded {nrow(dg_data)} observations")

if (is.null(species_codes)) {
  species_codes <- unique(dg_data$SPCD)
}

dg_data <- dg_data %>% filter(SPCD %in% species_codes)

# ============================================================================
# Load FVS Config for Priors
# ============================================================================

config_file <- file.path(project_root, "config", paste0(variant, ".json"))
config <- fromJSON(config_file)
fvs_params <- config$categories$other

param_map <- list(
  b1_ref = "DGLD", b2_ref = "DGDS", b3_ref = "DGSITE", b4_ref = "DGSLOP",
  b5_ref = "DGSLSQ", b6_ref = "DGSASP", b7_ref = "DGCASP", b8_ref = "DGEL",
  b9_ref = "DGELSQ", b10_ref = "DGCR", b11_ref = "DGCRSQ", b12_ref = "DGBAL",
  b13_ref = "DGPCCF"
)

prior_means <- sapply(param_map, function(param_name) {
  vals <- fvs_params[[param_name]]
  if (is.null(vals) || length(vals) == 0) return(0)
  mean(unlist(vals), na.rm = TRUE)
})

# ============================================================================
# Prepare Data
# ============================================================================

dg_model <- dg_data %>%
  mutate(species_idx = as.numeric(factor(SPCD)))

dg_model <- dg_model %>%
  mutate(
    ln_DBH_std = as.numeric(scale(ln_DBH)),
    SLOPE_std = as.numeric(scale(SLOPE)),
    ELEV_std = as.numeric(scale(ELEV)),
    DBH_sq_std = as.numeric(scale(DBH_sq)),
    ln_SI_std = as.numeric(scale(ln_SI)),
    SLOPE_sq = SLOPE_std^2,
    ELEV_sq = ELEV_std^2,
    SLOPE_SASP_std = as.numeric(scale(SLOPE_SASP)),
    SLOPE_CASP_std = as.numeric(scale(SLOPE_CASP)),
    CR = CR_pct / 100,
    CR_sq = CR^2,
    BAL_std = as.numeric(scale(BAL)),
    BA_std = as.numeric(scale(BA))
  )

dg_model <- dg_model %>%
  filter(
    is.finite(ln_DBH_std), is.finite(SLOPE_std), is.finite(ELEV_std),
    is.finite(ln_DDS), is.finite(ln_SI_std), is.finite(DBH_sq_std),
    is.finite(SLOPE_SASP_std), is.finite(SLOPE_CASP_std),
    !is.na(CR), is.finite(BAL_std), is.finite(BA_std)
  )

# Save standardization parameters
std_params <- list(
  ln_DBH_mean = mean(dg_model$ln_DBH, na.rm = TRUE),
  ln_DBH_sd = sd(dg_model$ln_DBH, na.rm = TRUE),
  SLOPE_mean = mean(dg_model$SLOPE, na.rm = TRUE),
  SLOPE_sd = sd(dg_model$SLOPE, na.rm = TRUE),
  ELEV_mean = mean(dg_model$ELEV, na.rm = TRUE),
  ELEV_sd = sd(dg_model$ELEV, na.rm = TRUE),
  DBH_sq_mean = mean(dg_model$DBH_sq, na.rm = TRUE),
  DBH_sq_sd = sd(dg_model$DBH_sq, na.rm = TRUE),
  ln_SI_mean = mean(dg_model$ln_SI, na.rm = TRUE),
  ln_SI_sd = sd(dg_model$ln_SI, na.rm = TRUE),
  SLOPE_SASP_mean = mean(dg_model$SLOPE_SASP, na.rm = TRUE),
  SLOPE_SASP_sd = sd(dg_model$SLOPE_SASP, na.rm = TRUE),
  SLOPE_CASP_mean = mean(dg_model$SLOPE_CASP, na.rm = TRUE),
  SLOPE_CASP_sd = sd(dg_model$SLOPE_CASP, na.rm = TRUE),
  BAL_mean = mean(dg_model$BAL, na.rm = TRUE),
  BAL_sd = sd(dg_model$BAL, na.rm = TRUE),
  BA_mean = mean(dg_model$BA, na.rm = TRUE),
  BA_sd = sd(dg_model$BA, na.rm = TRUE)
)
saveRDS(std_params, file.path(output_dir, "standardization_params.rds"))

# Aggressive subsampling: 5K stratified by species
if (nrow(dg_model) > MAX_OBS) {
  logger::log_info("Subsampling from {nrow(dg_model)} to {MAX_OBS} (stratified)")
  set.seed(42)
  subsample_prop <- min(1, MAX_OBS / nrow(dg_model) * 1.1)
  dg_model <- dg_model %>%
    group_by(species_idx) %>%
    slice_sample(prop = subsample_prop) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), MAX_OBS))
  logger::log_info("After subsampling: {nrow(dg_model)} observations")
}

dg_model <- dg_model %>%
  mutate(species_idx = as.integer(factor(species_idx)))

N <- nrow(dg_model)
N_species <- max(dg_model$species_idx)
logger::log_info("N = {N}, N_species = {N_species}")

prior_b0_species <- rep(prior_means["b1_ref"] * 0.5, N_species)

stan_data <- list(
  N = N, N_species = N_species,
  ln_DDS = dg_model$ln_DDS,
  ln_DBH = dg_model$ln_DBH_std, DBH_sq = dg_model$DBH_sq_std,
  ln_SI = dg_model$ln_SI_std,
  SLOPE = dg_model$SLOPE_std, SLOPE_sq = dg_model$SLOPE_sq,
  SLOPE_SASP = dg_model$SLOPE_SASP_std, SLOPE_CASP = dg_model$SLOPE_CASP_std,
  ELEV = dg_model$ELEV_std, ELEV_sq = dg_model$ELEV_sq,
  CR = dg_model$CR, CR_sq = dg_model$CR_sq,
  BAL = dg_model$BAL_std, BA = dg_model$BA_std,
  species_id = dg_model$species_idx,
  prior_b0 = prior_b0_species
)

# ============================================================================
# Compile Stan Model v2 (with sigma lower bounds)
# ============================================================================

stan_file <- file.path(stan_dir, "wykoff_dg_v2.stan")
if (!file.exists(stan_file)) stan_file <- file.path(stan_dir, "wykoff_dg.stan")

mod <- cmdstanr::cmdstan_model(stan_file, dir = output_dir)

# ============================================================================
# OLS for initial values
# ============================================================================

ols_fit <- lm(ln_DDS ~ ln_DBH_std + DBH_sq_std + ln_SI_std +
                SLOPE_std + I(SLOPE_std^2) + SLOPE_SASP_std + SLOPE_CASP_std +
                ELEV_std + I(ELEV_std^2) + CR + I(CR^2) + BAL_std + BA_std,
              data = dg_model)
ols_coef <- coef(ols_fit)
ols_sigma <- sigma(ols_fit)

init_fn <- function() {
  list(
    mu_b0 = unname(ols_coef["(Intercept)"]),
    sigma_b0 = 0.5,
    z_b0 = rep(0, N_species),
    b1 = unname(ols_coef["ln_DBH_std"]),
    b2 = unname(ols_coef["DBH_sq_std"]),
    b3 = unname(ols_coef["ln_SI_std"]),
    b4 = unname(ols_coef["SLOPE_std"]),
    b5 = unname(ols_coef["I(SLOPE_std^2)"]),
    b6 = unname(ols_coef["SLOPE_SASP_std"]),
    b7 = unname(ols_coef["SLOPE_CASP_std"]),
    b8 = unname(ols_coef["ELEV_std"]),
    b9 = unname(ols_coef["I(ELEV_std^2)"]),
    b10 = unname(ols_coef["CR"]),
    b11 = unname(ols_coef["I(CR^2)"]),
    b12 = unname(ols_coef["BAL_std"]),
    b13 = unname(ols_coef["BA_std"]),
    sigma = ols_sigma
  )
}

logger::log_info("OLS sigma: {round(ols_sigma, 3)} (this is the target for Bayesian sigma)")

# ============================================================================
# Stage 1: MAP
# ============================================================================

logger::log_info("Stage 1: MAP optimization...")

fit_map <- tryCatch({
  mod$optimize(data = stan_data, init = init_fn, algorithm = "lbfgs",
               iter = 10000, tol_rel_grad = 1e-8, refresh = 100)
}, error = function(e) { logger::log_warn("MAP failed: {e$message}"); NULL })

if (!is.null(fit_map)) {
  map_estimates <- fit_map$summary()
  write_csv(map_estimates, file.path(output_dir, "diameter_growth_map.csv"))

  map_vals <- setNames(map_estimates$estimate, map_estimates$variable)
  init_from_map <- function() {
    z_b0_names <- paste0("z_b0[", 1:N_species, "]")
    list(
      mu_b0 = unname(map_vals["mu_b0"]),
      sigma_b0 = max(unname(map_vals["sigma_b0"]), 0.02),
      z_b0 = unname(map_vals[z_b0_names]),
      b1 = unname(map_vals["b1"]), b2 = unname(map_vals["b2"]),
      b3 = unname(map_vals["b3"]), b4 = unname(map_vals["b4"]),
      b5 = unname(map_vals["b5"]), b6 = unname(map_vals["b6"]),
      b7 = unname(map_vals["b7"]), b8 = unname(map_vals["b8"]),
      b9 = unname(map_vals["b9"]), b10 = unname(map_vals["b10"]),
      b11 = unname(map_vals["b11"]), b12 = unname(map_vals["b12"]),
      b13 = unname(map_vals["b13"]),
      sigma = max(min(unname(map_vals["sigma"]), 2.0), 0.1)
    )
  }
  map_sigma <- unname(map_vals["sigma"])
  logger::log_info("MAP sigma: {round(map_sigma, 3)}")
} else {
  init_from_map <- init_fn
  map_sigma <- ols_sigma
}

# ============================================================================
# Stage 2: HMC (feasible at 5K observations)
# ============================================================================

logger::log_info("Stage 2: HMC sampling on {N} observations...")

sampling_method <- "none"

fit <- tryCatch({
  result <- mod$sample(
    data = stan_data,
    chains = HMC_CHAINS,
    parallel_chains = min(HMC_CHAINS, parallel::detectCores()),
    iter_warmup = HMC_WARMUP,
    iter_sampling = HMC_SAMPLING,
    adapt_delta = HMC_ADAPT_DELTA,
    max_treedepth = HMC_MAX_TREEDEPTH,
    refresh = 50,
    show_messages = FALSE,
    init = init_from_map
  )

  diag <- result$diagnostic_summary()
  n_divergent <- sum(diag$num_divergent)
  n_treedepth <- sum(diag$num_max_treedepth)
  total_transitions <- 4 * 500

  # Check sigma is reasonable (not diverged)
  sigma_draws <- result$draws("sigma", format = "matrix")
  sigma_median <- median(sigma_draws)
  logger::log_info("HMC sigma median: {round(sigma_median, 3)}")

  if (n_divergent < total_transitions * 0.10 &&
      n_treedepth < total_transitions * 0.50 &&
      sigma_median < 3.0) {
    sampling_method <<- "hmc"
    logger::log_info("HMC converged: {n_divergent} divergent, sigma={round(sigma_median,3)}")
    result
  } else {
    logger::log_warn("HMC poor: {n_divergent} divergent, sigma={round(sigma_median,3)}")
    NULL
  }
}, error = function(e) {
  logger::log_warn("HMC failed: {e$message}")
  NULL
})

# ============================================================================
# Stage 3: Constrained ADVI fallback
# ============================================================================

if (is.null(fit) || sampling_method == "none") {
  logger::log_info("Stage 3: Constrained ADVI...")

  # Use fullrank with eta tuning and OLS sigma as reference
  advi_strategies <- list(
    list(algorithm = "fullrank", iter = 50000, tol = 0.005, eta = 0.05, label = "fullrank_eta005"),
    list(algorithm = "fullrank", iter = 50000, tol = 0.01,  eta = 0.01, label = "fullrank_eta001"),
    list(algorithm = "meanfield", iter = 50000, tol = 0.005, eta = 0.05, label = "meanfield_eta005"),
    list(algorithm = "meanfield", iter = 50000, tol = 0.01,  eta = 0.01, label = "meanfield_eta001")
  )

  for (strat in advi_strategies) {
    logger::log_info("Trying ADVI: {strat$label}...")

    fit_vb <- tryCatch({
      result <- mod$variational(
        data = stan_data,
        init = init_from_map,
        algorithm = strat$algorithm,
        iter = strat$iter,
        tol_rel_obj = strat$tol,
        eta = strat$eta,
        adapt_engaged = FALSE,
        output_samples = 4000,
        refresh = 1000
      )

      # Check if sigma is reasonable
      sigma_draws <- result$draws("sigma", format = "matrix")
      sigma_median <- median(sigma_draws)
      logger::log_info("ADVI {strat$label} sigma median: {round(sigma_median, 3)}")

      if (sigma_median < 3.0) {
        result
      } else {
        logger::log_warn("ADVI {strat$label} sigma too high: {round(sigma_median, 3)}")
        NULL
      }
    }, error = function(e) {
      logger::log_warn("ADVI {strat$label} failed: {e$message}")
      NULL
    })

    if (!is.null(fit_vb)) {
      fit <- fit_vb
      sampling_method <- paste0("variational_", strat$label)
      break
    }
  }

  if (is.null(fit)) {
    sampling_method <- "map_only"
    logger::log_warn("All methods failed; using MAP estimates")
  }
}

logger::log_info("Fitting complete (method: {sampling_method})")

# ============================================================================
# Save Results (same format as original)
# ============================================================================

summ <- tryCatch({
  if (sampling_method == "map_only") {
    map_estimates %>%
      rename(median = estimate) %>%
      mutate(q5 = NA_real_, q95 = NA_real_, rhat = NA_real_, ess_bulk = NA_real_)
  } else {
    s <- fit$summary()
    if (!"rhat" %in% names(s)) s$rhat <- NA_real_
    if (!"ess_bulk" %in% names(s)) s$ess_bulk <- NA_real_
    if (!"q5" %in% names(s)) s$q5 <- NA_real_
    if (!"q95" %in% names(s)) s$q95 <- NA_real_
    s
  }
}, error = function(e) {
  logger::log_warn("Summary failed: {e$message}")
  map_estimates %>% rename(median = estimate) %>%
    mutate(q5 = NA_real_, q95 = NA_real_, rhat = NA_real_, ess_bulk = NA_real_)
})

if (sampling_method != "map_only") {
  draws_df <- tryCatch(fit$draws(format = "df"),
    error = function(e) NULL)
  if (!is.null(draws_df)) {
    saveRDS(draws_df, file.path(output_dir, "diameter_growth_samples.rds"))
  }
}

writeLines(sampling_method, file.path(output_dir, "diameter_growth_method.txt"))
write_csv(summ, file.path(output_dir, "diameter_growth_summary.csv"))

posterior_summary <- summ %>%
  select(variable, median, q5, q95, any_of(c("rhat", "ess_bulk"))) %>%
  rename(p50 = median, p05 = q5, p95 = q95) %>%
  mutate(
    rhat = if ("rhat" %in% names(.)) rhat else NA_real_,
    ess_bulk = if ("ess_bulk" %in% names(.)) ess_bulk else NA_real_,
    ci_width = p95 - p05,
    converged = ifelse(is.na(rhat), NA, rhat <= 1.01)
  )
write_csv(posterior_summary, file.path(output_dir, "diameter_growth_posterior.csv"))

cat("\n=========================================\n")
cat("DG Fitting Complete (HMC-small v3)\n")
cat("=========================================\n")
cat("Variant:", variant, "\n")
cat("Method:", sampling_method, "\n")
cat("N:", N, ", Species:", N_species, "\n")
cat("Output:", output_dir, "\n\n")

logger::log_info("DG HMC-small complete for {variant} (method: {sampling_method})")