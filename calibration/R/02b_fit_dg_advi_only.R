#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Wykoff DG Model (ADVI-only, v2)
# For variants where HMC is computationally intractable.
# Skips HMC entirely: MAP -> fullrank ADVI with relaxed settings.
#
# Changes from 02_fit_diameter_growth.R:
#   1. Uses wykoff_dg_v2.stan (sigma lower bounds to prevent zero-scale)
#   2. Skips HMC entirely (MAP -> ADVI)
#   3. Uses fullrank ADVI (captures posterior correlations, better than meanfield)
#   4. Reduces max observations to 15000 for faster convergence
#   5. Increases ADVI iterations to 100000 with relaxed tolerance
#   6. Adds eta tuning with multiple restarts
#
# Usage: Rscript 02b_fit_dg_advi_only.R --variant bc
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

variant <- "ca"  # Default
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

# Resolve project root from the environment first. The previous form passed the
# sys.frame() fallback as the Sys.getenv() default, which R evaluates eagerly and
# which crashes under `Rscript` (no `ofile` frame) even when FVS_PROJECT_ROOT is
# set. Evaluate the fallback only when the env var is absent.
project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (!nzchar(project_root)) {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
calibration_dir <- file.path(project_root, "calibration")
processed_data_dir <- file.path(calibration_dir, "data", "processed")
output_dir <- file.path(calibration_dir, "output", "variants", variant)
stan_dir <- file.path(calibration_dir, "stan")

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Reduced max obs for ADVI tractability
MAX_OBS <- as.integer(Sys.getenv("FVS_MAX_OBS", "15000"))

# Set up logging
log_file <- file.path(calibration_dir, "logs", paste0("02b_fit_dg_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting DG ADVI-only fitting for variant {variant}")
logger::log_info("Max observations: {MAX_OBS}")

# ============================================================================
# Load Data
# ============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

data_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (!file.exists(data_file)) {
  logger::log_error("Data file not found: {data_file}")
  stop("Data file not found: ", data_file)
}

dg_data <- read_csv(data_file, show_col_types = FALSE) %>%
  as_tibble()

logger::log_info("Loaded {nrow(dg_data)} observations")

# If species codes not provided, fit all species in the variant
if (is.null(species_codes)) {
  species_codes <- unique(dg_data$SPCD)
  logger::log_info("No species specified; will fit {length(species_codes)} species")
}

# Filter to requested species
dg_data <- dg_data %>%
  filter(SPCD %in% species_codes)

logger::log_info("After species filter: {nrow(dg_data)} observations")

# ============================================================================
# Load FVS Config for Priors
# ============================================================================

logger::log_info("Loading FVS parameters for priors...")

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
# Prepare Data for Stan Model
# ============================================================================

logger::log_info("Preparing data for Stan model...")

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

# Subsample: reduced from 30K to 15K for ADVI tractability
if (nrow(dg_model) > MAX_OBS) {
  logger::log_info("Subsampling from {nrow(dg_model)} to {MAX_OBS} (stratified by species)")
  set.seed(42)
  subsample_prop <- min(1, MAX_OBS / nrow(dg_model) * 1.1)
  dg_model <- dg_model %>%
    group_by(species_idx) %>%
    slice_sample(prop = subsample_prop) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), MAX_OBS))
  logger::log_info("After subsampling: {nrow(dg_model)} observations")
}

# Recompute consecutive indices
dg_model <- dg_model %>%
  mutate(species_idx = as.integer(factor(species_idx)))

N <- nrow(dg_model)
N_species <- max(dg_model$species_idx)

logger::log_info("N = {N}, N_species = {N_species}")

# Persist the dense species_idx -> SPCD crosswalk so downstream tooling (e.g.
# calibration/R/multipliers.R) can map b0[i] to the correct FVS species without
# replaying the stochastic, subsample-dependent data prep. Without this file the
# dense index is not recoverable after the processed data is regenerated.
dg_species_index <- dg_model %>%
  dplyr::distinct(species_idx, SPCD) %>%
  dplyr::arrange(species_idx)
utils::write.csv(dg_species_index,
                 file.path(output_dir, "diameter_growth_species_index.csv"),
                 row.names = FALSE)

prior_b0_species <- rep(prior_means["b1_ref"] * 0.5, N_species)

stan_data <- list(
  N = N,
  N_species = N_species,
  ln_DDS = dg_model$ln_DDS,
  ln_DBH = dg_model$ln_DBH_std,
  DBH_sq = dg_model$DBH_sq_std,
  ln_SI = dg_model$ln_SI_std,
  SLOPE = dg_model$SLOPE_std,
  SLOPE_sq = dg_model$SLOPE_sq,
  SLOPE_SASP = dg_model$SLOPE_SASP_std,
  SLOPE_CASP = dg_model$SLOPE_CASP_std,
  ELEV = dg_model$ELEV_std,
  ELEV_sq = dg_model$ELEV_sq,
  CR = dg_model$CR,
  CR_sq = dg_model$CR_sq,
  BAL = dg_model$BAL_std,
  BA = dg_model$BA_std,
  species_id = dg_model$species_idx,
  prior_b0 = prior_b0_species
)

# ============================================================================
# Compile Stan Model (v2 with sigma lower bounds)
# ============================================================================

logger::log_info("Compiling Stan model v2...")

stan_file <- file.path(stan_dir, "wykoff_dg_v2.stan")

# Fall back to original if v2 not available
if (!file.exists(stan_file)) {
  stan_file <- file.path(stan_dir, "wykoff_dg.stan")
  logger::log_warn("v2 Stan model not found, using original")
}

mod <- cmdstanr::cmdstan_model(stan_file, dir = output_dir)
logger::log_info("Stan model compiled successfully")

# ============================================================================
# Stage 0: OLS for initial values
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

logger::log_info("OLS intercept: {round(ols_coef[1], 3)}, sigma: {round(ols_sigma, 3)}")

# ============================================================================
# Stage 1: MAP optimization
# ============================================================================

logger::log_info("Stage 1: MAP optimization...")

fit_map <- tryCatch({
  mod$optimize(
    data = stan_data,
    init = init_fn,
    algorithm = "lbfgs",
    iter = 10000,
    tol_rel_grad = 1e-8,
    refresh = 100
  )
}, error = function(e) {
  logger::log_warn("MAP optimization failed: {e$message}")
  NULL
})

if (!is.null(fit_map)) {
  map_estimates <- fit_map$summary()
  logger::log_info("MAP complete. Log posterior: {round(fit_map$lp(), 2)}")

  map_file <- file.path(output_dir, "diameter_growth_map.csv")
  write_csv(map_estimates, map_file)

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
      sigma = max(unname(map_vals["sigma"]), 0.01)
    )
  }
} else {
  init_from_map <- init_fn
}

# ============================================================================
# Stage 2: ADVI (skip HMC entirely for these problematic variants)
# ============================================================================

logger::log_info("Stage 2: Fullrank ADVI (skipping HMC)...")

sampling_method <- "none"

# Strategy: try fullrank first, fall back to meanfield, try multiple eta values
advi_strategies <- list(
  list(algorithm = "fullrank", iter = 100000, tol = 0.005, eta = NULL, label = "fullrank_auto"),
  list(algorithm = "fullrank", iter = 100000, tol = 0.01,  eta = 0.1,  label = "fullrank_eta01"),
  list(algorithm = "meanfield", iter = 100000, tol = 0.005, eta = NULL, label = "meanfield_auto"),
  list(algorithm = "meanfield", iter = 100000, tol = 0.01,  eta = 0.01, label = "meanfield_eta001")
)

fit <- NULL
for (strat in advi_strategies) {
  logger::log_info("Trying ADVI strategy: {strat$label}...")

  vb_args <- list(
    data = stan_data,
    init = init_from_map,
    algorithm = strat$algorithm,
    iter = strat$iter,
    tol_rel_obj = strat$tol,
    output_samples = 4000,
    refresh = 1000
  )

  # Add eta if specified
  if (!is.null(strat$eta)) {
    vb_args$eta <- strat$eta
    vb_args$adapt_engaged <- FALSE
  }

  fit_vb <- tryCatch({
    do.call(mod$variational, vb_args)
  }, error = function(e) {
    logger::log_warn("ADVI strategy {strat$label} failed: {e$message}")
    NULL
  })

  if (!is.null(fit_vb)) {
    fit <- fit_vb
    sampling_method <- paste0("variational_", strat$label)
    logger::log_info("ADVI strategy {strat$label} succeeded")
    break
  }
}

# Last resort: use MAP estimates directly
if (is.null(fit)) {
  sampling_method <- "map_only"
  logger::log_warn("All ADVI strategies failed; using MAP estimates only")
}

logger::log_info("Model fitting complete (method: {sampling_method})")

# ============================================================================
# Extract and Save Results (same format as original script)
# ============================================================================

logger::log_info("Computing summaries (method: {sampling_method})...")

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
  logger::log_warn("Summary computation failed: {e$message}")
  if (!is.null(fit_map)) {
    map_estimates %>%
      rename(median = estimate) %>%
      mutate(q5 = NA_real_, q95 = NA_real_, rhat = NA_real_, ess_bulk = NA_real_)
  } else {
    tibble(variable = character(), median = numeric())
  }
})

# Save posterior samples
if (sampling_method != "map_only") {
  draws_df <- tryCatch(
    fit$draws(format = "df"),
    error = function(e) { logger::log_warn("Draw extraction failed: {e$message}"); NULL }
  )
  if (!is.null(draws_df)) {
    saveRDS(draws_df, file.path(output_dir, "diameter_growth_samples.rds"))
    logger::log_info("Saved posterior samples")
  }
}

# Save method used
writeLines(sampling_method, file.path(output_dir, "diameter_growth_method.txt"))

# Save summary (this is what triggers "completion" in the monitoring pipeline)
summary_file <- file.path(output_dir, "diameter_growth_summary.csv")
write_csv(summ, summary_file)
logger::log_info("Saved summary statistics to {summary_file}")

# Posterior medians and CIs
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
logger::log_info("Saved posterior summaries")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("=========================================\n")
cat("DG Model Fitting Complete (ADVI-only v2)\n")
cat("=========================================\n")
cat("Variant:", variant, "\n")
cat("Inference method:", sampling_method, "\n")
cat("Species fitted:", N_species, "\n")
cat("Total observations:", N, "\n")
cat("\nModel summary (first 20 parameters):\n")
tryCatch(print(head(summ, 20)), error = function(e) cat("(summary unavailable)\n"))
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("DG ADVI-only fitting complete for variant {variant} (method: {sampling_method})")