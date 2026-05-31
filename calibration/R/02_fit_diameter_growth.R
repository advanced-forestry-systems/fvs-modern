#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Wykoff Diameter Growth Model
# Bayesian hierarchical regression using CmdStanR and brms
#
# Usage: Rscript calibration/R/02_fit_diameter_growth.R --variant ca --species 122
# If no args provided, runs a quick test fit for variant ca
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

# Set up logging
log_file <- file.path(calibration_dir, "logs", paste0("02_fit_dg_", variant, ".log"))
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting diameter growth model fitting for variant {variant}")

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

# Extract growth parameters from config
# The config has WEIBB1 (intercept adjustment) and DG* parameters
fvs_params <- config$categories$other

# Map parameter names to our model
param_map <- list(
  b1_ref = "DGLD",           # ln(DBH) coefficient
  b2_ref = "DGDS",           # DBH^2 coefficient
  b3_ref = "DGSITE",         # ln(SI) coefficient
  b4_ref = "DGSLOP",         # SLOPE coefficient
  b5_ref = "DGSLSQ",         # SLOPE^2 coefficient
  b6_ref = "DGSASP",         # SLOPE*sin(ASP) coefficient
  b7_ref = "DGCASP",         # SLOPE*cos(ASP) coefficient
  b8_ref = "DGEL",           # ELEV coefficient
  b9_ref = "DGELSQ",         # ELEV^2 coefficient
  b10_ref = "DGCR",          # CR coefficient
  b11_ref = "DGCRSQ",        # CR^2 coefficient
  b12_ref = "DGBAL",         # BAL coefficient
  b13_ref = "DGPCCF"         # BA/CCF coefficient
)

# Extract reference values (use first species as reference)
prior_means <- sapply(param_map, function(param_name) {
  vals <- fvs_params[[param_name]]
  if (is.null(vals) || length(vals) == 0) {
    return(0)
  }
  # Use mean of all species' values
  mean(unlist(vals), na.rm = TRUE)
})

logger::log_info("Extracted prior means from FVS config")

# ============================================================================
# Prepare Data for Stan Model
# ============================================================================

logger::log_info("Preparing data for Stan model...")

# Create numeric indices for grouping variables
dg_model <- dg_data %>%
  # Create factor levels for species
  mutate(species_idx = as.numeric(factor(SPCD)))

# Standardize ALL continuous predictors for stable HMC sampling
dg_model <- dg_model %>%
  mutate(
    # Standardize core predictors
    ln_DBH_std = as.numeric(scale(ln_DBH)),
    SLOPE_std = as.numeric(scale(SLOPE)),
    ELEV_std = as.numeric(scale(ELEV)),
    DBH_sq_std = as.numeric(scale(DBH_sq)),
    ln_SI_std = as.numeric(scale(ln_SI)),

    # Quadratic terms (on standardized scale)
    SLOPE_sq = SLOPE_std^2,
    ELEV_sq = ELEV_std^2,

    # Standardize interaction terms
    SLOPE_SASP_std = as.numeric(scale(SLOPE_SASP)),
    SLOPE_CASP_std = as.numeric(scale(SLOPE_CASP)),

    # Crown ratio: scale to 0 to 1
    CR = CR_pct / 100,
    CR_sq = CR^2,

    # Standardize competition variables
    BAL_std = as.numeric(scale(BAL)),
    BA_std = as.numeric(scale(BA))
  )

# Remove any rows with NA in key predictors after transformation
dg_model <- dg_model %>%
  filter(
    is.finite(ln_DBH_std), is.finite(SLOPE_std), is.finite(ELEV_std),
    is.finite(ln_DDS), is.finite(ln_SI_std), is.finite(DBH_sq_std),
    is.finite(SLOPE_SASP_std), is.finite(SLOPE_CASP_std),
    !is.na(CR), is.finite(BAL_std), is.finite(BA_std)
  )

# Save standardization parameters for back transformation
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

# Subsample for computational tractability (species stratified to preserve rare species)
max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(dg_model) > max_n) {
  logger::log_info("Subsampling from {nrow(dg_model)} to {max_n} observations (stratified by species)")
  set.seed(42)
  subsample_prop <- min(1, max_n / nrow(dg_model) * 1.1)
  dg_model <- dg_model %>%
    group_by(species_idx) %>%
    slice_sample(prop = subsample_prop) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), max_n))
  logger::log_info("After subsampling: {nrow(dg_model)} observations")
}

# Recompute consecutive indices after filtering (gaps cause Stan indexing errors)
dg_model <- dg_model %>%
  mutate(species_idx = as.integer(factor(species_idx)))

# Extract Stan data list
N <- nrow(dg_model)
N_species <- max(dg_model$species_idx)

logger::log_info("N = {N}, N_species = {N_species}")

# Prior: use mean FVS intercept value for each species
prior_b0_species <- rep(prior_means["b1_ref"] * 0.5, N_species)

stan_data <- list(
  N = N,
  N_species = N_species,

  # Response
  ln_DDS = dg_model$ln_DDS,

  # Predictors (all standardized for stable HMC)
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

  # Group indices
  species_id = dg_model$species_idx,

  # Priors
  prior_b0 = prior_b0_species
)

# ============================================================================
# Compile Stan Model
# ============================================================================

logger::log_info("Compiling Stan model...")

stan_file <- file.path(stan_dir, "wykoff_dg.stan")

tryCatch({
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    dir = output_dir  # Use variant dir for compiled model
  )
  logger::log_info("Stan model compiled successfully")
}, error = function(e) {
  logger::log_error("Failed to compile Stan model: {e$message}")
  stop(e)
})

# ============================================================================
# Fit Model
# ============================================================================

logger::log_info("Fitting diameter growth model (multi-strategy)...")

# ---- Stage 0: OLS for initial values ----
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

# ---- Stage 1: MAP optimization (fast, always works) ----
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
  logger::log_info("MAP optimization converged. Log posterior: {round(fit_map$lp(), 2)}")

  # Save MAP estimates
  map_file <- file.path(output_dir, "diameter_growth_map.csv")
  write_csv(map_estimates, map_file)
  logger::log_info("Saved MAP estimates to {map_file}")

  # Build init from MAP for subsequent sampling
  map_vals <- setNames(map_estimates$estimate, map_estimates$variable)
  init_from_map <- function() {
    z_b0_names <- paste0("z_b0[", 1:N_species, "]")
    list(
      mu_b0 = unname(map_vals["mu_b0"]),
      sigma_b0 = max(unname(map_vals["sigma_b0"]), 0.01),
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

# ---- Stage 2: Try HMC sampling with MAP initialization ----
logger::log_info("Stage 2: HMC sampling with MAP initialization...")

sampling_method <- "none"  # track which method produced the final fit

fit <- tryCatch({
  result <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = min(4, parallel::detectCores()),
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 100,
    show_messages = FALSE,
    init = init_from_map
  )

  # Check if sampling actually converged
  diag <- result$diagnostic_summary()
  n_divergent <- sum(diag$num_divergent)
  n_treedepth <- sum(diag$num_max_treedepth)
  total_transitions <- 4 * 1000  # chains * iter_sampling

  # Accept if < 50% treedepth violations and < 5% divergences
  if (n_treedepth < total_transitions * 0.5 &&
      n_divergent < total_transitions * 0.05) {
    sampling_method <<- "hmc"
    logger::log_info("HMC sampling converged: {n_divergent} divergent, {n_treedepth} treedepth")
    result
  } else {
    logger::log_warn("HMC sampling poor: {n_divergent} divergent, {n_treedepth}/{total_transitions} treedepth")
    logger::log_info("Falling back to variational inference...")
    NULL
  }
}, error = function(e) {
  logger::log_warn("HMC sampling failed: {e$message}")
  logger::log_info("Falling back to variational inference...")
  NULL
})

# ---- Stage 3: Variational inference fallback ----
if (is.null(fit) || sampling_method == "none") {
  logger::log_info("Stage 3: Variational inference (ADVI)...")

  fit_vb <- tryCatch({
    mod$variational(
      data = stan_data,
      init = init_from_map,
      algorithm = "meanfield",
      iter = 50000,
      tol_rel_obj = 0.001,
      output_samples = 4000,
      refresh = 500
    )
  }, error = function(e) {
    logger::log_warn("Variational inference failed: {e$message}")
    NULL
  })

  if (!is.null(fit_vb)) {
    fit <- fit_vb
    sampling_method <- "variational"
    logger::log_info("Variational inference complete")
  } else {
    # Last resort: use MAP estimates directly (no uncertainty)
    sampling_method <- "map_only"
    logger::log_warn("All sampling methods failed; using MAP estimates only")
  }
}

logger::log_info("Model fitting complete (method: {sampling_method})")

# ============================================================================
# Check Convergence Diagnostics
# ============================================================================

logger::log_info("Computing convergence diagnostics (method: {sampling_method})...")

if (sampling_method == "hmc") {
  diagnostics <- fit$diagnostic_summary()
  logger::log_info("HMC diagnostics: divergent={sum(diagnostics$num_divergent)}, treedepth={sum(diagnostics$num_max_treedepth)}")
}

# Get summary statistics
summ <- tryCatch({
  if (sampling_method == "map_only") {
    # Use MAP estimates as summary
    map_estimates %>%
      rename(median = estimate) %>%
      mutate(q5 = NA_real_, q95 = NA_real_, rhat = NA_real_, ess_bulk = NA_real_)
  } else {
    s <- fit$summary()
    # Ensure required columns exist (variational may lack rhat/ess_bulk)
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

# Check Rhat and ESS (only meaningful for HMC)
if (sampling_method == "hmc" && "rhat" %in% names(summ)) {
  rhat_issues <- summ %>%
    filter(rhat > 1.01 & !is.na(rhat))

  if (nrow(rhat_issues) > 0) {
    logger::log_warn("Found {nrow(rhat_issues)} parameters with Rhat > 1.01")
    logger::log_info("Parameters with poor convergence:")
    logger::log_info("{paste(head(rhat_issues$variable, 10), collapse = ', ')}")
  } else {
    logger::log_info("All parameters converged (Rhat <= 1.01)")
  }
}

# ============================================================================
# Extract and Save Posterior Samples
# ============================================================================

logger::log_info("Extracting posterior samples (method: {sampling_method})...")

if (sampling_method != "map_only") {
  # Get draws in tibble format
  draws_df <- tryCatch(
    fit$draws(format = "df"),
    error = function(e) { logger::log_warn("Draw extraction failed: {e$message}"); NULL }
  )

  if (!is.null(draws_df)) {
    draws_file <- file.path(output_dir, "diameter_growth_samples.rds")
    saveRDS(draws_df, draws_file)
    logger::log_info("Saved posterior samples to {draws_file}")
  }
} else {
  draws_df <- NULL
  logger::log_info("MAP only: no posterior samples to save")
}

# Save method used
writeLines(sampling_method, file.path(output_dir, "diameter_growth_method.txt"))

# Save summary
summary_file <- file.path(output_dir, "diameter_growth_summary.csv")
write_csv(summ, summary_file)
logger::log_info("Saved summary statistics to {summary_file}")

# ============================================================================
# Extract Posterior Medians and Credible Intervals
# ============================================================================

logger::log_info("Computing posterior summaries...")

# Posterior medians and 95% CIs
# Use any_of() for rhat/ess_bulk since variational inference does not produce these
posterior_summary <- summ %>%
  select(variable, median, q5, q95, any_of(c("rhat", "ess_bulk"))) %>%
  rename(
    p50 = median,
    p05 = q5,
    p95 = q95
  ) %>%
  mutate(
    rhat = if ("rhat" %in% names(.)) rhat else NA_real_,
    ess_bulk = if ("ess_bulk" %in% names(.)) ess_bulk else NA_real_,
    ci_width = p95 - p05,
    converged = ifelse(is.na(rhat), NA, rhat <= 1.01)
  )

# Save posterior summary
posterior_file <- file.path(output_dir, "diameter_growth_posterior.csv")
write_csv(posterior_summary, posterior_file)
logger::log_info("Saved posterior summaries to {posterior_file}")

# ============================================================================
# Generate Diagnostic Plots
# ============================================================================

logger::log_info("Generating diagnostic plots...")

# Traceplot for key parameters (HMC only)
if (sampling_method == "hmc") {
  tryCatch({
    key_params <- c("mu_b0", "sigma_b0", "b1", "b3", "b10", "b12", "b13", "sigma")
    pdf(file.path(output_dir, "diameter_growth_traceplots.pdf"), width = 12, height = 8)
    print(bayesplot::mcmc_trace(
      fit$draws(),
      pars = intersect(key_params, summ$variable),
      facet_args = list(nrow = 3)
    ))
    dev.off()
    logger::log_info("Saved traceplots to PDF")
  }, error = function(e) {
    logger::log_warn("Traceplot generation failed: {e$message}")
  })
}

# Posterior predictive check (if draws available)
if (sampling_method != "map_only") {
  tryCatch({
    y_rep <- fit$draws("ln_DDS_rep", format = "matrix")
    M <- min(N, 1000)
    ppc_file <- file.path(output_dir, "diameter_growth_ppc.pdf")
    pdf(ppc_file, width = 10, height = 6)
    print(bayesplot::ppc_dens_overlay(
      y = stan_data$ln_DDS[1:M],
      yrep = y_rep[1:100, ]
    ))
    dev.off()
    logger::log_info("Saved posterior predictive check to {ppc_file}")
  }, error = function(e) {
    logger::log_warn("PPC plot failed: {e$message}")
  })
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("=========================================\n")
cat("Diameter Growth Model Fitting Complete\n")
cat("=========================================\n")
cat("Variant:", variant, "\n")
cat("Inference method:", sampling_method, "\n")
cat("Species fitted:", N_species, "\n")
cat("Total observations:", N, "\n")
cat("\nModel summary (first 20 parameters):\n")
tryCatch(print(head(summ, 20)), error = function(e) cat("(summary unavailable)\n"))
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Diameter growth model fitting complete for variant {variant} (method: {sampling_method})")
