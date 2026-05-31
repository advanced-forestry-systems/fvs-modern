#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Mortality Model
# Bayesian logistic regression for tree survival probability with multi-strategy inference
#
# Model: P(dead) = logit^-1(b0 + b1*DBH + b2*DBH^2 + b3*BAL + b4*CR + b5*SI + b6*BA)
# With random effects for species
#
# Inference strategies:
#   1. HMC (CmdStanR backend) with reduced parameters
#   2. Variational inference (meanfield) as fallback
#
# Usage: Rscript calibration/R/04_fit_mortality.R --variant ca
#

library(tidyverse)
library(data.table)
library(brms)
library(cmdstanr)
library(jsonlite)
library(posterior)
library(bayesplot)
library(tidybayes)
library(logger)

# ============================================================================
# Parse Command Line Arguments
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
variant <- "ca"

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variant <- args[i + 1]
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

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- file.path(calibration_dir, "logs", paste0("04_fit_mortality_", variant, ".log"))
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting mortality model fitting for variant {variant}")

# ============================================================================
# Load Data
# ============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

# Try dedicated mortality.csv first, fall back to diameter_growth.csv
mortality_file <- file.path(processed_data_dir, variant, "mortality.csv")
diameter_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (file.exists(mortality_file)) {
  logger::log_info("Using dedicated mortality.csv file")
  data_file <- mortality_file
} else if (file.exists(diameter_file)) {
  logger::log_info("Mortality.csv not found, falling back to diameter_growth.csv")
  data_file <- diameter_file
} else {
  logger::log_error("Neither mortality.csv nor diameter_growth.csv found in {file.path(processed_data_dir, variant)}")
  stop("Data files not found")
}

# Load and prepare data
mortality_data <- read_csv(data_file, show_col_types = FALSE) %>%
  as_tibble()

# Harmonize column names across file formats
# mortality.csv: DIA, died
# diameter_growth.csv: DIA_t1, DIA_t2, STATUSCD_t2
if ("DIA" %in% names(mortality_data) & !("DIA_t1" %in% names(mortality_data))) {
  mortality_data <- mortality_data %>% rename(DIA_t1 = DIA)
}

# Create survival indicator from whatever columns are available
if ("survived" %in% names(mortality_data)) {
  # Already has survived column
} else if ("died" %in% names(mortality_data)) {
  mortality_data <- mortality_data %>% mutate(survived = !as.logical(died))
} else if ("DIA_t2" %in% names(mortality_data)) {
  mortality_data <- mortality_data %>% mutate(survived = !is.na(DIA_t2) & DIA_t2 > 0)
} else {
  stop("Cannot determine survival status: no 'survived', 'died', or 'DIA_t2' column")
}

# Handle CR_pct: in mortality.csv it is stored as CR * 100 (e.g., 35% = 3500)
# If CR_pct is missing but CR exists, create it
if (!("CR_pct" %in% names(mortality_data)) & "CR" %in% names(mortality_data)) {
  mortality_data <- mortality_data %>% mutate(CR_pct = CR * 100)
}

# Standardize predictors
mortality_data <- mortality_data %>%
  mutate(
    DBH_std = scale(DIA_t1)[, 1],
    CR_std = scale(CR_pct / 100)[, 1],
    BAL_std = scale(BAL)[, 1],
    SI_std = scale(log(pmax(SI, 1)))[, 1],
    BA_std = scale(BA)[, 1],
    SPCD = as.factor(SPCD)
  ) %>%
  filter(!is.na(survived), !is.na(DBH_std), !is.na(BAL_std),
         !is.na(CR_std), !is.na(SI_std), !is.na(BA_std))

logger::log_info("Loaded {nrow(mortality_data)} tree observations")
logger::log_info("Mortality rate: {round(mean(!mortality_data$survived) * 100, 2)}%")

# Subsample for computational tractability (stratified by species)
max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(mortality_data) > max_n) {
  logger::log_info("Subsampling from {nrow(mortality_data)} to {max_n} observations (stratified by species)")
  set.seed(42)
  mortality_data <- mortality_data %>%
    group_by(SPCD) %>%
    slice_sample(prop = min(1, max_n / nrow(mortality_data) * 1.1)) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), max_n))
  logger::log_info("After subsampling: {nrow(mortality_data)} observations")
}

# ============================================================================
# Multi-Strategy Mortality Model Fitting
# ============================================================================

# Define model formula
model_formula <- survived ~ DBH_std + I(DBH_std^2) + BAL_std + CR_std + SI_std + BA_std +
  (1 | SPCD)

fit_mort <- NULL
inference_method <- "none"

# ============================================================================
# Strategy 1: HMC with CmdStanR Backend (Primary)
# ============================================================================

logger::log_info("Attempting Strategy 1: HMC (CmdStanR backend) with reduced parameters...")

fit_result <- tryCatch(
  {
    fit <- brm(
      model_formula,
      data = mortality_data,
      family = bernoulli(link = "logit"),
      chains = 4,
      iter = 1500,
      warmup = 500,
      cores = parallel::detectCores(),
      backend = "cmdstanr",
      refresh = 0,
      control = list(adapt_delta = 0.90, max_treedepth = 12)
    )
    list(fit = fit, success = TRUE, error = NULL)
  },
  error = function(e) {
    logger::log_warn("Strategy 1 (HMC) failed: {as.character(e)}")
    list(fit = NULL, success = FALSE, error = as.character(e))
  }
)

if (fit_result$success) {
  fit_mort <- fit_result$fit
  inference_method <- "HMC_CmdStanR"
  logger::log_info("Strategy 1 succeeded: using HMC with CmdStanR")
} else {
  logger::log_warn("Strategy 1 failed, attempting Strategy 2...")

  # ============================================================================
  # Strategy 2: Variational Inference (Fallback)
  # ============================================================================

  logger::log_info("Attempting Strategy 2: Variational Inference (meanfield)...")

  vi_result <- tryCatch(
    {
      fit <- brm(
        model_formula,
        data = mortality_data,
        family = bernoulli(link = "logit"),
        algorithm = "meanfield",
        cores = parallel::detectCores(),
        backend = "cmdstanr",
        refresh = 0,
      )
      list(fit = fit, success = TRUE, error = NULL)
    },
    error = function(e) {
      logger::log_error("Strategy 2 (Variational) also failed: {as.character(e)}")
      list(fit = NULL, success = FALSE, error = as.character(e))
    }
  )

  if (vi_result$success) {
    fit_mort <- vi_result$fit
    inference_method <- "VariationalInference_meanfield"
    logger::log_info("Strategy 2 succeeded: using Variational Inference")
  } else {
    logger::log_error("Both inference strategies failed")
    stop("Mortality model fitting failed with both HMC and variational inference")
  }
}

logger::log_info("Mortality model fitting complete using method: {inference_method}")

# ============================================================================
# Check Convergence
# ============================================================================

logger::log_info("Checking convergence diagnostics...")

summ_mort <- tryCatch(
  posterior_summary(fit_mort),
  error = function(e) {
    logger::log_warn("posterior_summary failed: {e$message}. Using summary() instead.")
    as.data.frame(fixef(fit_mort))
  }
)

# Only check Rhat for HMC (variational inference returns NA for Rhat)
if (inference_method == "HMC_CmdStanR" & "Rhat" %in% colnames(summ_mort)) {
  rhat_vals <- summ_mort[, "Rhat"]
  n_bad_rhat <- sum(rhat_vals > 1.01, na.rm = TRUE)
  if (n_bad_rhat > 0) {
    logger::log_warn("{n_bad_rhat} parameters have Rhat > 1.01")
  }
}

# ============================================================================
# Extract and Save Results
# ============================================================================

logger::log_info("Extracting posterior samples and preparing outputs...")

# Get posterior draws
draws_mort <- as_draws_df(fit_mort)

# Save raw posterior samples
draws_file <- file.path(output_dir, "mortality_samples.rds")
saveRDS(draws_mort, draws_file)
logger::log_info("Saved posterior samples to {draws_file}")

# Use brms summary which is more reliable than posterior_summary
brms_summ <- summary(fit_mort)

# Extract fixed effects and random effects into a combined data frame
fixed_df <- as_tibble(brms_summ$fixed, rownames = "variable")
# Standardize column names across brms versions
names(fixed_df) <- gsub("^Estimate$", "p50", names(fixed_df))
names(fixed_df) <- gsub("^l-95% CI$", "p025", names(fixed_df))
names(fixed_df) <- gsub("^u-95% CI$", "p975", names(fixed_df))

summary_file <- file.path(output_dir, "mortality_summary.csv")
write_csv(fixed_df, summary_file)
logger::log_info("Saved summary to {summary_file}")

# Build posterior summary from draws for more control
# Use explicit posterior:: namespace to avoid rstan masking issues
# posterior::summarise_draws with quantile() inherits names like "5%" instead of "p05"
# so we use unname() to force our chosen column names
draws_summ <- tryCatch({
  posterior::summarise_draws(draws_mort,
    p50 = function(x) unname(median(x)),
    p05 = function(x) unname(quantile(x, 0.05)),
    p95 = function(x) unname(quantile(x, 0.95)),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
}, error = function(e) {
  logger::log_warn("summarise_draws failed: {e$message}. Using manual summary.")
  draws_mat <- as.matrix(draws_mort)
  tibble(
    variable = colnames(draws_mat),
    p50 = apply(draws_mat, 2, function(x) unname(median(x))),
    p05 = apply(draws_mat, 2, function(x) unname(quantile(x, 0.05))),
    p95 = apply(draws_mat, 2, function(x) unname(quantile(x, 0.95))),
    rhat = NA_real_,
    ess_bulk = NA_real_
  )
})

# Normalize column names: posterior 1.6+ may still rename to "5%" / "95%"
names(draws_summ) <- gsub("^5%$", "p05", names(draws_summ))
names(draws_summ) <- gsub("^95%$", "p95", names(draws_summ))
if (!"variable" %in% names(draws_summ)) {
  names(draws_summ)[1] <- "variable"
}

posterior_summary_df <- draws_summ %>%
  mutate(
    ci_width = p95 - p05,
    converged = if_else(!is.na(rhat), rhat < 1.05, NA)
  )

posterior_file <- file.path(output_dir, "mortality_posterior.csv")
write_csv(posterior_summary_df, posterior_file)
logger::log_info("Saved posterior summary to {posterior_file}")

# Prepare point estimates (posterior median)
map_df <- posterior_summary_df %>%
  select(variable, estimate = p50)

map_file <- file.path(output_dir, "mortality_map.csv")
write_csv(map_df, map_file)
logger::log_info("Saved point estimates to {map_file}")

# ============================================================================
# Posterior Predictive Checks
# ============================================================================

logger::log_info("Generating diagnostic plots...")

# PPC for binary outcome (only for HMC, skip for variational)
if (inference_method == "HMC_CmdStanR") {
  tryCatch(
    {
      pdf(file.path(output_dir, "mortality_ppc.pdf"), width = 10, height = 6)
      print(pp_check(fit_mort, ndraws = 100, type = "error_binned"))
      dev.off()
      logger::log_info("Generated PPC plot")
    },
    error = function(e) {
      logger::log_warn("Could not generate PPC plot: {as.character(e)}")
    }
  )
}

# Conditional effects (works for both inference methods)
tryCatch(
  {
    pdf(file.path(output_dir, "mortality_effects.pdf"), width = 12, height = 8)
    print(plot(conditional_effects(fit_mort)))
    dev.off()
    logger::log_info("Generated conditional effects plot")
  },
  error = function(e) {
    logger::log_warn("Could not generate conditional effects plot: {as.character(e)}")
  }
)

logger::log_info("Saved diagnostic plots")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("======================================\n")
cat("Mortality Model Fitting Complete\n")
cat("======================================\n")
cat("Variant:", variant, "\n")
cat("Total observations:", nrow(mortality_data), "\n")
cat("Mortality rate:", round(mean(!mortality_data$survived) * 100, 2), "%\n")
cat("Inference method:", inference_method, "\n")
cat("\nModel:\n")
print(fit_mort)
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Mortality model fitting complete with {inference_method} inference")
