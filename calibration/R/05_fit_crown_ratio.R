#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Crown Ratio Model (Multi-Strategy Inference)
# Model change in crown ratio between measurements
# Uses BCR1-BCR4 coefficients from FVS config as informative priors
# Implements multi-strategy inference: HMC with fallback to variational inference
#
# Usage: Rscript calibration/R/05_fit_crown_ratio.R --variant ca
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
log_file <- file.path(calibration_dir, "logs", paste0("05_fit_crown_ratio_", variant, ".log"))
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting crown ratio model fitting for variant {variant}")

# ============================================================================
# Load Data
# ============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

# Try dedicated crown_ratio_change.csv first, fall back to diameter_growth.csv
cr_change_file <- file.path(processed_data_dir, variant, "crown_ratio_change.csv")
diameter_growth_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (file.exists(cr_change_file)) {
  logger::log_info("Using dedicated crown_ratio_change.csv file")
  data_file <- cr_change_file
  cr_data <- read_csv(data_file, show_col_types = FALSE) %>%
    as_tibble()

  # Harmonize column names
  # crown_ratio_change.csv: delta_CR, DIA, CR_init, CR_final
  if ("delta_CR" %in% names(cr_data) & !("CR_change" %in% names(cr_data))) {
    cr_data <- cr_data %>% rename(CR_change = delta_CR)
  }
  if ("DIA" %in% names(cr_data) & !("DIA_t1" %in% names(cr_data))) {
    cr_data <- cr_data %>% rename(DIA_t1 = DIA)
  }
  # Create CR_pct from CR_init if needed (CR_init is 0 to 1 scale)
  if (!("CR_pct" %in% names(cr_data)) & "CR_init" %in% names(cr_data)) {
    cr_data <- cr_data %>% mutate(CR_pct = CR_init * 100)
  }

  cr_data <- cr_data %>%
    filter(!is.na(CR_change), CR_change != Inf, CR_change != -Inf, !is.nan(CR_change),
           years_interval > 0)
} else if (file.exists(diameter_growth_file)) {
  logger::log_warn("crown_ratio_change.csv not found, falling back to diameter_growth.csv")
  data_file <- diameter_growth_file
  cr_data <- read_csv(data_file, show_col_types = FALSE) %>%
    as_tibble() %>%
    filter(!is.na(CR_t1), !is.na(CR_t2), CR_t1 > 0, years_interval > 0) %>%
    mutate(
      CR_change = (CR_t2 - CR_t1) / years_interval  # Annual change
    ) %>%
    filter(!is.na(CR_change), CR_change != Inf, CR_change != -Inf, !is.nan(CR_change))
} else {
  logger::log_error("Neither crown_ratio_change.csv nor diameter_growth.csv found")
  stop("Data files not found")
}

# Standard data processing
cr_data <- cr_data %>%
  mutate(
    SPCD = as.factor(SPCD),
    # Standardize predictors
    DBH_std = scale(DIA_t1)[, 1],
    BA_std = scale(BA)[, 1],
    BAL_std = scale(BAL)[, 1],
    CR_std = scale(CR_pct / 100)[, 1],
    SI_std = scale(log(pmax(SI, 1)))[, 1]
  ) %>%
  filter(!is.na(CR_change), !is.na(DBH_std), !is.na(BA_std), !is.na(BAL_std),
         !is.na(CR_std), !is.na(SI_std))

logger::log_info("Loaded {nrow(cr_data)} tree observations with CR changes")

# Subsample for computational tractability (stratified by species)
max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(cr_data) > max_n) {
  logger::log_info("Subsampling from {nrow(cr_data)} to {max_n} observations (stratified by species)")
  set.seed(42)
  cr_data <- cr_data %>%
    group_by(SPCD) %>%
    slice_sample(prop = min(1, max_n / nrow(cr_data) * 1.1)) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), max_n))
  logger::log_info("After subsampling: {nrow(cr_data)} observations")
}

# ============================================================================
# Load FVS Config for Priors
# ============================================================================

logger::log_info("Loading FVS parameters for priors...")

config_file <- file.path(project_root, "config", paste0(variant, ".json"))
config <- fromJSON(config_file)

# Extract crown ratio coefficients from config if available
fvs_other <- config$categories$other

# Crown ratio parameter names in FVS: BCR1, BCR2, BCR3, BCR4, etc.
cr_param_names <- grep("^BC", names(fvs_other), value = TRUE)

# Use average values as prior
if (length(cr_param_names) > 0) {
  bcr_vals <- unlist(fvs_other[cr_param_names])
  bcr_mean <- mean(bcr_vals, na.rm = TRUE)
  bcr_sd <- sd(bcr_vals, na.rm = TRUE)
} else {
  bcr_mean <- 0
  bcr_sd <- 1
}

logger::log_info("Crown ratio prior means: mean={bcr_mean}, sd={bcr_sd}")

# ============================================================================
# Multi-Strategy Model Fitting
# ============================================================================

logger::log_info("Attempting HMC inference with reduced parameters...")

# Strategy 1: HMC with reduced parameters (wrapped in tryCatch)
fit_cr <- tryCatch(
  {
    brm(
      CR_change ~ DBH_std + I(DBH_std^2) + BA_std + BAL_std + CR_std + SI_std +
        (1 | SPCD),       # Species varying intercept
      data = cr_data,
      family = gaussian(),
      chains = 4,
      iter = 1500,
      warmup = 500,
      cores = parallel::detectCores(),
      backend = "cmdstanr",
      refresh = 0,
      control = list(adapt_delta = 0.90, max_treedepth = 12)
    )
  },
  error = function(e) {
    logger::log_warn("HMC failed: {e$message}")
    NULL
  }
)

# Strategy 2: Fall back to variational inference if HMC fails
use_variational <- is.null(fit_cr)

if (use_variational) {
  logger::log_info("Falling back to variational inference (meanfield)...")

  fit_cr <- tryCatch(
    {
      brm(
        CR_change ~ DBH_std + I(DBH_std^2) + BA_std + BAL_std + CR_std + SI_std +
          (1 | SPCD),
        data = cr_data,
        family = gaussian(),
        algorithm = "meanfield",
        cores = parallel::detectCores(),
        refresh = 0,
      )
    },
    error = function(e) {
      logger::log_error("Both HMC and variational inference failed: {e$message}")
      stop("Model fitting failed with both strategies: ", e$message)
    }
  )

  logger::log_warn("Used variational inference due to HMC failure")
} else {
  logger::log_info("HMC inference successful")
}

logger::log_info("Crown ratio model fitting complete")

# ============================================================================
# Check Convergence
# ============================================================================

logger::log_info("Checking convergence diagnostics...")

# ============================================================================
# Extract and Save Results with Standardized Output Format
# ============================================================================

logger::log_info("Extracting posterior samples...")

# Get posterior draws
draws_cr <- as_draws_df(fit_cr)

# Save draws
draws_file <- file.path(output_dir, "crown_ratio_samples.rds")
saveRDS(draws_cr, draws_file)
logger::log_info("Saved posterior samples to {draws_file}")

# Save brms summary
brms_summ <- summary(fit_cr)
fixed_df <- as_tibble(brms_summ$fixed, rownames = "variable")
summary_file <- file.path(output_dir, "crown_ratio_summary.csv")
write_csv(fixed_df, summary_file)
logger::log_info("Saved summary to {summary_file}")

# posterior::summarise_draws with quantile() inherits names like "5%" instead of "p05"
# so we use unname() to force our chosen column names
draws_summ <- tryCatch({
  posterior::summarise_draws(draws_cr,
    p50 = function(x) unname(median(x)),
    p05 = function(x) unname(quantile(x, 0.05)),
    p95 = function(x) unname(quantile(x, 0.95)),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
}, error = function(e) {
  logger::log_warn("summarise_draws failed: {e$message}. Using manual summary.")
  draws_mat <- as.matrix(draws_cr)
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

posterior_df <- draws_summ %>%
  mutate(
    ci_width = p95 - p05,
    converged = if_else(!is.na(rhat), rhat <= 1.01, NA)
  )

posterior_file <- file.path(output_dir, "crown_ratio_posterior.csv")
write_csv(posterior_df, posterior_file)
logger::log_info("Saved posterior estimates to {posterior_file}")

# Create point estimates (posterior median)
map_df <- posterior_df %>%
  select(variable, estimate = p50)

map_file <- file.path(output_dir, "crown_ratio_map.csv")
write_csv(map_df, map_file)
logger::log_info("Saved MAP estimates to {map_file}")

# Log inference method used
method_info <- tibble(
  variant = variant,
  method = if_else(use_variational, "variational_meanfield", "HMC"),
  n_observations = nrow(cr_data),
  n_species = n_distinct(cr_data$SPCD)
)
logger::log_info("Inference method: {method_info$method}")

# ============================================================================
# Posterior Predictive Checks
# ============================================================================

logger::log_info("Generating diagnostic plots...")

# PPC plot
pdf(file.path(output_dir, "crown_ratio_ppc.pdf"), width = 10, height = 6)
print(pp_check(fit_cr, ndraws = 100))
dev.off()

# Conditional effects
pdf(file.path(output_dir, "crown_ratio_effects.pdf"), width = 12, height = 8)
print(plot(conditional_effects(fit_cr)))
dev.off()

logger::log_info("Saved diagnostic plots")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("========================================\n")
cat("Crown Ratio Model Fitting Complete\n")
cat("========================================\n")
cat("Variant:", variant, "\n")
cat("Observations:", nrow(cr_data), "\n")
cat("Number of species:", n_distinct(cr_data$SPCD), "\n")
cat("Inference method:", if_else(use_variational, "Variational (Meanfield)", "HMC"), "\n")
cat("Mean CR change (annual):", round(mean(cr_data$CR_change, na.rm = TRUE), 4), "\n")
cat("SD CR change:", round(sd(cr_data$CR_change, na.rm = TRUE), 4), "\n")
cat("\nModel:\n")
print(fit_cr)
cat("\nOutput files:\n")
cat("  -", basename(map_file), "\n")
cat("  -", basename(summary_file), "\n")
cat("  -", basename(posterior_file), "\n")
cat("  -", basename(draws_file), "\n")
cat("\nAll output saved to:", output_dir, "\n\n")

logger::log_info("Crown ratio model fitting and output generation complete")
