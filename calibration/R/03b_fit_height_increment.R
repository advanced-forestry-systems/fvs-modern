#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: Fit Height Increment Model
# For variants with explicit HG parameters (IE, CI, KT, BC, WS, EM)
#
# Variants WITHOUT HG parameters derive height growth from the H-D
# relationship and diameter growth prediction, so script 03
# (height-diameter) effectively calibrates their height growth.
#
# Model: ln(HTG) = species intercept + species*ln(DBH) + ln(HT) +
#         CR + SI + BAL + BA + SLOPE + cos(ASP) + species_group*ln(HT)
#         + plot random effect
#
# Multi-strategy inference:
#   1. MAP optimization (always runs first)
#   2. HMC sampling with MAP initialization (if feasible)
#   3. Variational inference fallback (if HMC fails)
#
# Usage: Rscript calibration/R/03b_fit_height_increment.R --variant ie
# =============================================================================

library(tidyverse)
library(data.table)
library(cmdstanr)
library(jsonlite)
library(posterior)
library(bayesplot)
library(tidybayes)
library(parallel)
library(logger)

# =============================================================================
# Parse Command Line Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
variant <- "ie"  # Default to IE which has HG params

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variant <- args[i + 1]
    }
  }
}

# =============================================================================
# Configuration
# =============================================================================

project_root <- Sys.getenv("FVS_PROJECT_ROOT",
                             normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE))
calibration_dir <- file.path(project_root, "calibration")
processed_data_dir <- file.path(calibration_dir, "data", "processed")
output_dir <- file.path(calibration_dir, "output", "variants", variant)
stan_dir <- file.path(calibration_dir, "stan")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Logging
log_file <- file.path(calibration_dir, "logs",
                       paste0("03b_fit_htinc_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting height increment model for variant {variant}")

# =============================================================================
# Load FVS Config (for priors if available)
# =============================================================================

config_file <- file.path(project_root, "config", paste0(variant, ".json"))
has_config <- file.exists(config_file)
config <- if (has_config) fromJSON(config_file) else NULL

# Check for explicit HG parameter keys in the config (used as priors)
hg_param_names <- if (has_config && !is.null(config$categories$other)) {
  grep("^HG", names(config$categories$other), value = TRUE)
} else {
  character(0)
}

if (length(hg_param_names) > 0) {
  logger::log_info("Found {length(hg_param_names)} HG parameter arrays: {paste(hg_param_names, collapse = ', ')}")
  logger::log_info("Using FVS HG parameters as informative priors")
} else {
  logger::log_info("No explicit HG parameters in config; using weakly informative priors")
  logger::log_info("Height increment model will be fit from FIA remeasurement data directly")
}

# =============================================================================
# Load Data
# =============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

# Try dedicated height_growth.csv first, fall back to diameter_growth.csv
hg_file <- file.path(processed_data_dir, variant, "height_growth.csv")
dg_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (file.exists(hg_file)) {
  logger::log_info("Using dedicated height_growth.csv file")
  data_file <- hg_file
} else if (file.exists(dg_file)) {
  logger::log_info("height_growth.csv not found, falling back to diameter_growth.csv")
  data_file <- dg_file
} else {
  logger::log_error("No data file found for height increment fitting")
  stop("Data file not found")
}

raw_data <- read_csv(data_file, show_col_types = FALSE) %>%
  as_tibble()

# Harmonize column names
# height_growth.csv: DIA_t1, HT_t1, HT_t2, HT_annual, CR_pct
# diameter_growth.csv: DIA_t1, HT_t1, HT_t2, CR_pct (same names, good)

# Compute height increment from remeasurement pairs
# Need both HT_t1 and HT_t2 to calculate HTG
hg_data <- raw_data %>%
  filter(
    !is.na(HT_t1), !is.na(HT_t2),
    HT_t1 > 4.5, HT_t2 > 4.5,        # Above breast height
    !is.na(DIA_t1), DIA_t1 > 0,
    !is.na(years_interval), years_interval > 0,
    !is.na(CR_pct), CR_pct > 0, CR_pct <= 100,  # Filter NA/invalid crown ratio
    !is.na(SI), SI > 0,                           # Site index present
    !is.na(BAL),                                   # BAL present
    !is.na(BA), BA > 0,                            # Stand basal area present
    !is.na(SLOPE)                                  # Slope present
  ) %>%
  mutate(
    # Annual height increment
    HTG_annual = (HT_t2 - HT_t1) / years_interval,
    # Keep only positive growth (exclude measurement errors and dying trees)
    # Small negative values can occur from measurement error
    HTG_annual = pmax(HTG_annual, 0.01)
  ) %>%
  filter(
    HTG_annual > 0,
    HTG_annual < 10  # Sanity cap: no tree grows >10 ft/yr (or m/yr)
  ) %>%
  mutate(
    ln_HTG = log(HTG_annual),
    ln_DBH = log(DIA_t1),
    ln_HT = log(HT_t1),
    CR_prop = CR_pct / 100,         # Convert to proportion
    SI_scaled = SI / 100,           # Scale for numerical stability
    BAL_scaled = BAL / 100,
    BA_scaled = BA / 100,
    SLOPE_prop = SLOPE / 100,       # Scale slope percent
    CASP = if ("SLOPE_CASP" %in% names(raw_data)) {
      ifelse(SLOPE > 0, SLOPE_CASP / SLOPE, 0)
    } else {
      0  # Default when aspect data not available
    },
    SPCD = as.integer(SPCD)
  )

logger::log_info("Prepared {nrow(hg_data)} height increment observations")
logger::log_info("Mean annual height increment: {round(mean(hg_data$HTG_annual), 2)}")
logger::log_info("Height increment range: {round(min(hg_data$HTG_annual), 2)} to {round(max(hg_data$HTG_annual), 2)}")

# =============================================================================
# Create Species Group Mapping
# =============================================================================

# FVS uses species groups for the height growth curve shape parameters
# The number of groups matches the HGHC array length (or default to ~5 groups)
fvs_other <- if (has_config && !is.null(config$categories$other)) config$categories$other else list()
n_spgrp <- if (!is.null(fvs_other$HGHC)) length(fvs_other$HGHC) else min(5, length(unique(hg_data$SPCD)))

# Map species to species groups
# FVS typically uses ISPSPE or similar mapping; we approximate with ordinal groups
species_list <- sort(unique(hg_data$SPCD))
N_species <- length(species_list)

# Create species group assignment (cycle through groups if more species than groups)
species_to_group <- tibble(
  SPCD = species_list,
  species_idx = seq_along(species_list),
  spgrp_idx = ((seq_along(species_list) - 1) %% n_spgrp) + 1
)

hg_data <- hg_data %>%
  left_join(species_to_group, by = "SPCD")

# Subsample for computational tractability (stratified by species)
max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(hg_data) > max_n) {
  logger::log_info("Subsampling from {nrow(hg_data)} to {max_n} observations (stratified by species)")
  set.seed(42)
  subsample_prop <- min(1, max_n / nrow(hg_data) * 1.1)
  hg_data <- hg_data %>%
    group_by(species_idx) %>%
    slice_sample(prop = subsample_prop) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), max_n))
  # Recompute consecutive indices after subsampling
  hg_data <- hg_data %>%
    mutate(
      species_idx = as.integer(factor(species_idx)),
      spgrp_idx = as.integer(factor(spgrp_idx))
    )
  N_species <- max(hg_data$species_idx)
  n_spgrp <- max(hg_data$spgrp_idx)
  logger::log_info("After subsampling: {nrow(hg_data)} observations")
}

logger::log_info("Species: {N_species}, Species groups: {n_spgrp}")

# Persist the dense species_idx -> SPCD crosswalk so downstream tooling
# (calibration/R/multipliers.R) can map b0[i] to the correct FVS species without
# replaying the stochastic, subsample-dependent data prep.
hi_species_index <- hg_data %>%
  dplyr::distinct(species_idx, SPCD) %>%
  dplyr::arrange(species_idx)
utils::write.csv(hi_species_index,
                 file.path(output_dir, "height_increment_species_index.csv"),
                 row.names = FALSE)

# =============================================================================
# Extract Priors from FVS Config
# =============================================================================

logger::log_info("Extracting priors from FVS config...")

# fvs_other was already set above when determining n_spgrp

# HGLD: ln(DBH) coefficient by species (N_species length in FVS)
prior_hgld <- if ("HGLD" %in% names(fvs_other)) {
  vals <- unlist(fvs_other$HGLD)
  # Pad or truncate to match our species count
  rep_len(vals, N_species)
} else {
  rep(0, N_species)
}

# HGHC: height growth curve shape by species group
prior_hghc <- if ("HGHC" %in% names(fvs_other)) {
  rep_len(unlist(fvs_other$HGHC), n_spgrp)
} else {
  rep(1, n_spgrp)
}

# HGLDD: ln(DBH) coefficient by species group
prior_hgldd <- if ("HGLDD" %in% names(fvs_other)) {
  rep_len(unlist(fvs_other$HGLDD), n_spgrp)
} else {
  rep(0.5, n_spgrp)
}

# HGH2: HT^2 coefficient by species group
prior_hgh2 <- if ("HGH2" %in% names(fvs_other)) {
  rep_len(unlist(fvs_other$HGH2), n_spgrp)
} else {
  rep(0, n_spgrp)
}

logger::log_info("Priors extracted: HGLD mean={round(mean(prior_hgld), 4)}, HGHC mean={round(mean(prior_hghc), 4)}")

# =============================================================================
# Prepare Stan Data
# =============================================================================

logger::log_info("Preparing data for Stan model...")

N <- nrow(hg_data)

stan_data <- list(
  N = N,
  N_species = N_species,
  N_spgrp = n_spgrp,

  # Response
  ln_HTG = hg_data$ln_HTG,

  # Predictors
  ln_DBH = hg_data$ln_DBH,
  ln_HT = hg_data$ln_HT,
  CR = hg_data$CR_prop,
  SI = hg_data$SI_scaled,
  BAL = hg_data$BAL_scaled,
  BA = hg_data$BA_scaled,
  SLOPE = hg_data$SLOPE_prop,
  CASP = hg_data$CASP,

  # Group indices
  species_id = hg_data$species_idx,
  spgrp_id = hg_data$spgrp_idx,

  # Priors from FVS config
  prior_hgld = prior_hgld,
  prior_hghc = prior_hghc,
  prior_hgldd = prior_hgldd,
  prior_hgh2 = prior_hgh2
)

logger::log_info("Stan data: N={N}, N_species={N_species}, N_spgrp={n_spgrp}")

# =============================================================================
# Compile Stan Model
# =============================================================================

logger::log_info("Compiling Stan model...")

stan_file <- file.path(stan_dir, "height_increment.stan")

tryCatch({
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    dir = output_dir
  )
  logger::log_info("Stan model compiled successfully")
}, error = function(e) {
  logger::log_error("Failed to compile Stan model: {e$message}")
  stop(e)
})

# =============================================================================
# Fit Model: Multi-Strategy Inference
# =============================================================================

logger::log_info("Fitting Bayesian height increment model (multi-strategy)...")

# ---- Stage 1: MAP optimization (fast, always works) ----
logger::log_info("Stage 1: MAP optimization...")

fit_map <- tryCatch({
  mod$optimize(
    data = stan_data,
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
  # Wrap summary() in tryCatch: optimization may return an object that

  # cannot produce draws (e.g., line search failure)
  map_estimates <- tryCatch({
    fit_map$summary()
  }, error = function(e) {
    logger::log_warn("MAP summary extraction failed: {e$message}")
    NULL
  })

  if (is.null(map_estimates)) {
    fit_map <- NULL  # Mark as failed so downstream uses fallback
  }
}

if (!is.null(fit_map) && !is.null(map_estimates)) {
  logger::log_info("MAP optimization converged. Log posterior: {round(fit_map$lp(), 2)}")

  # Save MAP estimates
  map_file <- file.path(output_dir, "height_increment_map.csv")
  write_csv(map_estimates, map_file)
  logger::log_info("Saved MAP estimates to {map_file}")

  # Build init from MAP for subsequent sampling
  map_vals <- setNames(map_estimates$estimate, map_estimates$variable)
  init_from_map <- function() {
    # Extract MAP values for initialization
    # b0 and b1 are vectors[N_species], gamma_shape is vector[N_spgrp]
    b0_names <- paste0("b0[", seq_len(N_species), "]")
    b1_names <- paste0("b1[", seq_len(N_species), "]")
    gs_names <- paste0("gamma_shape[", seq_len(n_spgrp), "]")
    list(
      b0 = unname(map_vals[b0_names]),
      b1 = unname(map_vals[b1_names]),
      b2 = unname(map_vals["b2"]),
      b3 = unname(map_vals["b3"]),
      b4 = unname(map_vals["b4"]),
      b5 = unname(map_vals["b5"]),
      b6 = unname(map_vals["b6"]),
      b7 = unname(map_vals["b7"]),
      b8 = unname(map_vals["b8"]),
      gamma_shape = unname(map_vals[gs_names]),
      mu_b0 = unname(map_vals["mu_b0"]),
      tau_b0 = max(unname(map_vals["tau_b0"]), 0.01),
      mu_b1 = unname(map_vals["mu_b1"]),
      tau_b1 = max(unname(map_vals["tau_b1"]), 0.01),
      sigma = max(unname(map_vals["sigma"]), 0.01)
    )
  }
} else {
  init_from_map <- NULL
}

# ---- Stage 2: Try HMC sampling with MAP initialization ----
logger::log_info("Stage 2: HMC sampling with MAP initialization...")

sampling_method <- "none"  # track which method produced the final fit

fit <- tryCatch({
  result <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = min(4, parallel::detectCores()),
    iter_warmup = 500,
    iter_sampling = 1000,
    adapt_delta = 0.90,
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

# =============================================================================
# Check Convergence
# =============================================================================

logger::log_info("Computing convergence diagnostics (method: {sampling_method})...")

if (sampling_method == "hmc") {
  diagnostics <- fit$diagnostic_summary()
  logger::log_info("HMC diagnostics: divergent={sum(diagnostics$num_divergent)}, treedepth={sum(diagnostics$num_max_treedepth)}")
}

# Get summary statistics
summ <- tryCatch({
  if (sampling_method == "map_only") {
    # Use MAP estimates as summary
    if (!is.null(map_estimates) && nrow(map_estimates) > 0) {
      map_estimates %>%
        rename(median = estimate) %>%
        mutate(q5 = NA_real_, q95 = NA_real_, rhat = NA_real_, ess_bulk = NA_real_)
    } else {
      logger::log_warn("No MAP estimates available; writing empty summary")
      tibble(variable = character(), median = numeric(),
             q5 = numeric(), q95 = numeric(), rhat = numeric(), ess_bulk = numeric())
    }
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
  if (!is.null(map_estimates) && nrow(map_estimates) > 0) {
    map_estimates %>%
      rename(median = estimate) %>%
      mutate(q5 = NA_real_, q95 = NA_real_, rhat = NA_real_, ess_bulk = NA_real_)
  } else {
    tibble(variable = character(), median = numeric(),
           q5 = numeric(), q95 = numeric(), rhat = numeric(), ess_bulk = numeric())
  }
})

# Check Rhat and ESS (only meaningful for HMC)
if (sampling_method == "hmc" && "rhat" %in% names(summ)) {
  rhat_issues <- summ %>%
    filter(rhat > 1.01 & !is.na(rhat))

  if (nrow(rhat_issues) > 0) {
    logger::log_warn("Found {nrow(rhat_issues)} parameters with Rhat > 1.01")
    logger::log_info("Parameters: {paste(rhat_issues$variable, collapse = ', ')}")
  } else {
    logger::log_info("All parameters converged (Rhat <= 1.01)")
  }
}

# =============================================================================
# Extract and Save Posterior Samples
# =============================================================================

logger::log_info("Extracting posterior samples (method: {sampling_method})...")

if (sampling_method != "map_only") {
  # Get draws in tibble format
  draws_df <- tryCatch(
    fit$draws(format = "df"),
    error = function(e) { logger::log_warn("Draw extraction failed: {e$message}"); NULL }
  )

  if (!is.null(draws_df)) {
    draws_file <- file.path(output_dir, "height_increment_samples.rds")
    saveRDS(draws_df, draws_file)
    logger::log_info("Saved posterior samples to {draws_file}")
  }
} else {
  draws_df <- NULL
  logger::log_info("MAP only: no posterior samples to save")
}

# Save method used
writeLines(sampling_method, file.path(output_dir, "height_increment_method.txt"))

# Save summary
summary_file <- file.path(output_dir, "height_increment_summary.csv")
write_csv(summ, summary_file)
logger::log_info("Saved summary statistics to {summary_file}")

# =============================================================================
# Extract Posterior Medians and Credible Intervals
# =============================================================================

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
posterior_file <- file.path(output_dir, "height_increment_posterior.csv")
write_csv(posterior_summary, posterior_file)
logger::log_info("Saved posterior summaries to {posterior_file}")

# =============================================================================
# Diagnostic Plots
# =============================================================================

logger::log_info("Generating diagnostic plots...")

# Traceplots (HMC only)
if (sampling_method == "hmc") {
  key_params <- c("b0[1]", "b1", "b2", "b3", "b4", "b5", "sigma",
                   "gamma_shape[1]")
  pdf(file.path(output_dir, "height_increment_traceplots.pdf"), width = 12, height = 8)
  tryCatch({
    print(bayesplot::mcmc_trace(
      fit$draws(),
      pars = intersect(key_params, summ$variable),
      facet_args = list(nrow = 3)
    ))
  }, error = function(e) {
    logger::log_warn("Traceplot generation failed: {e$message}")
  })
  dev.off()
}

# Posterior predictive check
if (sampling_method != "map_only") {
  tryCatch({
    y_rep <- fit$draws("ln_HTG_rep", format = "matrix")
    ppc_file <- file.path(output_dir, "height_increment_ppc.pdf")
    pdf(ppc_file, width = 10, height = 6)
    print(bayesplot::ppc_dens_overlay(
      y = stan_data$ln_HTG,
      yrep = y_rep[1:100, ]
    ))
    dev.off()
    logger::log_info("Saved posterior predictive check")
  }, error = function(e) {
    logger::log_warn("PPC plot failed: {e$message}")
  })
}

# =============================================================================
# Summary Report
# =============================================================================

cat("\n")
cat("==========================================\n")
cat("Height Increment Model Fitting Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("Inference method:", sampling_method, "\n")
cat("Species fitted:", N_species, "\n")
cat("Species groups:", n_spgrp, "\n")
cat("Total observations:", N, "\n")
cat("Mean annual height increment:", round(mean(hg_data$HTG_annual), 2), "\n")
cat("\nHG parameters in config:", paste(hg_param_names, collapse = ", "), "\n")
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Height increment model fitting complete for variant {variant} (method: {sampling_method})")
