#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Height-Diameter Model
# Chapman-Richards equation: H = 4.5 + a * (1 - exp(-b * DBH))^c
# Multi-strategy inference: MAP optimization → HMC sampling → Variational inference (ADVI)
#
# Usage: Rscript calibration/R/03_fit_height_diameter.R --variant ca
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
library(parallel)

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
log_file <- file.path(calibration_dir, "logs", paste0("03_fit_htdbh_", variant, ".log"))
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting height-diameter model fitting for variant {variant} (multi-strategy inference)")

# ============================================================================
# Load Data
# ============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

# Try dedicated height_diameter.csv first (larger dataset), fall back to diameter_growth.csv
hd_file <- file.path(processed_data_dir, variant, "height_diameter.csv")
dg_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (file.exists(hd_file)) {
  logger::log_info("Using dedicated height_diameter.csv file")
  data_file <- hd_file
} else if (file.exists(dg_file)) {
  logger::log_info("height_diameter.csv not found, falling back to diameter_growth.csv")
  data_file <- dg_file
} else {
  logger::log_error("No data file found for height diameter fitting")
  stop("Data files not found")
}

raw_data <- read_csv(data_file, show_col_types = FALSE) %>%
  as_tibble()

# Harmonize column names across file formats
# height_diameter.csv: DIA, HT
# diameter_growth.csv: DIA_t1, HT_t1
if ("DIA" %in% names(raw_data) & !("DIA_t1" %in% names(raw_data))) {
  raw_data <- raw_data %>% rename(DIA_t1 = DIA)
}
if ("HT" %in% names(raw_data) & !("HT_t1" %in% names(raw_data))) {
  raw_data <- raw_data %>% rename(HT_t1 = HT)
}

dg_data <- raw_data %>%
  # Select unique tree measurements (take the first measurement if multiple)
  {
    if ("tree_id" %in% names(.)) {
      group_by(., tree_id) %>% slice(1) %>% ungroup()
    } else {
      .
    }
  } %>%
  # Remove missing heights
  filter(!is.na(HT_t1), HT_t1 > 4.5, !is.na(DIA_t1), DIA_t1 > 0) %>%
  mutate(
    SPCD = as.factor(SPCD),
    # Scale predictors for better model fitting
    DBH_scaled = DIA_t1 / 20,  # Center around typical diameter
    SI_scaled = pmax(SI, 1) / 50  # Center around typical SI
  )

logger::log_info("Loaded {nrow(dg_data)} tree-level observations")

# Subsample for computational tractability (stratified by species)
max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(dg_data) > max_n) {
  logger::log_info("Subsampling from {nrow(dg_data)} to {max_n} observations (stratified by species)")
  set.seed(42)
  dg_data <- dg_data %>%
    group_by(SPCD) %>%
    slice_sample(prop = min(1, max_n / nrow(dg_data) * 1.1)) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), max_n))
  logger::log_info("After subsampling: {nrow(dg_data)} observations")
}

n_obs <- nrow(dg_data)
n_species <- n_distinct(dg_data$SPCD)

# ============================================================================
# Load FVS Config for Priors
# ============================================================================

logger::log_info("Loading FVS parameters for priors...")

config_file <- file.path(project_root, "config", paste0(variant, ".json"))
config <- fromJSON(config_file)

# Chapman Richards parameters: H = b0 + b1 * (1 - exp(-b2 * DBH))^b3
# In FVS these may be parameterized differently
# We'll use weakly informative priors based on typical values

# ============================================================================
# Fit Height-Diameter Model: Multi-Strategy Inference
# ============================================================================

logger::log_info("Fitting height-diameter model with multi-strategy inference...")
logger::log_info("Model: H = 4.5 + a * (1 - exp(-b * DBH_scaled))^c with species random effects")

# Chapman-Richards equation using brms nonlinear model
# H = 4.5 + a * (1 - exp(-b * DBH))^c

# Specify priors for nonlinear parameters (use nlpar, not coef, for brms nl models)
# Chapman-Richards requires a > 0, b > 0, c > 0 to avoid NaN from (1-exp(-b*x))^c
# Using lb (lower bound) to constrain parameter space
priors_htdbh <- c(
  prior(normal(25, 10), nlpar = "a", class = "b", lb = 1),    # Asymptotic height (must be > 0)
  prior(normal(0.05, 0.02), nlpar = "b", class = "b", lb = 0.001), # Rate (must be > 0)
  prior(normal(1, 0.5), nlpar = "c", class = "b", lb = 0.1),  # Shape (must be > 0)
  prior(exponential(1), class = "sigma")                       # Error scale
)

# Custom initialization from MAP estimates (if available)
# Random inits from (-2, 2) cause NaN because species random effects can push
# b negative -> (1-exp(-b*DBH))^c = NaN. Using MAP values ensures safe start.
init_fun <- function() {
  list(
    b_a = array(runif(1, 15, 35)),        # a ~ 15 to 35 ft asymptote
    b_b = array(runif(1, 0.02, 0.08)),    # b ~ small positive rate
    b_c = array(runif(1, 0.5, 1.5))       # c ~ shape near 1
  )
}

# Build init from MAP estimates (will be set after MAP stage completes)
build_map_inits <- function(map_fit, stan_data_list) {
  # Extract unconstrained MAP parameter values for all chains
  tryCatch({
    # CmdStanR optimize uses $mle() not $draws()
    map_params <- tryCatch(
      as.list(map_fit$mle()),
      error = function(e) {
        # Fallback: try summary
        s <- map_fit$summary()
        setNames(as.list(s$estimate), s$variable)
      }
    )
    # Remove lp__ if present
    map_params[["lp__"]] <- NULL
    # Return same init for all 4 chains (with small jitter)
    lapply(1:4, function(i) {
      jittered <- lapply(map_params, function(x) {
        if (is.numeric(x) && length(x) > 0) {
          x * runif(length(x), 0.95, 1.05)
        } else {
          x
        }
      })
      jittered
    })
  }, error = function(e) {
    logger::log_warn("Could not build MAP inits: {e$message}")
    NULL
  })
}

# Build the model formula once (reusable across strategies)
model_formula <- bf(
  HT_t1 ~ 4.5 + a * (1 - exp(-b * DBH_scaled))^c,
  a ~ 1 + (1 | SPCD),      # Species varying intercept
  b ~ 1 + (1 | SPCD),      # Species varying rate
  c ~ 1,                   # Shared shape parameter
  nl = TRUE
)

sampling_method <- "none"  # Track which method produces final fit
fit_htdbh <- NULL
fit_map <- NULL

# ---- Stage 1: MAP optimization (always fast) ----
logger::log_info("Stage 1: MAP optimization via brms...")

tryCatch({
  # Compile model first (chains=0 just compiles, doesn't sample)
  fit_compiled <- brm(
    model_formula,
    data = dg_data,
    family = gaussian(),
    prior = priors_htdbh,
    chains = 0,
    backend = "cmdstanr",
  )

  logger::log_info("Model compiled successfully")

  # Extract stan model and data for MAP optimization
  stan_model <- fit_compiled$fit
  stan_data_list <- make_standata(model_formula, data = dg_data, family = gaussian(), prior = priors_htdbh)

  # Run MAP optimization
  fit_map <- tryCatch({
    opt_result <- stan_model$optimize(
      data = stan_data_list,
      algorithm = "lbfgs",
      iter = 10000,
      tol_rel_grad = 1e-8,
      refresh = 100
    )
    opt_result
  }, error = function(e) {
    logger::log_warn("MAP optimization failed: {e$message}")
    NULL
  })

  if (!is.null(fit_map)) {
    map_estimates <- tryCatch(fit_map$summary(), error = function(e) NULL)
    lp_val <- tryCatch(fit_map$lp(), error = function(e) NA)
    logger::log_info("MAP optimization converged. Log posterior: {round(lp_val, 2)}")

    # Save MAP estimates
    if (!is.null(map_estimates)) {
      map_file <- file.path(output_dir, "height_diameter_map.csv")
      write_csv(map_estimates, map_file)
      logger::log_info("Saved MAP estimates to {map_file}")
    }
  }
}, error = function(e) {
  logger::log_warn("Model compilation or MAP failed: {e$message}")
})

# ---- Stage 2: HMC sampling (reduced iterations for faster convergence) ----
logger::log_info("Stage 2: HMC sampling with reduced iterations...")

# Use MAP estimates for init if available (most robust for nonlinear models)
hmc_init <- if (!is.null(fit_map)) {
  map_inits <- build_map_inits(fit_map, stan_data_list)
  if (!is.null(map_inits)) {
    logger::log_info("Using MAP estimates as HMC starting values")
    map_inits
  } else {
    logger::log_info("MAP init extraction failed, using init=0 (unconstrained zero)")
    0
  }
} else {
  logger::log_info("No MAP fit available, using init=0")
  0
}

fit_htdbh <- tryCatch({
  result <- brm(
    model_formula,
    data = dg_data,
    family = gaussian(),
    prior = priors_htdbh,
    chains = 4,
    iter = 1500,               # Reduced from 2000
    warmup = 500,              # Reduced from 1000
    cores = parallel::detectCores(),
    backend = "cmdstanr",
    refresh = 100,
    init = hmc_init,           # MAP-based or unconstrained-zero initialization
    control = list(adapt_delta = 0.90, max_treedepth = 12)
  )

  # Accept fit if sampling completed (brms validates internally)
  sampling_method <<- "hmc"
  logger::log_info("HMC sampling completed successfully")
  result
}, error = function(e) {
  logger::log_warn("HMC sampling failed: {e$message}")
  logger::log_info("Falling back to variational inference...")
  NULL
})

# ---- Stage 3: Variational inference fallback (ADVI) ----
if (is.null(fit_htdbh) || sampling_method == "none") {
  logger::log_info("Stage 3: Variational inference (ADVI)...")

  fit_htdbh <- tryCatch({
    brm(
      model_formula,
      data = dg_data,
      family = gaussian(),
      prior = priors_htdbh,
      chains = 1,
      iter = 50000,
      algorithm = "meanfield",
      backend = "cmdstanr",
      refresh = 500,
      init = 0  # Start VI from prior means (avoids NaN at random inits)
    )
  }, error = function(e) {
    logger::log_warn("Variational inference failed: {e$message}")
    NULL
  })

  if (!is.null(fit_htdbh)) {
    sampling_method <<- "variational"
    logger::log_info("Variational inference complete")
  } else {
    if (!is.null(fit_map)) {
      sampling_method <<- "map_only"
      logger::log_warn("All sampling methods failed; using MAP estimates only")
    } else {
      logger::log_error("All fitting methods failed; model fitting cannot continue")
      stop("Model fitting failed with all strategies")
    }
  }
}

logger::log_info("Height-diameter model fitting complete (method: {sampling_method})")

# ============================================================================
# Check Convergence Diagnostics
# ============================================================================

logger::log_info("Computing convergence diagnostics (method: {sampling_method})...")

# Get summary statistics
summ_htdbh <- tryCatch({
  if (sampling_method == "map_only") {
    # Use MAP estimates as summary
    map_estimates %>%
      rename(median = estimate) %>%
      mutate(q5 = NA_real_, q95 = NA_real_, Rhat = NA_real_, Bulk_ESS = NA_real_)
  } else {
    s <- posterior_summary(fit_htdbh)
    # Ensure required columns exist (variational may lack Rhat/Bulk_ESS)
    if (!"Rhat" %in% colnames(s)) s <- cbind(s, Rhat = NA_real_)
    if (!"Bulk_ESS" %in% colnames(s)) s <- cbind(s, Bulk_ESS = NA_real_)
    if (!"Q5" %in% colnames(s)) s <- cbind(s, Q5 = NA_real_)
    if (!"Q95" %in% colnames(s)) s <- cbind(s, Q95 = NA_real_)
    as_tibble(s, rownames = "variable")
  }
}, error = function(e) {
  logger::log_warn("Summary computation failed: {e$message}")
  if (!is.null(fit_map)) {
    map_estimates %>%
      rename(median = estimate) %>%
      mutate(q5 = NA_real_, q95 = NA_real_, Rhat = NA_real_, Bulk_ESS = NA_real_)
  } else {
    tibble(variable = character(), median = numeric())
  }
})

# Check Rhat and ESS (only meaningful for HMC)
if (sampling_method == "hmc" && "Rhat" %in% names(summ_htdbh)) {
  rhat_issues <- summ_htdbh %>%
    filter(Rhat > 1.01 & !is.na(Rhat))

  if (nrow(rhat_issues) > 0) {
    logger::log_warn("Found {nrow(rhat_issues)} parameters with Rhat > 1.01")
    logger::log_info("Parameters with poor convergence: {paste(head(rhat_issues$variable, 5), collapse = ', ')}")
  } else {
    logger::log_info("All parameters converged (Rhat <= 1.01)")
  }
}

# ============================================================================
# Extract and Save Posterior Samples
# ============================================================================

logger::log_info("Extracting posterior samples (method: {sampling_method})...")

if (sampling_method != "map_only") {
  tryCatch({
    # Get posterior draws
    draws_htdbh <- as_draws_df(fit_htdbh)

    # Save draws
    draws_file <- file.path(output_dir, "height_diameter_samples.rds")
    saveRDS(draws_htdbh, draws_file)
    logger::log_info("Saved posterior samples to {draws_file}")
  }, error = function(e) {
    logger::log_warn("Failed to extract draws: {e$message}")
  })
} else {
  logger::log_info("MAP only: no posterior samples to save")
}

# Save method used
writeLines(sampling_method, file.path(output_dir, "height_diameter_method.txt"))

# Save summary statistics
summary_file <- file.path(output_dir, "height_diameter_summary.csv")
write_csv(summ_htdbh, summary_file)
logger::log_info("Saved summary statistics to {summary_file}")

# ============================================================================
# Extract Posterior Medians and Credible Intervals
# ============================================================================

logger::log_info("Computing posterior summaries...")

# Standardize column names and compute posterior summaries
# brms uses: Q2.5, Q97.5 for quantiles; rstan uses: Q5, Q95
# posterior_summary() returns: Estimate, Est.Error, Q2.5, Q97.5, Rhat, Bulk_ESS, Tail_ESS
col_names <- tolower(gsub("\\.", "_", colnames(summ_htdbh)))

posterior_summary <- summ_htdbh %>%
  as_tibble(rownames = "variable") %>%
  # Map various quantile column names to standard p05, p50, p95
  mutate(
    p50 = case_when(
      "estimate" %in% names(.) ~ estimate,
      "median" %in% names(.) ~ median,
      "q50" %in% names(.) ~ q50,
      TRUE ~ NA_real_
    ),
    p05 = case_when(
      "q5" %in% names(.) ~ q5,
      "q2_5" %in% names(.) ~ q2_5,
      TRUE ~ NA_real_
    ),
    p95 = case_when(
      "q95" %in% names(.) ~ q95,
      "q97_5" %in% names(.) ~ q97_5,
      TRUE ~ NA_real_
    ),
    rhat = if_else("rhat" %in% names(.), rhat, NA_real_),
    ess_bulk = if_else("bulk_ess" %in% names(.), bulk_ess, NA_real_),
    ci_width = p95 - p05,
    converged = if_else(is.na(rhat), NA, rhat <= 1.01)
  ) %>%
  select(variable, p50, p05, p95, rhat, ess_bulk, ci_width, converged)

# Save posterior summary
posterior_file <- file.path(output_dir, "height_diameter_posterior.csv")
write_csv(posterior_summary, posterior_file)
logger::log_info("Saved posterior summaries to {posterior_file}")

# ============================================================================
# Generate Diagnostic Plots
# ============================================================================

logger::log_info("Generating diagnostic plots...")

# PPC plot (only for non-MAP fits)
if (sampling_method != "map_only") {
  tryCatch({
    ppc_file <- file.path(output_dir, "height_diameter_ppc.pdf")
    pdf(ppc_file, width = 10, height = 6)
    print(pp_check(fit_htdbh, ndraws = 100))
    dev.off()
    logger::log_info("Saved posterior predictive check to {ppc_file}")
  }, error = function(e) {
    logger::log_warn("PPC plot generation failed: {e$message}")
  })
}

# Conditional effects plots (only for non-MAP fits)
if (sampling_method != "map_only") {
  tryCatch({
    effects_file <- file.path(output_dir, "height_diameter_effects.pdf")
    pdf(effects_file, width = 12, height = 8)
    plot(conditional_effects(fit_htdbh, effects = "DBH_scaled:SPCD"))
    dev.off()
    logger::log_info("Saved conditional effects to {effects_file}")
  }, error = function(e) {
    logger::log_warn("Effects plot generation failed: {e$message}")
  })
}

logger::log_info("Diagnostic plots complete")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("========================================\n")
cat("Height-Diameter Model Fitting Complete\n")
cat("========================================\n")
cat("Variant:", variant, "\n")
cat("Inference method:", sampling_method, "\n")
cat("Species fitted:", n_species, "\n")
cat("Total observations:", n_obs, "\n")
cat("\nModel summary (first 20 parameters):\n")
tryCatch(print(head(posterior_summary, 20)), error = function(e) cat("(summary unavailable)\n"))
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Height-diameter model fitting complete for variant {variant} (method: {sampling_method})")
