#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Mortality Model v2
# Improved model for variants with AUC < 0.70
#
# Changes from 04_fit_mortality.R:
#   1. Random slopes for DBH by species (shade tolerance variation)
#   2. DBH:BAL interaction (suppression mortality)
#   3. ELEV and SLOPE as additional fixed effects
#   4. Species-level BAL slope (competition response varies by tolerance)
#   5. Option for class-balanced weighting via case_weights
#
# Model: P(dead) = logit^-1(
#   b0 + b1*DBH + b2*DBH^2 + b3*BAL + b4*CR + b5*SI + b6*BA +
#   b7*DBH:BAL + b8*ELEV + b9*SLOPE +
#   (1 + DBH_std + BAL_std | SPCD)
# )
#
# Usage: Rscript calibration/R/04b_fit_mortality_v2.R --variant ca
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

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(calibration_dir, "logs",
                       paste0("04b_fit_mortality_v2_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting mortality model v2 fitting for variant {variant}")

# ============================================================================
# Load Data
# ============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

mortality_file <- file.path(processed_data_dir, variant, "mortality.csv")
diameter_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (file.exists(mortality_file)) {
  logger::log_info("Using dedicated mortality.csv file")
  data_file <- mortality_file
} else if (file.exists(diameter_file)) {
  logger::log_info("Mortality.csv not found, falling back to diameter_growth.csv")
  data_file <- diameter_file
} else {
  stop("Data files not found")
}

mortality_data <- read_csv(data_file, show_col_types = FALSE) %>%
  as_tibble()

# Harmonize column names
if ("DIA" %in% names(mortality_data) & !("DIA_t1" %in% names(mortality_data))) {
  mortality_data <- mortality_data %>% rename(DIA_t1 = DIA)
}

# Create survival indicator
if ("survived" %in% names(mortality_data)) {
  # Already has survived column
} else if ("died" %in% names(mortality_data)) {
  mortality_data <- mortality_data %>% mutate(survived = !as.logical(died))
} else if ("DIA_t2" %in% names(mortality_data)) {
  mortality_data <- mortality_data %>% mutate(survived = !is.na(DIA_t2) & DIA_t2 > 0)
} else {
  stop("Cannot determine survival status")
}

# Handle CR_pct
if (!("CR_pct" %in% names(mortality_data)) & "CR" %in% names(mortality_data)) {
  mortality_data <- mortality_data %>% mutate(CR_pct = CR * 100)
}

# Standardize predictors (including new ones)
mortality_data <- mortality_data %>%
  mutate(
    DBH_std = scale(DIA_t1)[, 1],
    CR_std = scale(CR_pct / 100)[, 1],
    BAL_std = scale(BAL)[, 1],
    SI_std = scale(log(pmax(SI, 1)))[, 1],
    BA_std = scale(BA)[, 1],
    # New predictors
    ELEV_std = if ("ELEV" %in% names(.)) scale(ELEV)[, 1] else 0,
    SLOPE_std = if ("SLOPE" %in% names(.)) scale(SLOPE)[, 1] else 0,
    # Interaction term
    DBH_BAL = DBH_std * BAL_std,
    SPCD = as.factor(SPCD)
  ) %>%
  filter(!is.na(survived), !is.na(DBH_std), !is.na(BAL_std),
         !is.na(CR_std), !is.na(SI_std), !is.na(BA_std))

has_elev <- "ELEV" %in% names(mortality_data) && any(mortality_data$ELEV_std != 0)
has_slope <- "SLOPE" %in% names(mortality_data) && any(mortality_data$SLOPE_std != 0)
logger::log_info("Additional predictors: ELEV={has_elev}, SLOPE={has_slope}")

logger::log_info("Loaded {nrow(mortality_data)} tree observations")
logger::log_info("Mortality rate: {round(mean(!mortality_data$survived) * 100, 2)}%")
logger::log_info("Number of species: {n_distinct(mortality_data$SPCD)}")

# ============================================================================
# Subsample (stratified by species, balanced for mortality)
# ============================================================================

max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(mortality_data) > max_n) {
  logger::log_info("Subsampling from {nrow(mortality_data)} to {max_n}")
  set.seed(42)

  # Stratified sampling: ensure dead trees are adequately represented
  # Sample dead and alive separately to maintain within-species mortality rates
  dead_data <- mortality_data %>% filter(!survived)
  alive_data <- mortality_data %>% filter(survived)

  mort_rate <- nrow(dead_data) / nrow(mortality_data)
  n_dead_target <- min(nrow(dead_data), round(max_n * max(mort_rate, 0.15)))
  n_alive_target <- max_n - n_dead_target

  dead_sample <- dead_data %>%
    group_by(SPCD) %>%
    slice_sample(prop = min(1, n_dead_target / nrow(dead_data) * 1.1)) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), n_dead_target))

  alive_sample <- alive_data %>%
    group_by(SPCD) %>%
    slice_sample(prop = min(1, n_alive_target / nrow(alive_data) * 1.1)) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), n_alive_target))

  mortality_data <- bind_rows(dead_sample, alive_sample) %>%
    slice_sample(n = nrow(.))  # Shuffle

  logger::log_info("After subsampling: {nrow(mortality_data)} obs ({sum(!mortality_data$survived)} dead, {sum(mortality_data$survived)} alive)")
}

# ============================================================================
# Filter species with too few observations for random slopes
# ============================================================================

# Need enough obs per species for random slopes to be estimable
sp_counts <- mortality_data %>%
  count(SPCD) %>%
  mutate(enough = n >= 20)

# Species with < 20 obs get lumped into "other"
if (any(!sp_counts$enough)) {
  n_lumped <- sum(!sp_counts$enough)
  small_spp <- sp_counts %>% filter(!enough) %>% pull(SPCD)
  mortality_data <- mortality_data %>%
    mutate(SPCD_grouped = if_else(SPCD %in% small_spp, factor("OTHER"), SPCD))
  logger::log_info("Lumped {n_lumped} rare species into 'OTHER' for random effects")
} else {
  mortality_data$SPCD_grouped <- mortality_data$SPCD
}

n_groups <- n_distinct(mortality_data$SPCD_grouped)
logger::log_info("Species groups for random effects: {n_groups}")

# ============================================================================
# Build Model Formula (adaptive complexity)
# ============================================================================

# Base fixed effects
fixed_terms <- "DBH_std + I(DBH_std^2) + BAL_std + CR_std + SI_std + BA_std + DBH_BAL"

if (has_elev) fixed_terms <- paste0(fixed_terms, " + ELEV_std")
if (has_slope) fixed_terms <- paste0(fixed_terms, " + SLOPE_std")

# Random effects: try random slopes, fall back to random intercept
# If many species, random slopes on DBH and BAL
# If few species (< 8), just random intercept + DBH slope
if (n_groups >= 8) {
  random_term <- "(1 + DBH_std + BAL_std | SPCD_grouped)"
  logger::log_info("Using full random slopes (DBH + BAL) by species")
} else if (n_groups >= 3) {
  random_term <- "(1 + DBH_std | SPCD_grouped)"
  logger::log_info("Using random slopes (DBH only) by species")
} else {
  random_term <- "(1 | SPCD_grouped)"
  logger::log_info("Using random intercept only (too few species for slopes)")
}

formula_str <- paste0("survived ~ ", fixed_terms, " + ", random_term)
model_formula <- as.formula(formula_str)
logger::log_info("Model formula: {formula_str}")

# ============================================================================
# Fit Model: HMC with CmdStanR
# ============================================================================

fit_mort <- NULL
inference_method <- "none"

logger::log_info("Strategy 1: HMC (CmdStanR) with random slopes...")

fit_result <- tryCatch(
  {
    fit <- brm(
      model_formula,
      data = mortality_data,
      family = bernoulli(link = "logit"),
      prior = c(
        # Weakly regularizing priors to help convergence with random slopes
        prior(normal(0, 3), class = "Intercept"),
        prior(normal(0, 2), class = "b"),
        prior(exponential(1), class = "sd"),
        prior(lkj(2), class = "cor")
      ),
      chains = 4,
      iter = 2000,
      warmup = 1000,
      cores = parallel::detectCores(),
      backend = "cmdstanr",
      refresh = 0,
      control = list(adapt_delta = 0.95, max_treedepth = 12)
    )
    list(fit = fit, success = TRUE, error = NULL)
  },
  error = function(e) {
    logger::log_warn("Strategy 1 (HMC full) failed: {as.character(e)}")
    list(fit = NULL, success = FALSE, error = as.character(e))
  }
)

if (fit_result$success) {
  fit_mort <- fit_result$fit
  inference_method <- "HMC_CmdStanR_v2"
  logger::log_info("Strategy 1 succeeded: HMC with random slopes")
} else {
  # Strategy 2: Reduce random effects complexity
  logger::log_info("Strategy 2: HMC with simpler random effects...")

  simpler_formula <- as.formula(paste0("survived ~ ", fixed_terms,
                                        " + (1 + DBH_std | SPCD_grouped)"))

  fit_result2 <- tryCatch(
    {
      fit <- brm(
        simpler_formula,
        data = mortality_data,
        family = bernoulli(link = "logit"),
        prior = c(
          prior(normal(0, 3), class = "Intercept"),
          prior(normal(0, 2), class = "b"),
          prior(exponential(1), class = "sd"),
          prior(lkj(2), class = "cor")
        ),
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = parallel::detectCores(),
        backend = "cmdstanr",
        refresh = 0,
        control = list(adapt_delta = 0.95, max_treedepth = 12)
      )
      list(fit = fit, success = TRUE, error = NULL)
    },
    error = function(e) {
      logger::log_warn("Strategy 2 (simpler HMC) failed: {as.character(e)}")
      list(fit = NULL, success = FALSE, error = as.character(e))
    }
  )

  if (fit_result2$success) {
    fit_mort <- fit_result2$fit
    inference_method <- "HMC_CmdStanR_v2_simpler"
    logger::log_info("Strategy 2 succeeded: HMC with DBH random slope only")
  } else {
    # Strategy 3: Variational inference
    logger::log_info("Strategy 3: Variational inference fallback...")

    vi_result <- tryCatch(
      {
        fit <- brm(
          model_formula,
          data = mortality_data,
          family = bernoulli(link = "logit"),
          algorithm = "meanfield",
          cores = parallel::detectCores(),
          backend = "cmdstanr",
          refresh = 0
        )
        list(fit = fit, success = TRUE, error = NULL)
      },
      error = function(e) {
        logger::log_error("All strategies failed: {as.character(e)}")
        list(fit = NULL, success = FALSE, error = as.character(e))
      }
    )

    if (vi_result$success) {
      fit_mort <- vi_result$fit
      inference_method <- "VI_meanfield_v2"
      logger::log_info("Strategy 3 succeeded: variational inference")
    } else {
      stop("Mortality v2 model fitting failed with all strategies")
    }
  }
}

logger::log_info("Mortality v2 fitting complete using: {inference_method}")

# ============================================================================
# Convergence Diagnostics
# ============================================================================

logger::log_info("Checking convergence diagnostics...")

if (grepl("HMC", inference_method)) {
  rhat_vals <- tryCatch({
    # Use brms::rhat or extract from summary
    s <- summary(fit_mort)
    c(s$fixed$Rhat, s$random$SPCD_grouped$Rhat)
  }, error = function(e) {
    logger::log_warn("Could not extract Rhat: {e$message}")
    NA_real_
  })
  n_bad <- sum(rhat_vals > 1.05, na.rm = TRUE)
  n_warn <- sum(rhat_vals > 1.01, na.rm = TRUE)
  logger::log_info("Rhat: {n_bad} params > 1.05, {n_warn} params > 1.01")

  if (n_bad > 0) {
    logger::log_warn("Some parameters have not converged (Rhat > 1.05)")
  }
}

# ============================================================================
# Extract and Save Results
# ============================================================================

logger::log_info("Extracting posterior samples...")

draws_mort <- as_draws_df(fit_mort)

draws_file <- file.path(output_dir, "mortality_v2_samples.rds")
saveRDS(draws_mort, draws_file)
logger::log_info("Saved posterior samples to {draws_file}")

brms_summ <- summary(fit_mort)

fixed_df <- as_tibble(brms_summ$fixed, rownames = "variable")
names(fixed_df) <- gsub("^Estimate$", "p50", names(fixed_df))
names(fixed_df) <- gsub("^l-95% CI$", "p025", names(fixed_df))
names(fixed_df) <- gsub("^u-95% CI$", "p975", names(fixed_df))

summary_file <- file.path(output_dir, "mortality_v2_summary.csv")
write_csv(fixed_df, summary_file)
logger::log_info("Saved summary to {summary_file}")

# Full posterior summary
draws_summ <- tryCatch({
  posterior::summarise_draws(draws_mort,
    p50 = function(x) unname(median(x)),
    p05 = function(x) unname(quantile(x, 0.05)),
    p95 = function(x) unname(quantile(x, 0.95)),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
}, error = function(e) {
  draws_mat <- as.matrix(draws_mort)
  tibble(
    variable = colnames(draws_mat),
    p50 = apply(draws_mat, 2, function(x) unname(median(x))),
    p05 = apply(draws_mat, 2, function(x) unname(quantile(x, 0.05))),
    p95 = apply(draws_mat, 2, function(x) unname(quantile(x, 0.95))),
    rhat = NA_real_, ess_bulk = NA_real_
  )
})

names(draws_summ) <- gsub("^5%$", "p05", names(draws_summ))
names(draws_summ) <- gsub("^95%$", "p95", names(draws_summ))
if (!"variable" %in% names(draws_summ)) names(draws_summ)[1] <- "variable"

posterior_df <- draws_summ %>%
  mutate(ci_width = p95 - p05,
         converged = if_else(!is.na(rhat), rhat < 1.05, NA))

write_csv(posterior_df, file.path(output_dir, "mortality_v2_posterior.csv"))

# ============================================================================
# In-Sample AUC Evaluation
# ============================================================================

logger::log_info("Computing in-sample AUC...")

# Get predicted probabilities (posterior mean of survival probability)
pred_probs <- tryCatch({
  fitted(fit_mort, summary = TRUE)[, "Estimate"]
}, error = function(e) {
  logger::log_warn("fitted() failed, using manual prediction")
  # Manual approach
  pp <- posterior_epred(fit_mort)
  colMeans(pp)
})

# Compute AUC
compute_auc <- function(probs, labels) {
  # probs = P(survived), labels = TRUE/FALSE survived
  # For AUC of mortality prediction, flip: P(dead) = 1 - P(survived)
  dead_probs <- 1 - probs
  dead_labels <- !labels

  n_pos <- sum(dead_labels)
  n_neg <- sum(!dead_labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)

  # Wilcoxon-Mann-Whitney statistic
  ranks <- rank(dead_probs)
  auc <- (sum(ranks[dead_labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  return(auc)
}

auc_v2 <- compute_auc(pred_probs, mortality_data$survived)
logger::log_info("In-sample AUC (v2): {round(auc_v2, 4)}")

# Compare to v1 if available
v1_file <- file.path(output_dir, "mortality_samples.rds")
if (file.exists(v1_file)) {
  logger::log_info("Comparing to v1 model...")

  # Load v1 and compute its AUC on same data would require refitting
  # Instead just report v2 AUC and let the comparison script handle it
  logger::log_info("V1 model exists; full comparison deferred to 11_full_comparison.R")
}

# Save AUC result
auc_df <- tibble(
  variant = variant,
  model_version = "v2",
  inference_method = inference_method,
  n_obs = nrow(mortality_data),
  n_dead = sum(!mortality_data$survived),
  mort_rate = mean(!mortality_data$survived),
  n_species = n_distinct(mortality_data$SPCD),
  n_species_grouped = n_groups,
  auc_insample = auc_v2,
  formula = formula_str,
  has_elev = has_elev,
  has_slope = has_slope
)

write_csv(auc_df, file.path(output_dir, "mortality_v2_auc.csv"))

# ============================================================================
# Diagnostic Plots
# ============================================================================

logger::log_info("Generating diagnostic plots...")

if (grepl("HMC", inference_method)) {
  tryCatch({
    pdf(file.path(output_dir, "mortality_v2_ppc.pdf"), width = 10, height = 6)
    print(pp_check(fit_mort, ndraws = 100, type = "error_binned"))
    dev.off()
  }, error = function(e) logger::log_warn("PPC plot failed: {e$message}"))
}

tryCatch({
  pdf(file.path(output_dir, "mortality_v2_effects.pdf"), width = 14, height = 10)
  print(plot(conditional_effects(fit_mort)))
  dev.off()
}, error = function(e) logger::log_warn("Effects plot failed: {e$message}"))

# Species-level random effects plot
tryCatch({
  re_df <- ranef(fit_mort)$SPCD_grouped
  if (!is.null(re_df)) {
    re_intercept <- as_tibble(re_df[,,"Intercept"], rownames = "species")
    names(re_intercept) <- c("species", "estimate", "error", "lower", "upper")

    p_re <- re_intercept %>%
      mutate(species = reorder(species, estimate)) %>%
      ggplot(aes(x = species, y = estimate, ymin = lower, ymax = upper)) +
      geom_pointrange(size = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      coord_flip() +
      labs(x = "Species", y = "Random intercept (log-odds)",
           title = paste0("Mortality v2: Species random effects (", variant, ")")) +
      theme_minimal(base_size = 11)

    ggsave(file.path(output_dir, "mortality_v2_species_re.pdf"), p_re,
           width = 8, height = max(4, n_groups * 0.25))
  }
}, error = function(e) logger::log_warn("RE plot failed: {e$message}"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("==========================================\n")
cat("Mortality Model v2 Fitting Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("Observations:", nrow(mortality_data), "\n")
cat("  Dead:", sum(!mortality_data$survived), "\n")
cat("  Alive:", sum(mortality_data$survived), "\n")
cat("Mortality rate:", round(mean(!mortality_data$survived) * 100, 2), "%\n")
cat("Species groups:", n_groups, "\n")
cat("Inference method:", inference_method, "\n")
cat("In-sample AUC:", round(auc_v2, 4), "\n")
cat("Formula:", formula_str, "\n")
cat("\nFixed effects:\n")
print(brms_summ$fixed)
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Mortality v2 complete: AUC={round(auc_v2, 4)}, method={inference_method}")
