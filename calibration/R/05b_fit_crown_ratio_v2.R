#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: Fit Crown Ratio Model v2
# Improved model: predict CR level (time 2) instead of CR change
#
# Changes from 05_fit_crown_ratio.R:
#   1. Response is CR_t2 (crown ratio at remeasurement) instead of delta_CR
#   2. CR_t1 included as strong predictor (autoregressive structure)
#   3. Beta regression family (CR is bounded 0-1)
#   4. Random slopes for DBH by species
#   5. Additional predictors: relative DBH within stand
#
# Model: CR_t2 ~ beta(mu, phi)
#   logit(mu) = b0 + b1*CR_t1 + b2*DBH + b3*DBH^2 + b4*BAL + b5*BA +
#               b6*SI + b7*relDBH + (1 + DBH_std | SPCD)
#
# Usage: Rscript calibration/R/05b_fit_crown_ratio_v2.R --variant ca
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
                       paste0("05b_fit_crown_ratio_v2_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting crown ratio model v2 fitting for variant {variant}")

# ============================================================================
# Load Data
# ============================================================================

logger::log_info("Loading FIA data for variant {variant}...")

# Need time 1 and time 2 CR values, so use diameter_growth.csv or
# crown_ratio_change.csv which has both
cr_change_file <- file.path(processed_data_dir, variant, "crown_ratio_change.csv")
diameter_growth_file <- file.path(processed_data_dir, variant, "diameter_growth.csv")

if (file.exists(cr_change_file)) {
  logger::log_info("Using crown_ratio_change.csv")
  cr_data <- read_csv(cr_change_file, show_col_types = FALSE) %>% as_tibble()

  # Harmonize columns
  if ("DIA" %in% names(cr_data) & !("DIA_t1" %in% names(cr_data))) {
    cr_data <- cr_data %>% rename(DIA_t1 = DIA)
  }
  # Need CR at time 1 and time 2
  if ("CR_init" %in% names(cr_data) & "CR_final" %in% names(cr_data)) {
    cr_data <- cr_data %>%
      rename(CR_t1_raw = CR_init, CR_t2_raw = CR_final)
  }
} else if (file.exists(diameter_growth_file)) {
  logger::log_info("Using diameter_growth.csv (has CR_t1, CR_t2)")
  cr_data <- read_csv(diameter_growth_file, show_col_types = FALSE) %>% as_tibble()
  # diameter_growth.csv has CR_t1 and CR_t2 (already on 0-1 or 0-100 scale)
  if ("CR_t1" %in% names(cr_data) & "CR_t2" %in% names(cr_data)) {
    cr_data <- cr_data %>%
      rename(CR_t1_raw = CR_t1, CR_t2_raw = CR_t2)
  }
} else {
  stop("Data files not found")
}

# Ensure CR is on 0 to 1 scale
cr_data <- cr_data %>%
  mutate(
    CR_t1 = if_else(CR_t1_raw > 1, CR_t1_raw / 100, CR_t1_raw),
    CR_t2 = if_else(CR_t2_raw > 1, CR_t2_raw / 100, CR_t2_raw)
  ) %>%
  # Clamp to (0, 1) exclusive for beta regression
  mutate(
    CR_t1 = pmin(pmax(CR_t1, 0.005), 0.995),
    CR_t2 = pmin(pmax(CR_t2, 0.005), 0.995)
  ) %>%
  filter(!is.na(CR_t1), !is.na(CR_t2),
         !is.na(DIA_t1), DIA_t1 > 0,
         !is.na(BAL), !is.na(BA), !is.na(SI))

# Compute relative DBH (tree size relative to stand average)
cr_data <- cr_data %>%
  group_by(PLT_CN) %>%
  mutate(
    mean_dia_plot = mean(DIA_t1, na.rm = TRUE),
    rel_DBH = DIA_t1 / mean_dia_plot
  ) %>%
  ungroup()

# Standardize predictors
cr_data <- cr_data %>%
  mutate(
    SPCD = as.factor(SPCD),
    CR_t1_std = scale(CR_t1)[, 1],
    DBH_std = scale(DIA_t1)[, 1],
    BA_std = scale(BA)[, 1],
    BAL_std = scale(BAL)[, 1],
    SI_std = scale(log(pmax(SI, 1)))[, 1],
    relDBH_std = scale(rel_DBH)[, 1]
  ) %>%
  filter(!is.na(CR_t1_std), !is.na(relDBH_std))

logger::log_info("Loaded {nrow(cr_data)} tree observations with paired CR measurements")
logger::log_info("CR_t1 range: {round(min(cr_data$CR_t1), 3)} to {round(max(cr_data$CR_t1), 3)}")
logger::log_info("CR_t2 range: {round(min(cr_data$CR_t2), 3)} to {round(max(cr_data$CR_t2), 3)}")

# ============================================================================
# Subsample
# ============================================================================

max_n <- as.integer(Sys.getenv("FVS_MAX_OBS", "30000"))
if (nrow(cr_data) > max_n) {
  logger::log_info("Subsampling from {nrow(cr_data)} to {max_n}")
  set.seed(42)
  cr_data <- cr_data %>%
    group_by(SPCD) %>%
    slice_sample(prop = min(1, max_n / nrow(cr_data) * 1.1)) %>%
    ungroup() %>%
    slice_sample(n = min(nrow(.), max_n))
  logger::log_info("After subsampling: {nrow(cr_data)} observations")
}

# Lump rare species
sp_counts <- cr_data %>% count(SPCD) %>% mutate(enough = n >= 20)
if (any(!sp_counts$enough)) {
  small_spp <- sp_counts %>% filter(!enough) %>% pull(SPCD)
  cr_data <- cr_data %>%
    mutate(SPCD_grouped = if_else(SPCD %in% small_spp, factor("OTHER"), SPCD))
  logger::log_info("Lumped {length(small_spp)} rare species into OTHER")
} else {
  cr_data$SPCD_grouped <- cr_data$SPCD
}

n_groups <- n_distinct(cr_data$SPCD_grouped)
logger::log_info("Species groups: {n_groups}")

# ============================================================================
# Build Formula
# ============================================================================

fixed_terms <- "CR_t1_std + DBH_std + I(DBH_std^2) + BA_std + BAL_std + SI_std + relDBH_std"

if (n_groups >= 8) {
  random_term <- "(1 + DBH_std | SPCD_grouped)"
} else if (n_groups >= 3) {
  random_term <- "(1 | SPCD_grouped)"
} else {
  random_term <- "1"  # No random effects
}

formula_str <- paste0("CR_t2 ~ ", fixed_terms, " + ", random_term)
model_formula <- as.formula(formula_str)
logger::log_info("Model formula: {formula_str}")

# ============================================================================
# Fit Model: Beta regression via brms
# ============================================================================

fit_cr <- NULL
inference_method <- "none"

# Strategy 1: HMC with beta family
logger::log_info("Strategy 1: HMC beta regression...")

fit_result <- tryCatch(
  {
    fit <- brm(
      model_formula,
      data = cr_data,
      family = Beta(link = "logit", link_phi = "log"),
      prior = c(
        prior(normal(0, 3), class = "Intercept"),
        prior(normal(0, 2), class = "b"),
        prior(exponential(1), class = "sd")
      ),
      chains = 4,
      iter = 2000,
      warmup = 1000,
      cores = parallel::detectCores(),
      backend = "cmdstanr",
      refresh = 0,
      control = list(adapt_delta = 0.95, max_treedepth = 12)
    )
    list(fit = fit, success = TRUE)
  },
  error = function(e) {
    logger::log_warn("Strategy 1 (HMC beta) failed: {as.character(e)}")
    list(fit = NULL, success = FALSE)
  }
)

if (fit_result$success) {
  fit_cr <- fit_result$fit
  inference_method <- "HMC_beta_v2"
  logger::log_info("Strategy 1 succeeded: HMC beta regression")
} else {
  # Strategy 2: Gaussian on logit-transformed CR
  logger::log_info("Strategy 2: HMC Gaussian on logit(CR_t2)...")

  cr_data$logit_CR_t2 <- qlogis(cr_data$CR_t2)
  formula_gauss <- as.formula(gsub("CR_t2", "logit_CR_t2", formula_str))

  fit_result2 <- tryCatch(
    {
      fit <- brm(
        formula_gauss,
        data = cr_data,
        family = gaussian(),
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = parallel::detectCores(),
        backend = "cmdstanr",
        refresh = 0,
        control = list(adapt_delta = 0.90)
      )
      list(fit = fit, success = TRUE)
    },
    error = function(e) {
      logger::log_warn("Strategy 2 (Gaussian logit) failed: {as.character(e)}")
      list(fit = NULL, success = FALSE)
    }
  )

  if (fit_result2$success) {
    fit_cr <- fit_result2$fit
    inference_method <- "HMC_gaussian_logit_v2"
    logger::log_info("Strategy 2 succeeded: Gaussian on logit(CR)")
  } else {
    # Strategy 3: Variational beta
    logger::log_info("Strategy 3: Variational beta regression...")

    vi_result <- tryCatch(
      {
        fit <- brm(
          model_formula,
          data = cr_data,
          family = Beta(link = "logit", link_phi = "log"),
          algorithm = "meanfield",
          cores = parallel::detectCores(),
          backend = "cmdstanr",
          refresh = 0
        )
        list(fit = fit, success = TRUE)
      },
      error = function(e) {
        logger::log_error("All strategies failed: {as.character(e)}")
        list(fit = NULL, success = FALSE)
      }
    )

    if (vi_result$success) {
      fit_cr <- vi_result$fit
      inference_method <- "VI_beta_v2"
    } else {
      stop("Crown ratio v2 model fitting failed with all strategies")
    }
  }
}

logger::log_info("Crown ratio v2 fitting complete: {inference_method}")

# ============================================================================
# Extract and Save Results
# ============================================================================

logger::log_info("Extracting posterior...")

draws_cr <- as_draws_df(fit_cr)
saveRDS(draws_cr, file.path(output_dir, "crown_ratio_v2_samples.rds"))

brms_summ <- summary(fit_cr)
fixed_df <- as_tibble(brms_summ$fixed, rownames = "variable")
write_csv(fixed_df, file.path(output_dir, "crown_ratio_v2_summary.csv"))

draws_summ <- tryCatch({
  posterior::summarise_draws(draws_cr,
    p50 = function(x) unname(median(x)),
    p05 = function(x) unname(quantile(x, 0.05)),
    p95 = function(x) unname(quantile(x, 0.95)),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
}, error = function(e) {
  draws_mat <- as.matrix(draws_cr)
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

write_csv(draws_summ, file.path(output_dir, "crown_ratio_v2_posterior.csv"))

# ============================================================================
# In-Sample R² Evaluation
# ============================================================================

logger::log_info("Computing in-sample R²...")

# Get predicted CR_t2
pred_cr <- tryCatch({
  fitted(fit_cr, summary = TRUE)[, "Estimate"]
}, error = function(e) {
  pp <- posterior_epred(fit_cr)
  colMeans(pp)
})

# If we used logit-transformed Gaussian, back-transform
if (inference_method == "HMC_gaussian_logit_v2") {
  pred_cr <- plogis(pred_cr)
}

# R² on original scale
ss_res <- sum((cr_data$CR_t2 - pred_cr)^2)
ss_tot <- sum((cr_data$CR_t2 - mean(cr_data$CR_t2))^2)
r2_v2 <- 1 - ss_res / ss_tot

# RMSE
rmse_v2 <- sqrt(mean((cr_data$CR_t2 - pred_cr)^2))

# Bias
bias_v2 <- mean(pred_cr - cr_data$CR_t2)

logger::log_info("In-sample R² (v2, level): {round(r2_v2, 4)}")
logger::log_info("In-sample RMSE (v2): {round(rmse_v2, 4)}")
logger::log_info("In-sample bias (v2): {round(bias_v2, 5)}")

# Compare to v1 if available (v1 predicted change, so R² not directly comparable)
# But we can report both
eval_df <- tibble(
  variant = variant,
  model_version = "v2_level",
  inference_method = inference_method,
  n_obs = nrow(cr_data),
  n_species = n_distinct(cr_data$SPCD),
  n_species_grouped = n_groups,
  response = "CR_t2",
  family = if_else(grepl("beta", inference_method, ignore.case = TRUE),
                   "beta", "gaussian_logit"),
  r2_insample = r2_v2,
  rmse_insample = rmse_v2,
  bias_insample = bias_v2,
  formula = formula_str
)

write_csv(eval_df, file.path(output_dir, "crown_ratio_v2_eval.csv"))

# ============================================================================
# Diagnostic Plots
# ============================================================================

logger::log_info("Generating diagnostic plots...")

# Observed vs predicted
p_obspred <- tibble(
  observed = cr_data$CR_t2,
  predicted = pred_cr
) %>%
  ggplot(aes(x = observed, y = predicted)) +
  geom_hex(bins = 50) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    x = "Observed CR (time 2)",
    y = "Predicted CR (time 2)",
    title = paste0("Crown Ratio v2: Observed vs Predicted (", variant, ")"),
    subtitle = paste0("R² = ", round(r2_v2, 3),
                      ", RMSE = ", round(rmse_v2, 3),
                      ", Bias = ", round(bias_v2, 4))
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "crown_ratio_v2_obspred.pdf"), p_obspred,
       width = 8, height = 7)

tryCatch({
  pdf(file.path(output_dir, "crown_ratio_v2_effects.pdf"), width = 14, height = 10)
  print(plot(conditional_effects(fit_cr)))
  dev.off()
}, error = function(e) logger::log_warn("Effects plot failed: {e$message}"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat("==========================================\n")
cat("Crown Ratio Model v2 Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("Observations:", nrow(cr_data), "\n")
cat("Species groups:", n_groups, "\n")
cat("Inference method:", inference_method, "\n")
cat("Response: CR level (time 2), not change\n")
cat("Family:", if_else(grepl("beta", inference_method, ignore.case = TRUE),
                       "Beta", "Gaussian (logit-transformed)"), "\n")
cat("In-sample R²:", round(r2_v2, 4), "\n")
cat("In-sample RMSE:", round(rmse_v2, 4), "\n")
cat("In-sample bias:", round(bias_v2, 5), "\n")
cat("Formula:", formula_str, "\n")
cat("\nFixed effects:\n")
print(brms_summ$fixed)
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Crown ratio v2 complete: R2={round(r2_v2, 4)}, method={inference_method}")
