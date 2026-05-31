#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: Stand-Level Density Parameters
#
# Calibrates SDIMAX, BAMAX, and the self-thinning exponent using FIA data.
# Combines two estimation approaches:
#
# 1. Boundary line analysis via quantile regression (tau = 0.97-0.99)
#    to estimate the maximum size-density relationship.
#    This is the standard approach (Zhang et al. 2005, Weiskittel et al. 2011).
#
# 2. Bayesian hierarchical model with species level partial pooling
#    for SDIMAX and the self-thinning slope, using FVS current values
#    as informative priors.
#
# The key biological model is Reineke's (1933) self-thinning rule:
#    ln(TPA) = a + b * ln(QMD)
# where b is approximately -1.605 (Reineke's universal constant) but
# varies by species (Weiskittel et al. 2009, Pretzsch & Biber 2005).
# At a reference QMD of 10 inches: SDI = TPA * (QMD/10)^(-b)
# SDIMAX = the maximum SDI a species can sustain.
#
# Usage: Rscript calibration/R/09_fit_stand_density.R --variant ne
# =============================================================================

library(tidyverse)
library(data.table)
library(quantreg)       # Quantile regression for boundary line
library(brms)
library(cmdstanr)
library(jsonlite)
library(posterior)
library(bayesplot)
library(logger)

# =============================================================================
# Parse Command Line Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
variant <- "ne"

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

project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (!nzchar(project_root)) {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
calibration_dir <- file.path(project_root, "calibration")
processed_dir <- file.path(calibration_dir, "data", "processed", variant)
output_dir <- file.path(calibration_dir, "output", "variants", variant)
config_dir <- file.path(project_root, "config")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(calibration_dir, "logs",
                       paste0("09_stand_density_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting stand-level density calibration for variant {variant}")

# =============================================================================
# Load Data
# =============================================================================

# Stand-level data
stand_file <- file.path(processed_dir, "stand_density.csv")
if (!file.exists(stand_file)) {
  logger::log_error("Stand data not found: {stand_file}")
  stop("Run 08_fetch_stand_data.R first")
}

stand_data <- read_csv(stand_file, show_col_types = FALSE) %>% as_tibble()
logger::log_info("Loaded {nrow(stand_data)} stand observations")

# Species-level SDI data
species_file <- file.path(processed_dir, "species_sdi.csv")
species_sdi <- if (file.exists(species_file)) {
  read_csv(species_file, show_col_types = FALSE) %>% as_tibble()
} else {
  NULL
}

# Load FVS config for priors
config_file <- file.path(config_dir, paste0(variant, ".json"))
config <- fromJSON(config_file)
n_species <- config$maxsp

# Identify which SDIMAX parameter this variant uses
sdi_param_name <- NULL
sdi_priors <- NULL
bamax_param_name <- NULL
bamax_priors <- NULL

for (cat_name in names(config$categories)) {
  cat_data <- config$categories[[cat_name]]
  for (k in names(cat_data)) {
    if (k %in% c("SDICON", "R5SDI", "R4SDI", "FMSDI")) {
      sdi_param_name <- k
      sdi_priors <- unlist(cat_data[[k]])
      logger::log_info("Found SDI parameter: {k} in category '{cat_name}' ({length(sdi_priors)} values)")
    }
    if (k %in% c("BAMAXA", "BAMAX1")) {
      bamax_param_name <- k
      bamax_priors <- unlist(cat_data[[k]])
      logger::log_info("Found BAMAX parameter: {k} in category '{cat_name}' ({length(bamax_priors)} values)")
    }
  }
}

# Get FIA species code mapping
fia_species <- unlist(config$categories$species_definitions$FIAJSP)

# Build species code to FVS index lookup
spcd_to_fvs_idx <- setNames(seq_along(fia_species), fia_species)

logger::log_info("SDI param: {sdi_param_name}, BAMAX param: {bamax_param_name}")

# =============================================================================
# Part 1: Quantile Regression Boundary Line for Self-Thinning
# =============================================================================

logger::log_info("=== Part 1: Quantile Regression Boundary Line ===")

# Prepare data for Reineke's relationship: ln(TPA) = a + b * ln(QMD)
# Use all stands but focus on high density stands for the boundary

stand_qr <- stand_data %>%
  filter(
    tpa > 50,               # Reasonable minimum density
    qmd > 3,                # Exclude very small diameter stands
    ba_ft2ac > 30,          # Established stands
    n_trees >= 10            # Adequate sample within plot
  ) %>%
  mutate(
    ln_tpa = log(tpa),
    ln_qmd = log(qmd),
    SPCD = leading_spcd
  )

logger::log_info("Quantile regression data: {nrow(stand_qr)} stands")

# --- Overall self-thinning boundary (all species pooled) ---

qr_taus <- c(0.95, 0.97, 0.99)
qr_results <- list()

for (tau in qr_taus) {
  qr_fit <- rq(ln_tpa ~ ln_qmd, tau = tau, data = stand_qr)
  qr_results[[as.character(tau)]] <- list(
    tau = tau,
    intercept = coef(qr_fit)[1],
    slope = coef(qr_fit)[2],
    sdimax_at_10 = exp(coef(qr_fit)[1] + coef(qr_fit)[2] * log(10))
  )
  logger::log_info("  tau={tau}: slope={round(coef(qr_fit)[2], 3)}, SDIMAX(10in)={round(qr_results[[as.character(tau)]]$sdimax_at_10)}")
}

# Use tau=0.97 as primary (standard practice)
primary_qr <- qr_results[["0.97"]]
overall_slope <- primary_qr$slope
overall_sdimax <- primary_qr$sdimax_at_10

logger::log_info("Overall boundary: slope={round(overall_slope, 3)}, SDIMAX={round(overall_sdimax)}")

# --- Species specific self-thinning boundaries ---

# Only fit for species with sufficient data
species_counts <- stand_qr %>%
  count(SPCD) %>%
  filter(n >= 30)  # Need at least 30 stands per species

logger::log_info("Species with >= 30 stands: {nrow(species_counts)}")

species_qr_results <- list()

for (sp in species_counts$SPCD) {
  sp_data <- stand_qr %>% filter(SPCD == sp)

  tryCatch({
    qr_sp <- rq(ln_tpa ~ ln_qmd, tau = 0.97, data = sp_data)

    species_qr_results[[as.character(sp)]] <- tibble(
      SPCD = sp,
      n_stands = nrow(sp_data),
      qr_intercept = coef(qr_sp)[1],
      qr_slope = coef(qr_sp)[2],
      qr_sdimax_10 = exp(coef(qr_sp)[1] + coef(qr_sp)[2] * log(10)),
      # Also compute using Reineke's exponent for comparison
      reineke_sdimax = max(sp_data$sdi) * 0.97  # Rough upper bound
    )
  }, error = function(e) {
    logger::log_warn("QR failed for species {sp}: {e$message}")
  })
}

species_qr_df <- bind_rows(species_qr_results)
logger::log_info("Fitted QR for {nrow(species_qr_df)} species")

# Save quantile regression results
write_csv(species_qr_df, file.path(output_dir, "species_qr_sdimax.csv"))

# =============================================================================
# Part 2: Bayesian Hierarchical SDIMAX Model
# =============================================================================

logger::log_info("=== Part 2: Bayesian Hierarchical SDIMAX Model ===")

# Bayesian approach: fit the self-thinning line as a hierarchical model
# with species-level intercept (-> SDIMAX) and shared slope
#
# Model: ln(TPA) ~ Normal(a[species] + b * ln(QMD), sigma)
# But we want the UPPER boundary, not the mean. Use a truncated
# approach or a two-part model.
#
# Approach: Use only high-density stands (SDI > 60% of observed max
# for each species) and fit a normal model. This approximates the
# boundary without requiring the full quantile regression in Stan.

# Filter to high density observations per species
stand_bayes <- stand_qr %>%
  group_by(SPCD) %>%
  mutate(
    sdi_rank = percent_rank(sdi),
    high_density = sdi_rank >= 0.60  # Top 40% of density observations
  ) %>%
  ungroup() %>%
  filter(high_density) %>%
  mutate(
    species_idx = as.integer(factor(SPCD)),
    plot_idx = as.integer(factor(PLT_CN))
  )

# Build FVS prior for each species
# Map species in data to FVS species index
stand_bayes <- stand_bayes %>%
  mutate(
    fvs_idx = spcd_to_fvs_idx[as.character(SPCD)]
  )

# Get prior SDIMAX values from config
if (!is.null(sdi_priors)) {
  prior_sdimax_species <- sapply(unique(stand_bayes$SPCD), function(sp) {
    idx <- spcd_to_fvs_idx[as.character(sp)]
    if (!is.na(idx) && idx <= length(sdi_priors)) {
      return(sdi_priors[idx])
    }
    return(mean(sdi_priors, na.rm = TRUE))
  })
  names(prior_sdimax_species) <- unique(stand_bayes$SPCD)
} else {
  # No SDI priors, use overall quantile regression estimate
  prior_sdimax_species <- setNames(
    rep(overall_sdimax, n_distinct(stand_bayes$SPCD)),
    unique(stand_bayes$SPCD)
  )
}

N_species_bayes <- n_distinct(stand_bayes$species_idx)

logger::log_info("Bayesian model: {nrow(stand_bayes)} high-density stands, {N_species_bayes} species")

# --- Fit with brms ---

# Create species-level prior SDIMAX as a variable in the data
# The intercept of ln(TPA) = a + b*ln(QMD) at QMD=10 gives SDIMAX
# ln(SDIMAX) = a + b*ln(10)
# So a = ln(SDIMAX) - b*ln(10)

# Center the predictor around reference diameter
stand_bayes$ln_qmd_centered <- stand_bayes$ln_qmd - log(10)

# Now intercept directly estimates ln(TPA) at QMD=10, i.e., ln(SDIMAX)

# Fit the hierarchical model
fit_sdi <- brm(
  ln_tpa ~ ln_qmd_centered + (1 + ln_qmd_centered | SPCD),
  data = stand_bayes,
  family = gaussian(),
  prior = c(
    # Intercept = ln(SDIMAX): typical SDIMAX is 300-800, so ln is 5.7-6.7
    prior(normal(6.2, 0.5), class = "Intercept"),
    # Slope: Reineke's is -1.605, usually between -1.3 and -2.0
    prior(normal(-1.605, 0.3), class = "b", coef = "ln_qmd_centered"),
    # Species variation in SDIMAX
    prior(exponential(2), class = "sd"),
    # Correlation between intercept and slope
    prior(lkj(2), class = "cor"),
    # Residual error
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  iter = 3000,
  warmup = 1000,
  cores = parallel::detectCores(),
  backend = "cmdstanr",
  refresh = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

logger::log_info("Bayesian SDIMAX model fitted")

# =============================================================================
# Extract Species-Level SDIMAX Posteriors
# =============================================================================

logger::log_info("Extracting species-level SDIMAX estimates...")

# Get species-level predictions at QMD=10 (ln_qmd_centered = 0)
# The species-level intercept IS ln(SDIMAX)

# Extract random effect estimates
ranef_sdi <- ranef(fit_sdi)$SPCD
fixef_sdi <- fixef(fit_sdi)

# Population-level (average) SDIMAX
pop_ln_sdimax <- fixef_sdi["Intercept", "Estimate"]
pop_slope <- fixef_sdi["ln_qmd_centered", "Estimate"]
pop_sdimax <- exp(pop_ln_sdimax)

logger::log_info("Population SDIMAX: {round(pop_sdimax)} (slope: {round(pop_slope, 3)})")

# Species-level SDIMAX
species_sdimax <- tibble(
  SPCD = as.integer(rownames(ranef_sdi[, , "Intercept"])),
  ln_sdimax_offset = ranef_sdi[, , "Intercept"][, "Estimate"],
  ln_sdimax = pop_ln_sdimax + ln_sdimax_offset,
  sdimax_bayes = exp(ln_sdimax),
  slope_offset = ranef_sdi[, , "ln_qmd_centered"][, "Estimate"],
  species_slope = pop_slope + slope_offset
)

# Merge with quantile regression results
species_sdimax <- species_sdimax %>%
  left_join(species_qr_df %>% select(SPCD, qr_sdimax_10, qr_slope, n_stands),
            by = "SPCD") %>%
  mutate(
    # Combine estimates: weighted average of QR and Bayesian
    # Give more weight to QR for species with lots of data
    w_qr = pmin(n_stands / 100, 0.7),  # Cap QR weight at 0.7
    w_bayes = 1 - w_qr,
    sdimax_combined = w_qr * qr_sdimax_10 + w_bayes * sdimax_bayes,
    slope_combined = w_qr * qr_slope + w_bayes * species_slope
  )

# Add FVS prior values for comparison
species_sdimax <- species_sdimax %>%
  mutate(
    fvs_sdimax_prior = sapply(SPCD, function(sp) {
      idx <- spcd_to_fvs_idx[as.character(sp)]
      if (!is.null(sdi_priors) && !is.na(idx) && idx <= length(sdi_priors)) {
        return(sdi_priors[idx])
      }
      return(NA_real_)
    })
  )

logger::log_info("Species SDIMAX range: {round(min(species_sdimax$sdimax_combined))} to {round(max(species_sdimax$sdimax_combined))}")

# Save species SDIMAX estimates
write_csv(species_sdimax, file.path(output_dir, "species_sdimax_calibrated.csv"))

# =============================================================================
# Part 3: BAMAX Calibration (for variants that use it)
# =============================================================================

if (!is.null(bamax_param_name)) {
  logger::log_info("=== Part 3: BAMAX Calibration ===")

  # For variants that use BAMAX instead of SDIMAX, calibrate
  # maximum basal area by species from FIA observed distributions

  # Use quantile regression at tau=0.97 for BA ~ species
  ba_qr_results <- stand_qr %>%
    group_by(SPCD) %>%
    filter(n() >= 20) %>%
    summarise(
      n_stands = n(),
      ba_p95 = quantile(ba_ft2ac, 0.95),
      ba_p97 = quantile(ba_ft2ac, 0.97),
      ba_p99 = quantile(ba_ft2ac, 0.99),
      ba_max_obs = max(ba_ft2ac),
      ba_mean = mean(ba_ft2ac),
      .groups = "drop"
    ) %>%
    mutate(
      # Use 97th percentile as BAMAX estimate
      bamax_calibrated = ba_p97,
      # Get FVS prior
      fvs_bamax_prior = sapply(SPCD, function(sp) {
        idx <- spcd_to_fvs_idx[as.character(sp)]
        if (!is.na(idx) && idx <= length(bamax_priors)) {
          return(bamax_priors[idx])
        }
        return(NA_real_)
      })
    )

  # Bayesian shrinkage: blend observed BAMAX with FVS prior
  ba_qr_results <- ba_qr_results %>%
    mutate(
      # Weight by sample size (more data = trust observation more)
      w_obs = pmin(n_stands / 100, 0.8),
      w_prior = 1 - w_obs,
      bamax_combined = ifelse(
        is.na(fvs_bamax_prior),
        bamax_calibrated,
        w_obs * bamax_calibrated + w_prior * fvs_bamax_prior
      )
    )

  write_csv(ba_qr_results, file.path(output_dir, "species_bamax_calibrated.csv"))
  logger::log_info("BAMAX calibrated for {nrow(ba_qr_results)} species")
}

# =============================================================================
# Part 4: Self-Thinning Exponent (Reineke's b) by Species
# =============================================================================

logger::log_info("=== Part 4: Self-Thinning Exponent by Species ===")

# The standard FVS assumption is b = -1.605 (Reineke 1933)
# But empirical evidence shows this varies by species
# (Pretzsch & Biber 2005, Weiskittel et al. 2009)
# FVS stores BNORML which contains species group exponents

# Extract species-level slopes from the Bayesian model
slope_summary <- species_sdimax %>%
  select(SPCD, n_stands,
         qr_slope, species_slope = species_slope, slope_combined) %>%
  mutate(
    # Convert to Reineke exponent (positive form used in SDI formula)
    reineke_b_qr = -qr_slope,
    reineke_b_bayes = -species_slope,
    reineke_b_combined = -slope_combined,
    # Compare to standard -1.605
    diff_from_reineke = reineke_b_combined - 1.605
  )

write_csv(slope_summary, file.path(output_dir, "self_thinning_slopes.csv"))
logger::log_info("Self-thinning slopes saved for {nrow(slope_summary)} species")

# =============================================================================
# Part 5: Diagnostic Figures
# =============================================================================

logger::log_info("Generating diagnostic figures...")

# Publication ready theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )

# --- Figure 1: Self-thinning boundary with data ---
p_boundary <- ggplot(stand_qr, aes(x = ln_qmd, y = ln_tpa)) +
  geom_point(alpha = 0.15, size = 0.5, color = "grey40") +
  # Quantile regression lines
  geom_abline(
    intercept = qr_results[["0.95"]]$intercept,
    slope = qr_results[["0.95"]]$slope,
    color = "#0072B2", linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    intercept = qr_results[["0.97"]]$intercept,
    slope = qr_results[["0.97"]]$slope,
    color = "#D55E00", linewidth = 1.2
  ) +
  geom_abline(
    intercept = qr_results[["0.99"]]$intercept,
    slope = qr_results[["0.99"]]$slope,
    color = "#009E73", linetype = "dotted", linewidth = 0.8
  ) +
  # Reineke's reference line
  geom_abline(
    intercept = log(overall_sdimax) + 1.605 * log(10),
    slope = -1.605,
    color = "black", linetype = "longdash", linewidth = 0.6
  ) +
  annotate("text", x = max(stand_qr$ln_qmd) - 0.3, y = max(stand_qr$ln_tpa) - 0.5,
           label = paste0("tau=0.97 slope: ", round(primary_qr$slope, 3),
                          "\nReineke slope: -1.605"),
           hjust = 1, size = 3.5) +
  labs(
    x = "ln(QMD, inches)",
    y = "ln(TPA)",
    title = paste0("Self-Thinning Boundary: ", config$variant_name, " (", variant, ")"),
    caption = "Orange = tau 0.97, Blue = tau 0.95, Green = tau 0.99, Black dashed = Reineke"
  ) +
  theme_pub

ggsave(file.path(output_dir, "self_thinning_boundary.pdf"), p_boundary,
       width = 10, height = 7)

# --- Figure 2: Species SDIMAX comparison (FVS prior vs calibrated) ---
if (nrow(species_sdimax) > 0 && any(!is.na(species_sdimax$fvs_sdimax_prior))) {
  p_sdimax_compare <- species_sdimax %>%
    filter(!is.na(fvs_sdimax_prior)) %>%
    ggplot(aes(x = fvs_sdimax_prior, y = sdimax_combined)) +
    geom_point(aes(size = n_stands), alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linewidth = 0.8) +
    scale_size_continuous(name = "N Stands", range = c(1, 6)) +
    labs(
      x = "FVS Original SDIMAX",
      y = "Calibrated SDIMAX (FIA)",
      title = paste0("SDIMAX Comparison: ", config$variant_name)
    ) +
    theme_pub

  ggsave(file.path(output_dir, "sdimax_prior_vs_calibrated.pdf"), p_sdimax_compare,
         width = 8, height = 7)
}

# --- Figure 3: Species-specific self-thinning slopes ---
if (nrow(slope_summary) > 0) {
  p_slopes <- slope_summary %>%
    filter(!is.na(reineke_b_combined)) %>%
    ggplot(aes(x = reorder(factor(SPCD), reineke_b_combined),
               y = reineke_b_combined)) +
    geom_point(aes(size = n_stands), alpha = 0.7) +
    geom_hline(yintercept = 1.605, linetype = "dashed", color = "#D55E00") +
    annotate("text", x = 1, y = 1.62, label = "Reineke (1.605)",
             hjust = 0, color = "#D55E00", size = 3.5) +
    coord_flip() +
    labs(
      x = "FIA Species Code",
      y = "Self-Thinning Exponent (Reineke b)",
      title = paste0("Species Self-Thinning Slopes: ", config$variant_name)
    ) +
    theme_pub

  ggsave(file.path(output_dir, "species_selfthinning_slopes.pdf"), p_slopes,
         width = 8, height = max(4, nrow(slope_summary) * 0.3))
}

# --- Figure 4: SDMD overlay ---
sdmd_data <- tibble(qmd = seq(2, 30, by = 0.5)) %>%
  mutate(
    max_tpa_calibrated = exp(primary_qr$intercept + primary_qr$slope * log(qmd)),
    max_tpa_reineke = overall_sdimax * (qmd / 10)^(-1.605),
    mgmt_upper = 0.55 * max_tpa_calibrated,
    mgmt_lower = 0.35 * max_tpa_calibrated
  )

p_sdmd <- ggplot() +
  geom_point(data = stand_qr, aes(x = qmd, y = tpa),
             alpha = 0.1, size = 0.5, color = "grey50") +
  geom_line(data = sdmd_data, aes(x = qmd, y = max_tpa_calibrated),
            color = "#D55E00", linewidth = 1.2) +
  geom_line(data = sdmd_data, aes(x = qmd, y = max_tpa_reineke),
            color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_line(data = sdmd_data, aes(x = qmd, y = mgmt_upper),
            color = "#0072B2", linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = sdmd_data, aes(x = qmd, y = mgmt_lower),
            color = "#0072B2", linetype = "dotdash", linewidth = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Quadratic Mean Diameter (inches)",
    y = "Trees per Acre",
    title = paste0("Stand Density Management Diagram: ", config$variant_name),
    caption = "Orange = calibrated boundary, Black dashed = Reineke, Blue = management zone"
  ) +
  theme_pub

ggsave(file.path(output_dir, "sdmd_calibrated.pdf"), p_sdmd, width = 10, height = 7)

logger::log_info("Diagnostic figures saved")

# =============================================================================
# Save Full Model Summary
# =============================================================================

# Save brms model summary
summ_sdi <- posterior_summary(fit_sdi)
summary_df <- as_tibble(summ_sdi, rownames = "variable")
write_csv(summary_df, file.path(output_dir, "stand_density_summary.csv"))

# Save brms draws for posterior integration
draws_sdi <- as_draws_df(fit_sdi)
saveRDS(draws_sdi, file.path(output_dir, "stand_density_samples.rds"))

# Save overall results
overall_results <- tibble(
  variant = variant,
  variant_name = config$variant_name,
  n_stands = nrow(stand_qr),
  n_species_fitted = nrow(species_sdimax),
  overall_sdimax_qr97 = overall_sdimax,
  overall_slope_qr97 = overall_slope,
  pop_sdimax_bayes = pop_sdimax,
  pop_slope_bayes = pop_slope,
  reineke_standard = -1.605,
  sdi_param_name = sdi_param_name %||% "none",
  bamax_param_name = bamax_param_name %||% "none"
)

write_csv(overall_results, file.path(output_dir, "stand_density_overall.csv"))

# =============================================================================
# Summary Report
# =============================================================================

cat("\n")
cat("==========================================\n")
cat("Stand-Level Density Calibration Complete\n")
cat("==========================================\n")
cat("Variant:", config$variant_name, "(", variant, ")\n")
cat("\nOverall self-thinning boundary (tau=0.97):\n")
cat("  Slope:", round(overall_slope, 3), "(Reineke = -1.605)\n")
cat("  SDIMAX at QMD=10:", round(overall_sdimax), "\n")
cat("\nBayesian hierarchical model:\n")
cat("  Population SDIMAX:", round(pop_sdimax), "\n")
cat("  Population slope:", round(pop_slope, 3), "\n")
cat("  Species fitted:", nrow(species_sdimax), "\n")
cat("  SDIMAX range:", round(min(species_sdimax$sdimax_combined)),
    "to", round(max(species_sdimax$sdimax_combined)), "\n")
if (!is.null(bamax_param_name)) {
  cat("\nBAMAX calibration:", nrow(ba_qr_results), "species\n")
}
cat("\nOutput saved to:", output_dir, "\n\n")

logger::log_info("Stand-level density calibration complete for variant {variant}")
