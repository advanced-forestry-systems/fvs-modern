#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration: SDIMAX Post-Hoc Fixes
#
# Two fixes:
#   1. Revert EC and EM SDIMAX to FVS defaults (negative correlations
#      mean calibration degraded predictions relative to defaults)
#   2. Raise per-species minimum stands threshold from 30 to 100
#      for all variants' quantile regression estimates
#
# For fix #2, species with < 100 stands fall back to population-level
# (all-species pooled) SDIMAX estimate, which is more stable.
#
# Usage: Rscript calibration/R/09b_fix_sdimax.R --variant ec
#        Rscript calibration/R/09b_fix_sdimax.R --all
#

library(tidyverse)
library(jsonlite)
library(logger)

# ============================================================================
# Parse Command Line Arguments
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
variants_to_fix <- NULL
do_all <- FALSE

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variants_to_fix <- args[i + 1]
    }
    if (args[i] == "--all") {
      do_all <- TRUE
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
output_dir <- file.path(calibration_dir, "output", "variants")
config_dir <- file.path(project_root, "config")

log_file <- file.path(calibration_dir, "logs", "09b_fix_sdimax.log")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting SDIMAX post-hoc fixes")

# Variants to revert completely to defaults (negative correlations)
revert_variants <- c("ec", "em")

# Minimum stands per species for reliable QR boundary estimate
MIN_STANDS_SPECIES <- 100

# ============================================================================
# Determine which variants to process
# ============================================================================

if (do_all) {
  # Find all variants that have SDIMAX calibrations
  all_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = FALSE)
  variants_with_sdi <- all_dirs[sapply(all_dirs, function(v) {
    file.exists(file.path(output_dir, v, "species_sdimax_calibrated.csv"))
  })]
  logger::log_info("Found {length(variants_with_sdi)} variants with SDIMAX calibrations")
} else if (!is.null(variants_to_fix)) {
  variants_with_sdi <- strsplit(variants_to_fix, ",")[[1]]
} else {
  # Default: just fix EC and EM
  variants_with_sdi <- revert_variants
}

# ============================================================================
# Process each variant
# ============================================================================

results <- list()

for (variant in variants_with_sdi) {
  logger::log_info("=== Processing variant: {variant} ===")

  var_output_dir <- file.path(output_dir, variant)
  sdi_file <- file.path(var_output_dir, "species_sdimax_calibrated.csv")
  config_file <- file.path(config_dir, paste0(variant, ".json"))

  if (!file.exists(sdi_file)) {
    logger::log_warn("No SDIMAX file for {variant}, skipping")
    next
  }

  # Load current calibration
  sdi_cal <- read_csv(sdi_file, show_col_types = FALSE)
  logger::log_info("Loaded {nrow(sdi_cal)} species SDIMAX estimates")

  # Load FVS config for default values
  config <- fromJSON(config_file)
  fia_species <- unlist(config$categories$species_definitions$FIAJSP)
  spcd_to_fvs_idx <- setNames(seq_along(fia_species), fia_species)

  # Find SDI parameter
  sdi_priors <- NULL
  for (cat_name in names(config$categories)) {
    cat_data <- config$categories[[cat_name]]
    for (k in names(cat_data)) {
      if (k %in% c("SDICON", "R5SDI", "R4SDI", "FMSDI")) {
        sdi_priors <- unlist(cat_data[[k]])
      }
    }
  }

  # Back up current calibration
  backup_file <- file.path(var_output_dir, "species_sdimax_calibrated.csv.v1_bak")
  if (!file.exists(backup_file)) {
    file.copy(sdi_file, backup_file)
    logger::log_info("Backed up original to {basename(backup_file)}")
  }

  # ========================================================================
  # Fix 1: Revert to defaults if variant has negative correlation
  # ========================================================================

  if (variant %in% revert_variants) {
    logger::log_info("REVERTING {variant} to FVS defaults (negative correlation)")

    if (!is.null(sdi_priors)) {
      sdi_cal <- sdi_cal %>%
        mutate(
          fvs_sdimax_default = sapply(SPCD, function(sp) {
            idx <- spcd_to_fvs_idx[as.character(sp)]
            if (!is.na(idx) && idx <= length(sdi_priors)) sdi_priors[idx]
            else NA_real_
          }),
          # Revert: use defaults as calibrated values
          sdimax_combined_orig = sdimax_combined,
          sdimax_combined = if_else(!is.na(fvs_sdimax_default),
                                    fvs_sdimax_default,
                                    sdimax_combined),
          sdimax_source = if_else(!is.na(fvs_sdimax_default),
                                  "reverted_to_default",
                                  "calibrated_original")
        )

      n_reverted <- sum(sdi_cal$sdimax_source == "reverted_to_default", na.rm = TRUE)
      logger::log_info("Reverted {n_reverted} of {nrow(sdi_cal)} species to defaults")
    } else {
      logger::log_warn("No SDI priors found in config for {variant}")
    }

    results[[variant]] <- tibble(
      variant = variant,
      action = "reverted_to_defaults",
      n_species = nrow(sdi_cal),
      reason = "negative_correlation"
    )

  } else {
    # ========================================================================
    # Fix 2: Apply higher species threshold for QR estimates
    # ========================================================================

    if ("n_stands" %in% names(sdi_cal)) {
      n_below_threshold <- sum(sdi_cal$n_stands < MIN_STANDS_SPECIES, na.rm = TRUE)

      if (n_below_threshold > 0) {
        logger::log_info("Applying {MIN_STANDS_SPECIES}-stand threshold: {n_below_threshold} species below threshold")

        # For species below threshold, increase weight of Bayesian (pooled) estimate
        # and decrease weight of QR (noisy)
        sdi_cal <- sdi_cal %>%
          mutate(
            sdimax_combined_orig = sdimax_combined,
            # Recalculate weights with higher threshold
            w_qr_new = if_else(n_stands >= MIN_STANDS_SPECIES,
                                pmin(n_stands / 200, 0.7),  # Stricter scaling
                                pmin(n_stands / 200, 0.3)),  # Much less QR weight
            w_bayes_new = 1 - w_qr_new,
            sdimax_combined = if_else(
              !is.na(qr_sdimax_10) & !is.na(sdimax_bayes),
              w_qr_new * qr_sdimax_10 + w_bayes_new * sdimax_bayes,
              sdimax_combined  # Keep as-is if components missing
            ),
            sdimax_source = if_else(
              n_stands < MIN_STANDS_SPECIES,
              "bayes_dominant",
              "qr_dominant"
            )
          )

        logger::log_info("Adjusted {n_below_threshold} species to Bayesian-dominant estimates")
      } else {
        logger::log_info("All species have >= {MIN_STANDS_SPECIES} stands, no adjustment needed")
        sdi_cal <- sdi_cal %>%
          mutate(sdimax_source = "qr_dominant",
                 sdimax_combined_orig = sdimax_combined)
      }
    } else {
      logger::log_warn("n_stands column not found, skipping threshold fix")
      sdi_cal <- sdi_cal %>%
        mutate(sdimax_source = "original",
               sdimax_combined_orig = sdimax_combined)
    }

    results[[variant]] <- tibble(
      variant = variant,
      action = "threshold_adjusted",
      n_species = nrow(sdi_cal),
      n_below_threshold = if ("n_stands" %in% names(sdi_cal))
        sum(sdi_cal$n_stands < MIN_STANDS_SPECIES, na.rm = TRUE)
        else NA_integer_,
      reason = "raise_min_stands_to_100"
    )
  }

  # ========================================================================
  # Save updated calibration
  # ========================================================================

  write_csv(sdi_cal, sdi_file)
  logger::log_info("Saved updated SDIMAX to {sdi_file}")

  # Also re-export posteriors to JSON
  logger::log_info("Re-exporting posteriors to JSON...")
  tryCatch({
    source_script <- file.path(calibration_dir, "R", "06_posterior_to_json.R")
    if (file.exists(source_script)) {
      system2("Rscript", c(source_script, "--variant", variant),
              stdout = TRUE, stderr = TRUE)
      logger::log_info("Posterior JSON updated for {variant}")
    }
  }, error = function(e) {
    logger::log_warn("JSON export failed for {variant}: {e$message}")
  })
}

# ============================================================================
# Summary
# ============================================================================

results_df <- bind_rows(results)
write_csv(results_df, file.path(calibration_dir, "output", "sdimax_fixes_log.csv"))

cat("\n")
cat("==========================================\n")
cat("SDIMAX Post-Hoc Fixes Complete\n")
cat("==========================================\n")
cat("Variants processed:", nrow(results_df), "\n")
cat("\nActions taken:\n")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("  %s: %s (%s)\n", r$variant, r$action, r$reason))
}
cat("\nResults saved to: output/sdimax_fixes_log.csv\n\n")

logger::log_info("SDIMAX fixes complete: {nrow(results_df)} variants processed")
