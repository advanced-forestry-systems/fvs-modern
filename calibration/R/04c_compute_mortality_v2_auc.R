#!/usr/bin/env Rscript
#
# Compute mortality v2 AUC from posterior draws (draws_df objects)
#
# The v2 mortality models are saved as draws_df posterior samples, NOT brms
# fit objects. This script reconstructs predictions by:
#   1. Loading posterior draws and extracting median fixed/random effects
#   2. Loading the FULL mortality dataset
#   3. Standardizing predictors (same as 04b fitting script)
#   4. Computing logit(P(survived)) from fixed + random effects
#   5. Computing AUC via Wilcoxon-Mann-Whitney
#
# Usage: Rscript 04c_compute_mortality_v2_auc.R --variant ca
#        Rscript 04c_compute_mortality_v2_auc.R --all

library(tidyverse)

# Parse args
args <- commandArgs(trailingOnly = TRUE)
run_all <- "--all" %in% args
variant <- "ca"

if (!run_all) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) variant <- args[i + 1]
  }
}

project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (!nzchar(project_root)) {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
calibration_dir <- file.path(project_root, "calibration")
processed_data_dir <- file.path(calibration_dir, "data", "processed")

# v2 variants (those that had v2 mortality refit)
v2_variants <- c("ca", "ie", "nc", "oc", "op", "sn", "so", "tt", "ws")

if (run_all) {
  variants_to_run <- v2_variants
} else {
  variants_to_run <- variant
}

# AUC function (Wilcoxon-Mann-Whitney)
compute_auc <- function(probs, labels) {
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(probs)
  auc <- (sum(ranks[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  return(auc)
}

results <- list()

for (v in variants_to_run) {
  cat("Processing variant:", v, "\n")
  output_dir <- file.path(calibration_dir, "output", "variants", v)
  model_file <- file.path(output_dir, "mortality_v2_samples.rds")

  if (!file.exists(model_file)) {
    cat("  No v2 model found, skipping\n")
    next
  }

  # Load posterior draws
  draws <- tryCatch(readRDS(model_file), error = function(e) {
    cat("  Error loading draws:", e$message, "\n")
    NULL
  })
  if (is.null(draws)) next

  draws_tbl <- as_tibble(draws)
  cat("  Loaded draws:", ncol(draws_tbl), "parameters,", nrow(draws_tbl), "samples\n")

  # Extract median fixed effects
  get_fe <- function(nm) {
    if (nm %in% names(draws_tbl)) median(draws_tbl[[nm]], na.rm = TRUE) else 0
  }

  b0 <- get_fe("b_Intercept")
  b_dbh <- get_fe("b_DBH_std")
  b_dbh2 <- get_fe("b_IDBH_stdE2")
  b_bal <- get_fe("b_BAL_std")
  b_cr <- get_fe("b_CR_std")
  b_si <- get_fe("b_SI_std")
  b_ba <- get_fe("b_BA_std")
  # v2 may also have random slopes on DBH (check for r_SPCD_grouped columns)
  b_dbh_bal <- get_fe("b_DBH_std:BAL_std")

  cat("  Fixed effects: Intercept=", round(b0, 3),
      " DBH=", round(b_dbh, 3), " DBH2=", round(b_dbh2, 3),
      " BAL=", round(b_bal, 3), " CR=", round(b_cr, 3),
      " SI=", round(b_si, 3), " BA=", round(b_ba, 3), "\n")

  # Extract species random effects (intercepts and potentially slopes)
  re_int_cols <- names(draws_tbl)[grepl("^r_SPCD_grouped\\[", names(draws_tbl)) &
                                   grepl("Intercept", names(draws_tbl))]
  re_slope_cols <- names(draws_tbl)[grepl("^r_SPCD_grouped\\[", names(draws_tbl)) &
                                     grepl("DBH_std", names(draws_tbl)) &
                                     !grepl("Intercept", names(draws_tbl))]

  # Also check for r_SPCD patterns (without _grouped)
  if (length(re_int_cols) == 0) {
    re_int_cols <- names(draws_tbl)[grepl("^r_SPCD\\[", names(draws_tbl)) &
                                     grepl("Intercept", names(draws_tbl))]
    re_slope_cols <- names(draws_tbl)[grepl("^r_SPCD\\[", names(draws_tbl)) &
                                       grepl("DBH_std", names(draws_tbl)) &
                                       !grepl("Intercept", names(draws_tbl))]
  }

  sp_re_int <- setNames(
    sapply(re_int_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD(_grouped)?\\[(.+),Intercept\\]", "\\2", re_int_cols)
  )

  sp_re_slope <- if (length(re_slope_cols) > 0) {
    setNames(
      sapply(re_slope_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
      gsub("r_SPCD(_grouped)?\\[(.+),DBH_std\\]", "\\2", re_slope_cols)
    )
  } else {
    NULL
  }

  cat("  Random effects:", length(sp_re_int), "species intercepts")
  if (!is.null(sp_re_slope)) cat(",", length(sp_re_slope), "species slopes")
  cat("\n")

  # Load FULL dataset
  data_file <- file.path(processed_data_dir, v, "mortality.csv")
  if (!file.exists(data_file)) {
    data_file <- file.path(processed_data_dir, v, "diameter_growth.csv")
  }
  if (!file.exists(data_file)) {
    cat("  No data file found, skipping\n")
    next
  }

  mort_data <- read_csv(data_file, show_col_types = FALSE)

  # Reconstruct survived column
  if ("survived" %in% names(mort_data)) {
    # good
  } else if ("died" %in% names(mort_data)) {
    mort_data <- mort_data %>% mutate(survived = !as.logical(died))
  } else if ("STATUSCD_t2" %in% names(mort_data)) {
    mort_data <- mort_data %>% mutate(survived = STATUSCD_t2 == 1)
  } else if ("DIA_t2" %in% names(mort_data)) {
    mort_data <- mort_data %>% mutate(survived = !is.na(DIA_t2) & DIA_t2 > 0)
  }

  # Standardize columns
  if ("DIA" %in% names(mort_data) & !("DIA_t1" %in% names(mort_data))) {
    mort_data <- mort_data %>% rename(DIA_t1 = DIA)
  }
  if (!("CR_pct" %in% names(mort_data)) & "CR" %in% names(mort_data)) {
    mort_data <- mort_data %>% mutate(CR_pct = CR * 100)
  }
  if (!("CR_pct" %in% names(mort_data)) & "CR_t1" %in% names(mort_data)) {
    mort_data <- mort_data %>% mutate(CR_pct = CR_t1)
  }

  # Filter complete cases
  required_raw <- c("DIA_t1", "BAL", "SI", "BA")
  if ("CR_pct" %in% names(mort_data)) required_raw <- c(required_raw, "CR_pct")
  mort_data <- mort_data %>%
    filter(!is.na(survived)) %>%
    filter(if_all(all_of(required_raw[required_raw %in% names(mort_data)]), ~ !is.na(.x)))

  if (nrow(mort_data) < 100) {
    cat("  Too few complete observations:", nrow(mort_data), "\n")
    next
  }

  # Create SPCD_grouped: lump rare species (same as fitting)
  sp_counts <- mort_data %>% count(SPCD) %>% arrange(desc(n))
  rare_threshold <- 20
  rare_sp <- sp_counts %>% filter(n < rare_threshold) %>% pull(SPCD)
  if (length(rare_sp) > 0) {
    mort_data$SPCD_grouped <- ifelse(mort_data$SPCD %in% rare_sp,
                                      "OTHER", as.character(mort_data$SPCD))
  } else {
    mort_data$SPCD_grouped <- as.character(mort_data$SPCD)
  }

  # Standardize predictors (same as 04b)
  mort_data <- mort_data %>%
    mutate(
      DBH_std = as.numeric(scale(DIA_t1)),
      BAL_std = as.numeric(scale(BAL)),
      CR_std  = if ("CR_pct" %in% names(.)) as.numeric(scale(CR_pct / 100)) else 0,
      SI_std  = as.numeric(scale(log(pmax(SI, 1)))),
      BA_std  = as.numeric(scale(BA))
    )

  cat("  Full dataset:", nrow(mort_data), "obs,", sum(!mort_data$survived), "dead\n")

  # Compute predictions using posterior medians
  mort_data <- mort_data %>%
    mutate(
      sp_int   = ifelse(SPCD_grouped %in% names(sp_re_int),
                         sp_re_int[SPCD_grouped], 0),
      sp_slope = if (!is.null(sp_re_slope)) {
                   ifelse(SPCD_grouped %in% names(sp_re_slope),
                          sp_re_slope[SPCD_grouped], 0)
                 } else { 0 },
      logit_surv = b0 + (b_dbh + sp_slope) * DBH_std + b_dbh2 * DBH_std^2 +
                   b_bal * BAL_std + b_cr * CR_std + b_si * SI_std +
                   b_ba * BA_std + b_dbh_bal * DBH_std * BAL_std + sp_int,
      pred_surv = plogis(logit_surv)
    )

  # AUC: P(dead) vs actual dead
  mort_probs <- 1 - mort_data$pred_surv
  mort_labels <- !mort_data$survived

  auc_val <- compute_auc(mort_probs, mort_labels)
  cat("  AUC (v2, full data):", round(auc_val, 4), "\n")

  results[[v]] <- tibble(
    variant = toupper(v),
    auc_v2 = auc_val,
    n_obs = nrow(mort_data),
    n_dead = sum(mort_labels, na.rm = TRUE),
    mort_rate = mean(mort_labels, na.rm = TRUE),
    n_species_re = length(sp_re_int),
    has_random_slopes = !is.null(sp_re_slope),
    note = "full_data_manual_prediction"
  )
}

# Combine and save
results_df <- bind_rows(results)
out_file <- file.path(calibration_dir, "output", "comparisons", "mortality_v2_auc_full.csv")
write_csv(results_df, out_file)
cat("\nResults saved to:", out_file, "\n")
print(results_df)
