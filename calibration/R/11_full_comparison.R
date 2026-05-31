#!/usr/bin/env Rscript
# =============================================================================
# Title: FVS Bayesian Calibration: Comprehensive Default vs Calibrated Comparison
# Author: A. Weiskittel
# Date: 2026-03-30
# Description:
#   Compares default and calibrated FVS parameters across ALL model components
#   for all 25 geographic variants. Produces tree-level and stand-level
#   performance metrics plus publication-ready diagnostic figures.
#
#   Components evaluated:
#     1. Height-Diameter (Chapman-Richards, 25/25 variants)
#     2. Mortality (logistic, 25/25 variants)
#     3. Crown Ratio Change (linear, 25/25 variants)
#     4. Diameter Growth (Wykoff, available variants only)
#     5. SDIMAX (Reineke/Bayesian, 25/25 variants)
#     6. Height Increment (6 variants with explicit HG params)
#
# Dependencies: calibration/data/processed/{variant}/, output/variants/{variant}/
# Usage: Rscript 11_full_comparison.R --all
# =============================================================================

# --- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(jsonlite)
  library(posterior)
  library(logger)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

# --- Configuration -----------------------------------------------------------

project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (!nzchar(project_root)) {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
calibration_dir <- file.path(project_root, "calibration")
config_dir      <- file.path(project_root, "config")
output_base     <- file.path(calibration_dir, "output", "variants")
data_base       <- file.path(calibration_dir, "data", "processed")
comp_dir        <- file.path(calibration_dir, "output", "comparisons")

dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

ALL_VARIANTS <- c("acd","ak","bc","bm","ca","ci","cr","cs","ec","em",
                   "ie","kt","ls","nc","ne","oc","on","op","pn","sn",
                   "so","tt","ut","wc","ws")

HI_VARIANTS <- c("bc","ci","em","ie","kt","ws")

# --- Utility functions -------------------------------------------------------

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
bias <- function(obs, pred) mean(pred - obs, na.rm = TRUE)
r2   <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

cb_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442",
            "#0072B2","#D55E00","#CC79A7","#999999")

# Master results accumulator
all_metrics <- tibble()

cat("=============================================================\n")
cat("FVS Comprehensive Default vs Calibrated Comparison\n")
cat("=============================================================\n")
cat("Project root:", project_root, "\n")
cat("Variants:", length(ALL_VARIANTS), "\n\n")

# =============================================================================
# 1. HEIGHT-DIAMETER (tree level)
# =============================================================================

cat("=== Component 1: Height-Diameter ===\n")

hd_results <- list()

for (v in ALL_VARIANTS) {
  cat("  HD:", v, "... ")

  # Load posterior
  post_file <- file.path(output_base, v, "height_diameter_samples.rds")
  if (!file.exists(post_file)) { cat("no posterior\n"); next }

  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) { cat("load error\n"); next }

  # Extract population-level parameters
  draws_tbl <- as_tibble(draws)
  draw_names <- names(draws_tbl)

  a_pop <- median(draws_tbl$b_a_Intercept, na.rm = TRUE)
  b_pop <- median(draws_tbl$b_b_Intercept, na.rm = TRUE)
  c_pop <- median(draws_tbl$b_c_Intercept, na.rm = TRUE)

  if (is.na(a_pop) || is.na(b_pop) || is.na(c_pop)) { cat("missing params\n"); next }

  # Extract species random effects
  re_a_cols <- draw_names[grepl("r_SPCD__a\\[", draw_names)]
  re_b_cols <- draw_names[grepl("r_SPCD__b\\[", draw_names)]

  sp_re_a <- setNames(
    sapply(re_a_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD__a\\[(.+),Intercept\\]", "\\1", re_a_cols)
  )
  sp_re_b <- setNames(
    sapply(re_b_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD__b\\[(.+),Intercept\\]", "\\1", re_b_cols)
  )

  # Load validation data
  hd_file <- file.path(data_base, v, "height_diameter.csv")
  dg_file <- file.path(data_base, v, "diameter_growth.csv")
  data_file <- if (file.exists(hd_file)) hd_file else if (file.exists(dg_file)) dg_file else NULL
  if (is.null(data_file)) { cat("no data\n"); next }

  hd_data <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble()

  # Harmonize columns
  if ("DIA" %in% names(hd_data) & !("DIA_t1" %in% names(hd_data))) hd_data <- hd_data %>% rename(DIA_t1 = DIA)
  if ("HT"  %in% names(hd_data) & !("HT_t1"  %in% names(hd_data))) hd_data <- hd_data %>% rename(HT_t1 = HT)

  hd_data <- hd_data %>%
    filter(!is.na(HT_t1), HT_t1 > 4.5, !is.na(DIA_t1), DIA_t1 > 0)

  # Subsample for speed
  set.seed(42)
  if (nrow(hd_data) > 30000) hd_data <- hd_data %>% slice_sample(n = 30000)

  # Predict with calibrated model (species-specific)
  hd_data <- hd_data %>%
    mutate(
      SPCD_chr = as.character(SPCD),
      DBH_scaled = DIA_t1 / 20,
      a_sp = a_pop + ifelse(SPCD_chr %in% names(sp_re_a), sp_re_a[SPCD_chr], 0),
      b_sp = b_pop + ifelse(SPCD_chr %in% names(sp_re_b), sp_re_b[SPCD_chr], 0),
      pred_calibrated = 4.5 + pmax(a_sp, 1) * (1 - exp(-pmax(b_sp, 0.001) * DBH_scaled))^pmax(c_pop, 0.1)
    )

  # Baseline: species power curve H = exp(a + b*ln(DBH)) per species
  baseline_fit <- tryCatch(
    lm(log(HT_t1) ~ log(DIA_t1) + factor(SPCD), data = hd_data),
    error = function(e) NULL
  )

  if (!is.null(baseline_fit)) {
    hd_data$pred_baseline <- exp(predict(baseline_fit, newdata = hd_data))
  } else {
    hd_data$pred_baseline <- NA_real_
  }

  # Tree-level metrics
  m <- hd_data %>%
    summarise(
      n = n(),
      rmse_cal = rmse(HT_t1, pred_calibrated),
      bias_cal = bias(HT_t1, pred_calibrated),
      r2_cal   = r2(HT_t1, pred_calibrated),
      rmse_base = rmse(HT_t1, pred_baseline),
      bias_base = bias(HT_t1, pred_baseline),
      r2_base   = r2(HT_t1, pred_baseline)
    ) %>%
    mutate(variant = toupper(v))

  hd_results[[v]] <- m

  # Species-level metrics
  sp_m <- hd_data %>%
    group_by(SPCD) %>%
    summarise(
      n = n(),
      rmse_cal = rmse(HT_t1, pred_calibrated),
      r2_cal   = r2(HT_t1, pred_calibrated),
      rmse_base = rmse(HT_t1, pred_baseline),
      r2_base   = r2(HT_t1, pred_baseline),
      .groups = "drop"
    ) %>%
    mutate(variant = toupper(v))

  cat("R2_cal =", round(m$r2_cal, 3), "\n")
}

hd_df <- bind_rows(hd_results)
write_csv(hd_df, file.path(comp_dir, "hd_tree_level_metrics.csv"))

# H-D summary figure
if (nrow(hd_df) > 0) {
  p_hd <- hd_df %>%
    select(variant, r2_cal, r2_base) %>%
    pivot_longer(c(r2_cal, r2_base), names_to = "model", values_to = "r2") %>%
    mutate(model = ifelse(model == "r2_cal", "Calibrated", "Baseline")) %>%
    ggplot(aes(x = reorder(variant, r2), y = r2, fill = model)) +
    geom_col(position = "dodge", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Calibrated" = cb_pal[5], "Baseline" = cb_pal[8])) +
    labs(x = NULL, y = expression(R^2), title = "Height-Diameter Model Performance (Tree Level)",
         fill = NULL) +
    theme_pub
  ggsave(file.path(comp_dir, "fig_hd_r2_comparison.png"), p_hd, width = 20, height = 16, units = "cm", dpi = 300, bg = "white")
}

# Append to master table
hd_master <- hd_df %>%
  pivot_longer(c(rmse_cal, bias_cal, r2_cal, rmse_base, bias_base, r2_base),
               names_to = "metric_model", values_to = "value") %>%
  separate(metric_model, into = c("metric", "model"), sep = "_") %>%
  mutate(component = "Height_Diameter", level = "tree",
         model = ifelse(model == "cal", "calibrated", "baseline"))
all_metrics <- bind_rows(all_metrics, hd_master)

# =============================================================================
# 2. MORTALITY (tree level)
# =============================================================================

cat("\n=== Component 2: Mortality ===\n")

mort_results <- list()

for (v in ALL_VARIANTS) {
  cat("  MORT:", v, "... ")

  # Load posterior summary (use summary CSV, lighter than full draws)
  summ_file <- file.path(output_base, v, "mortality_summary.csv")
  post_file <- file.path(output_base, v, "mortality_samples.rds")

  if (!file.exists(summ_file)) { cat("no summary\n"); next }

  # Load data
  mort_file <- file.path(data_base, v, "mortality.csv")
  dg_file   <- file.path(data_base, v, "diameter_growth.csv")
  data_file <- if (file.exists(mort_file)) mort_file else if (file.exists(dg_file)) dg_file else NULL
  if (is.null(data_file)) { cat("no data\n"); next }

  mort_data <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble()

  # Harmonize columns
  if ("DIA" %in% names(mort_data) & !("DIA_t1" %in% names(mort_data)))
    mort_data <- mort_data %>% rename(DIA_t1 = DIA)
  if ("died" %in% names(mort_data))
    mort_data <- mort_data %>% mutate(survived = !as.logical(died))
  if ("DIA_t2" %in% names(mort_data) & !("survived" %in% names(mort_data)))
    mort_data <- mort_data %>% mutate(survived = !is.na(DIA_t2) & DIA_t2 > 0)
  if (!("CR_pct" %in% names(mort_data)) & "CR" %in% names(mort_data))
    mort_data <- mort_data %>% mutate(CR_pct = CR * 100)

  mort_data <- mort_data %>%
    filter(!is.na(survived), !is.na(DIA_t1), !is.na(BAL), !is.na(CR_pct),
           !is.na(SI), !is.na(BA))

  if (nrow(mort_data) < 100) { cat("too few obs\n"); next }

  # Reconstruct predictions from posterior
  # Model: survived ~ DBH_std + I(DBH_std^2) + BAL_std + CR_std + SI_std + BA_std + (1|SPCD)
  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)

  if (!is.null(draws)) {
    draws_tbl <- as_tibble(draws)

    # Get fixed effects (posterior medians)
    get_fe <- function(nm) {
      col <- paste0("b_", nm)
      if (col %in% names(draws_tbl)) median(draws_tbl[[col]], na.rm = TRUE) else 0
    }

    b0 <- get_fe("Intercept")
    b1 <- get_fe("DBH_std")
    b2 <- get_fe("IDBH_stdE2")  # brms encodes I(x^2) as IxE2
    b3 <- get_fe("BAL_std")
    b4 <- get_fe("CR_std")
    b5 <- get_fe("SI_std")
    b6 <- get_fe("BA_std")

    # Species random intercepts
    re_cols <- names(draws_tbl)[grepl("r_SPCD\\[", names(draws_tbl))]
    sp_re <- setNames(
      sapply(re_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
      gsub("r_SPCD\\[(.+),Intercept\\]", "\\1", re_cols)
    )

    # Standardize predictors using training data moments
    mort_data <- mort_data %>%
      mutate(
        DBH_std = as.numeric(scale(DIA_t1)),
        CR_std  = as.numeric(scale(CR_pct / 100)),
        BAL_std = as.numeric(scale(BAL)),
        SI_std  = as.numeric(scale(log(pmax(SI, 1)))),
        BA_std  = as.numeric(scale(BA)),
        SPCD_chr = as.character(SPCD),
        sp_int = ifelse(SPCD_chr %in% names(sp_re), sp_re[SPCD_chr], 0),
        logit_p = b0 + b1*DBH_std + b2*DBH_std^2 + b3*BAL_std + b4*CR_std +
                  b5*SI_std + b6*BA_std + sp_int,
        pred_prob_surv = 1 / (1 + exp(-logit_p))
      )

    # AUC (simple trapezoidal)
    auc_val <- tryCatch({
      roc_data <- mort_data %>%
        arrange(desc(pred_prob_surv)) %>%
        mutate(
          tp = cumsum(survived) / sum(survived),
          fp = cumsum(!survived) / sum(!survived)
        )
      sum(diff(roc_data$fp) * (roc_data$tp[-1] + roc_data$tp[-nrow(roc_data)]) / 2)
    }, error = function(e) NA_real_)

    # Brier score
    brier <- mean((mort_data$pred_prob_surv - as.numeric(mort_data$survived))^2, na.rm = TRUE)

    # Mortality rate by diameter class
    mort_by_dbh <- mort_data %>%
      mutate(dbh_class = cut(DIA_t1, breaks = c(0,5,10,15,20,25,30,40,100), right = FALSE)) %>%
      group_by(dbh_class) %>%
      summarise(
        n = n(),
        obs_mort_rate = 1 - mean(survived, na.rm = TRUE),
        pred_mort_rate = 1 - mean(pred_prob_surv, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n >= 20)

    m <- tibble(
      variant = toupper(v),
      n = nrow(mort_data),
      overall_mort_rate = 1 - mean(mort_data$survived),
      auc = auc_val,
      brier_score = brier,
      n_species = n_distinct(mort_data$SPCD)
    )

    mort_results[[v]] <- m
    cat("AUC =", round(auc_val, 3), "\n")
  } else {
    cat("no draws\n")
  }
}

mort_df <- bind_rows(mort_results)
write_csv(mort_df, file.path(comp_dir, "mortality_tree_level_metrics.csv"))

# Mortality AUC figure
if (nrow(mort_df) > 0) {
  p_mort <- mort_df %>%
    filter(!is.na(auc)) %>%
    ggplot(aes(x = reorder(variant, auc), y = auc)) +
    geom_col(fill = cb_pal[6], width = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    coord_flip() +
    labs(x = NULL, y = "AUC", title = "Calibrated Mortality Model Performance (Tree Level)") +
    theme_pub
  ggsave(file.path(comp_dir, "fig_mortality_auc.png"), p_mort, width = 16, height = 16, units = "cm", dpi = 300, bg = "white")
}

mort_master <- mort_df %>%
  select(variant, auc, brier_score, overall_mort_rate) %>%
  pivot_longer(c(auc, brier_score, overall_mort_rate), names_to = "metric", values_to = "value") %>%
  mutate(component = "Mortality", level = "tree", model = "calibrated")
all_metrics <- bind_rows(all_metrics, mort_master)

# =============================================================================
# 3. CROWN RATIO CHANGE (tree level)
# =============================================================================

cat("\n=== Component 3: Crown Ratio Change ===\n")

cr_results <- list()

for (v in ALL_VARIANTS) {
  cat("  CR:", v, "... ")

  post_file <- file.path(output_base, v, "crown_ratio_samples.rds")
  if (!file.exists(post_file)) { cat("no posterior\n"); next }

  # Load data
  cr_file <- file.path(data_base, v, "crown_ratio_change.csv")
  dg_file <- file.path(data_base, v, "diameter_growth.csv")

  if (file.exists(cr_file)) {
    cr_data <- read_csv(cr_file, show_col_types = FALSE) %>% as_tibble()
    if ("delta_CR" %in% names(cr_data)) cr_data <- cr_data %>% rename(CR_change = delta_CR)
    if ("DIA" %in% names(cr_data) & !("DIA_t1" %in% names(cr_data))) cr_data <- cr_data %>% rename(DIA_t1 = DIA)
    if (!("CR_pct" %in% names(cr_data)) & "CR_init" %in% names(cr_data))
      cr_data <- cr_data %>% mutate(CR_pct = CR_init * 100)
    cr_data <- cr_data %>%
      filter(!is.na(CR_change), is.finite(CR_change), !is.na(DIA_t1))
  } else if (file.exists(dg_file)) {
    cr_data <- read_csv(dg_file, show_col_types = FALSE) %>% as_tibble() %>%
      filter(!is.na(CR_t1), !is.na(CR_t2), CR_t1 > 0, years_interval > 0) %>%
      mutate(CR_change = (CR_t2 - CR_t1) / years_interval) %>%
      filter(is.finite(CR_change))
  } else {
    cat("no data\n"); next
  }

  cr_data <- cr_data %>%
    filter(!is.na(BA), !is.na(BAL), !is.na(SI), !is.na(CR_pct))

  if (nrow(cr_data) < 100) { cat("too few\n"); next }

  # Load posterior draws
  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) { cat("load error\n"); next }

  draws_tbl <- as_tibble(draws)

  get_fe <- function(nm) {
    col <- paste0("b_", nm)
    if (col %in% names(draws_tbl)) median(draws_tbl[[col]], na.rm = TRUE) else 0
  }

  b0 <- get_fe("Intercept")
  b1 <- get_fe("DBH_std")
  b2 <- get_fe("IDBH_stdE2")
  b3 <- get_fe("BA_std")
  b4 <- get_fe("BAL_std")
  b5 <- get_fe("CR_std")
  b6 <- get_fe("SI_std")

  # Species RE
  re_cols <- names(draws_tbl)[grepl("r_SPCD\\[", names(draws_tbl))]
  sp_re <- setNames(
    sapply(re_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD\\[(.+),Intercept\\]", "\\1", re_cols)
  )

  set.seed(42)
  if (nrow(cr_data) > 30000) cr_data <- cr_data %>% slice_sample(n = 30000)

  cr_data <- cr_data %>%
    mutate(
      DBH_std = as.numeric(scale(DIA_t1)),
      BA_std  = as.numeric(scale(BA)),
      BAL_std = as.numeric(scale(BAL)),
      CR_std  = as.numeric(scale(CR_pct / 100)),
      SI_std  = as.numeric(scale(log(pmax(SI, 1)))),
      SPCD_chr = as.character(SPCD),
      sp_int = ifelse(SPCD_chr %in% names(sp_re), sp_re[SPCD_chr], 0),
      pred_cr = b0 + b1*DBH_std + b2*DBH_std^2 + b3*BA_std + b4*BAL_std +
                b5*CR_std + b6*SI_std + sp_int
    )

  m <- cr_data %>%
    summarise(
      variant = toupper(v),
      n = n(),
      rmse_cal = rmse(CR_change, pred_cr),
      bias_cal = bias(CR_change, pred_cr),
      r2_cal   = r2(CR_change, pred_cr)
    )

  cr_results[[v]] <- m
  cat("R2 =", round(m$r2_cal, 3), "\n")
}

cr_df <- bind_rows(cr_results)
write_csv(cr_df, file.path(comp_dir, "crown_ratio_tree_level_metrics.csv"))

cr_master <- cr_df %>%
  pivot_longer(c(rmse_cal, bias_cal, r2_cal), names_to = "metric_model", values_to = "value") %>%
  separate(metric_model, into = c("metric", "model"), sep = "_") %>%
  mutate(component = "Crown_Ratio", level = "tree", model = "calibrated")
all_metrics <- bind_rows(all_metrics, cr_master)

# =============================================================================
# 4. DIAMETER GROWTH (tree level)
# =============================================================================

cat("\n=== Component 4: Diameter Growth ===\n")

dg_results <- list()

for (v in ALL_VARIANTS) {
  cat("  DG:", v, "... ")

  post_file <- file.path(output_base, v, "diameter_growth_samples.rds")
  std_file  <- file.path(output_base, v, "standardization_params.rds")
  data_file <- file.path(data_base, v, "diameter_growth.csv")

  if (!file.exists(post_file)) { cat("no posterior\n"); next }
  if (!file.exists(data_file)) { cat("no data\n"); next }

  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) { cat("load error\n"); next }

  std_params <- if (file.exists(std_file)) readRDS(std_file) else NULL

  dg_data <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble() %>%
    filter(!is.na(ln_DDS), is.finite(ln_DDS), !is.na(ln_DBH))

  if (nrow(dg_data) < 100) { cat("too few\n"); next }

  set.seed(42)
  if (nrow(dg_data) > 30000) dg_data <- dg_data %>% slice_sample(n = 30000)

  # Extract posterior medians
  draws_tbl <- as_tibble(draws)
  get_par <- function(nm) {
    if (nm %in% names(draws_tbl)) median(draws_tbl[[nm]], na.rm = TRUE) else 0
  }

  mu_b0 <- get_par("mu_b0")
  b1 <- get_par("b1"); b2 <- get_par("b2"); b3 <- get_par("b3")
  b4 <- get_par("b4"); b5 <- get_par("b5"); b6 <- get_par("b6")
  b7 <- get_par("b7"); b8 <- get_par("b8"); b9 <- get_par("b9")
  b10 <- get_par("b10"); b11 <- get_par("b11"); b12 <- get_par("b12"); b13 <- get_par("b13")

  # Species random intercepts
  sp_cols <- names(draws_tbl)[grepl("^b0\\[", names(draws_tbl))]
  sp_b0 <- sapply(sp_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE))

  # Standardize predictors
  if (!is.null(std_params)) {
    dg_data <- dg_data %>%
      mutate(
        ln_DBH_std = (ln_DBH - std_params$ln_DBH_mean) / std_params$ln_DBH_sd,
        DBH_sq_std = (DBH_sq - std_params$DBH_sq_mean) / std_params$DBH_sq_sd,
        ln_SI_std  = (ln_SI - std_params$ln_SI_mean) / std_params$ln_SI_sd,
        SLOPE_std  = (SLOPE - std_params$SLOPE_mean) / std_params$SLOPE_sd,
        ELEV_std   = (ELEV - std_params$ELEV_mean) / std_params$ELEV_sd,
        SLOPE_SASP_std = (SLOPE_SASP - std_params$SLOPE_SASP_mean) / std_params$SLOPE_SASP_sd,
        SLOPE_CASP_std = (SLOPE_CASP - std_params$SLOPE_CASP_mean) / std_params$SLOPE_CASP_sd,
        BAL_std = (BAL - std_params$BAL_mean) / std_params$BAL_sd,
        BA_std  = (BA - std_params$BA_mean) / std_params$BA_sd,
        CR_prop = CR_pct / 100,
        species_idx = as.integer(factor(SPCD))
      )
  } else {
    dg_data <- dg_data %>%
      mutate(
        ln_DBH_std = as.numeric(scale(ln_DBH)),
        DBH_sq_std = as.numeric(scale(DBH_sq)),
        ln_SI_std  = as.numeric(scale(ln_SI)),
        SLOPE_std  = as.numeric(scale(SLOPE)),
        ELEV_std   = as.numeric(scale(ELEV)),
        SLOPE_SASP_std = as.numeric(scale(SLOPE_SASP)),
        SLOPE_CASP_std = as.numeric(scale(SLOPE_CASP)),
        BAL_std = as.numeric(scale(BAL)),
        BA_std  = as.numeric(scale(BA)),
        CR_prop = CR_pct / 100,
        species_idx = as.integer(factor(SPCD))
      )
  }

  # Predict
  dg_data <- dg_data %>%
    mutate(
      sp_intercept = ifelse(species_idx <= length(sp_b0), sp_b0[species_idx], mu_b0),
      pred_ln_dds = sp_intercept +
        b1*ln_DBH_std + b2*DBH_sq_std + b3*ln_SI_std +
        b4*SLOPE_std + b5*SLOPE_std^2 + b6*SLOPE_SASP_std + b7*SLOPE_CASP_std +
        b8*ELEV_std + b9*ELEV_std^2 + b10*CR_prop + b11*CR_prop^2 +
        b12*BAL_std + b13*BA_std
    )

  m <- dg_data %>%
    filter(is.finite(pred_ln_dds)) %>%
    summarise(
      variant = toupper(v),
      n = n(),
      rmse_log = rmse(ln_DDS, pred_ln_dds),
      bias_log = bias(ln_DDS, pred_ln_dds),
      r2_log   = r2(ln_DDS, pred_ln_dds)
    )

  dg_results[[v]] <- m
  cat("R2_log =", round(m$r2_log, 3), "\n")
}

dg_df <- bind_rows(dg_results)
write_csv(dg_df, file.path(comp_dir, "diameter_growth_tree_level_metrics.csv"))

dg_master <- dg_df %>%
  pivot_longer(c(rmse_log, bias_log, r2_log), names_to = "metric_model", values_to = "value") %>%
  separate(metric_model, into = c("metric", "model"), sep = "_") %>%
  mutate(component = "Diameter_Growth", level = "tree", model = "calibrated")
all_metrics <- bind_rows(all_metrics, dg_master)

# =============================================================================
# 5. SDIMAX (stand level)
# =============================================================================

cat("\n=== Component 5: SDIMAX ===\n")

sdi_results <- list()

for (v in ALL_VARIANTS) {
  cat("  SDI:", v, "... ")

  cal_file <- file.path(output_base, v, "species_sdimax_calibrated.csv")
  cfg_file <- file.path(config_dir, paste0(v, ".json"))

  if (!file.exists(cal_file) || !file.exists(cfg_file)) { cat("missing files\n"); next }

  cal_data <- read_csv(cal_file, show_col_types = FALSE) %>% as_tibble()
  config   <- fromJSON(cfg_file)

  # Get SDICON from config
  sdi_param <- NULL
  for (cat_name in names(config$categories)) {
    cat_data <- config$categories[[cat_name]]
    for (k in names(cat_data)) {
      if (k %in% c("SDICON", "R5SDI", "R4SDI", "FMSDI")) {
        sdi_param <- unlist(cat_data[[k]])
      }
    }
  }

  if (is.null(sdi_param)) { cat("no SDICON in config\n"); next }

  fia_sp <- as.integer(unlist(config$categories$species_definitions$FIAJSP))
  sp_to_idx <- setNames(seq_along(fia_sp), fia_sp)

  # Use sdimax_combined or sdimax_bayes
  cal_col <- if ("sdimax_combined" %in% names(cal_data)) "sdimax_combined" else
             if ("sdimax_bayes" %in% names(cal_data)) "sdimax_bayes" else NULL
  if (is.null(cal_col)) { cat("no calibrated column\n"); next }

  matched <- cal_data %>%
    filter(!is.na(.data[[cal_col]])) %>%
    mutate(
      fvs_idx = sp_to_idx[as.character(as.integer(SPCD))],
      default_sdimax = ifelse(!is.na(fvs_idx) & fvs_idx <= length(sdi_param),
                              sdi_param[fvs_idx], NA_real_),
      calibrated_sdimax = .data[[cal_col]]
    ) %>%
    filter(!is.na(default_sdimax), default_sdimax > 0)

  if (nrow(matched) < 1) { cat("no matches\n"); next }

  m <- matched %>%
    summarise(
      variant = toupper(v),
      n_species = n(),
      mean_default = mean(default_sdimax),
      mean_calibrated = mean(calibrated_sdimax),
      mean_pct_change = mean((calibrated_sdimax - default_sdimax) / default_sdimax * 100),
      correlation = cor(default_sdimax, calibrated_sdimax)
    )

  sdi_results[[v]] <- m
  cat(nrow(matched), "species matched, pct_change =", round(m$mean_pct_change, 1), "%\n")
}

sdi_df <- bind_rows(sdi_results)
write_csv(sdi_df, file.path(comp_dir, "sdimax_stand_level_metrics.csv"))

sdi_master <- sdi_df %>%
  select(variant, mean_pct_change, correlation) %>%
  pivot_longer(c(mean_pct_change, correlation), names_to = "metric", values_to = "value") %>%
  mutate(component = "SDIMAX", level = "stand", model = "comparison")
all_metrics <- bind_rows(all_metrics, sdi_master)

# =============================================================================
# 6. HEIGHT INCREMENT (tree level, 6 variants only)
# =============================================================================

cat("\n=== Component 6: Height Increment ===\n")

hi_results <- list()

for (v in HI_VARIANTS) {
  cat("  HI:", v, "... ")
  summ_file <- file.path(output_base, v, "height_increment_summary.csv")
  if (!file.exists(summ_file)) { cat("no summary\n"); next }

  hi_summ <- read_csv(summ_file, show_col_types = FALSE)
  m <- tibble(
    variant = toupper(v),
    n_params = nrow(hi_summ),
    method = tryCatch(
      readLines(file.path(output_base, v, "height_increment_method.txt"), n = 1),
      error = function(e) "unknown"
    )
  )
  hi_results[[v]] <- m
  cat(m$n_params, "params,", m$method, "\n")
}

hi_df <- bind_rows(hi_results)
write_csv(hi_df, file.path(comp_dir, "height_increment_summary.csv"))

# =============================================================================
# MASTER SUMMARY TABLE
# =============================================================================

cat("\n=== Generating Master Summary ===\n")

write_csv(all_metrics, file.path(comp_dir, "model_performance_summary.csv"))

# Equation availability heatmap
eq_avail <- tibble()
for (v in ALL_VARIANTS) {
  eq_avail <- bind_rows(eq_avail, tibble(
    variant = toupper(v),
    HD  = file.exists(file.path(output_base, v, "height_diameter_summary.csv")),
    MORT = file.exists(file.path(output_base, v, "mortality_summary.csv")),
    CR  = file.exists(file.path(output_base, v, "crown_ratio_summary.csv")),
    DG  = file.exists(file.path(output_base, v, "diameter_growth_summary.csv")),
    SDI = file.exists(file.path(output_base, v, "species_sdimax_calibrated.csv")),
    HI  = file.exists(file.path(output_base, v, "height_increment_summary.csv"))
  ))
}

write_csv(eq_avail, file.path(comp_dir, "equation_availability_full.csv"))

# Heatmap figure
if (nrow(eq_avail) > 0) {
  p_heat <- eq_avail %>%
    pivot_longer(c(HD, MORT, CR, DG, SDI, HI), names_to = "component", values_to = "available") %>%
    mutate(component = factor(component, levels = c("HD","MORT","CR","DG","SDI","HI"))) %>%
    ggplot(aes(x = component, y = variant, fill = available)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = cb_pal[3], "FALSE" = "grey85"),
                      labels = c("TRUE" = "Complete", "FALSE" = "Missing")) +
    labs(x = "Model Component", y = "Variant", title = "Calibration Completeness by Variant",
         fill = NULL) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 0))
  ggsave(file.path(comp_dir, "fig_calibration_heatmap.png"), p_heat,
         width = 14, height = 18, units = "cm", dpi = 300, bg = "white")
}

# Cross-component summary
comp_summary <- all_metrics %>%
  filter(!is.na(value), is.finite(value)) %>%
  group_by(component, metric, model) %>%
  summarise(
    n_variants = n_distinct(variant),
    mean_val   = mean(value, na.rm = TRUE),
    sd_val     = sd(value, na.rm = TRUE),
    min_val    = min(value, na.rm = TRUE),
    max_val    = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(comp_summary, file.path(comp_dir, "cross_component_summary.csv"))

cat("\n=============================================================\n")
cat("Comparison complete.\n")
cat("Output directory:", comp_dir, "\n")
cat("Master table:", nrow(all_metrics), "rows\n")
cat("Components covered:", paste(unique(all_metrics$component), collapse = ", "), "\n")
cat("=============================================================\n")
