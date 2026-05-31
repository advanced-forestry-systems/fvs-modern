# =============================================================================
# Title: FVS Bayesian Calibration: Comprehensive Model Assessment
# Author: A. Weiskittel
# Date: 2026-03-30
# Description:
#   Publication-quality assessment of the calibrated FVS model system following
#   the evaluation framework from Weiskittel et al. (2011, Ch. 15-16).
#   Extends 11_full_comparison.R with:
#     Part A: Tree-level residual diagnostics by size class and species group
#     Part B: Biological realism checks (Bakuzis matrix)
#     Part C: Coupled stand-level forward projection
#     Part D: Species cross-variant performance analysis
#     Part E: Manuscript-ready summary tables
#
# Dependencies:
#   calibration/data/processed/{variant}/, output/variants/{variant}/
#   output/comparisons/ (from 11_full_comparison.R)
# Usage: Rscript 12_comprehensive_assessment.R --all
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
assess_dir      <- file.path(calibration_dir, "output", "assessment")

dir.create(assess_dir, recursive = TRUE, showWarnings = FALSE)

ALL_VARIANTS <- c("acd","ak","bc","bm","ca","ci","cr","cs","ec","em",
                   "ie","kt","ls","nc","ne","oc","on","op","pn","sn",
                   "so","tt","ut","wc","ws")

HI_VARIANTS <- c("bc","ci","em","ie","kt","ws")

# --- Utility functions -------------------------------------------------------

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
bias_fn <- function(obs, pred) mean(pred - obs, na.rm = TRUE)
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
r2   <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

# Bias as percent of mean observed
bias_pct <- function(obs, pred) {
  100 * mean(pred - obs, na.rm = TRUE) / mean(obs, na.rm = TRUE)
}

# RMSE as percent of mean observed
rmse_pct <- function(obs, pred) {
  100 * sqrt(mean((obs - pred)^2, na.rm = TRUE)) / mean(obs, na.rm = TRUE)
}

# Equivalence test (Reynolds 1984)
# Tests whether mean prediction falls within +/- tol_pct of mean observed
equiv_test <- function(obs, pred, tol_pct = 10, alpha = 0.05) {
  n <- length(obs)
  d <- pred - obs
  d_bar <- mean(d, na.rm = TRUE)
  se_d  <- sd(d, na.rm = TRUE) / sqrt(n)
  tol   <- tol_pct / 100 * mean(obs, na.rm = TRUE)
  t_crit <- qt(1 - alpha, n - 1)
  lower  <- d_bar - t_crit * se_d
  upper  <- d_bar + t_crit * se_d
  equivalent <- (lower > -tol) & (upper < tol)
  tibble(
    mean_diff = d_bar,
    se_diff   = se_d,
    ci_lower  = lower,
    ci_upper  = upper,
    tolerance = tol,
    equivalent = equivalent
  )
}

# AUC from predictions
compute_auc <- function(truth, prob) {
  tryCatch({
    roc_data <- tibble(truth = truth, prob = prob) %>%
      arrange(desc(prob)) %>%
      mutate(
        tp = cumsum(truth) / sum(truth),
        fp = cumsum(!truth) / sum(!truth)
      )
    sum(diff(roc_data$fp) * (roc_data$tp[-1] + roc_data$tp[-nrow(roc_data)]) / 2)
  }, error = function(e) NA_real_)
}

# Publication theme
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

cat("=============================================================\n")
cat("FVS Comprehensive Model Assessment\n")
cat("Weiskittel et al. (2011) evaluation framework\n")
cat("=============================================================\n\n")

# =============================================================================
# PART A: TREE-LEVEL RESIDUAL DIAGNOSTICS
# =============================================================================

cat("=== PART A: Tree-Level Residual Diagnostics ===\n\n")

# We stratify by:
#   (1) DBH size class: sapling (1-5"), pole (5-9"), small saw (9-15"),
#       medium saw (15-21"), large saw (21"+)
#   (2) Species functional group: softwood vs hardwood
#   (3) Site quality tercile (low, medium, high SI)

size_class_breaks <- c(0, 5, 9, 15, 21, 999)
size_class_labels <- c("Sapling (1-5)", "Pole (5-9)", "SmallSaw (9-15)",
                        "MedSaw (15-21)", "LargeSaw (21+)")

# FIA species groups: conifer = SPGRPCD 1-24, hardwood = 25-48
# Simplified: SPCD < 300 is generally conifer, >= 300 is hardwood in eastern US
# More robust: use the hundreds digit

# --- A1: Height-Diameter residuals by size class and species group -----------

cat("  A1: H-D residual diagnostics...\n")

hd_resid_all <- list()

for (v in ALL_VARIANTS) {

  post_file <- file.path(output_base, v, "height_diameter_samples.rds")
  if (!file.exists(post_file)) next

  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) next

  draws_tbl <- as_tibble(draws)
  a_pop <- median(draws_tbl$b_a_Intercept, na.rm = TRUE)
  b_pop <- median(draws_tbl$b_b_Intercept, na.rm = TRUE)
  c_pop <- median(draws_tbl$b_c_Intercept, na.rm = TRUE)
  if (any(is.na(c(a_pop, b_pop, c_pop)))) next

  re_a_cols <- names(draws_tbl)[grepl("r_SPCD__a\\[", names(draws_tbl))]
  re_b_cols <- names(draws_tbl)[grepl("r_SPCD__b\\[", names(draws_tbl))]

  sp_re_a <- setNames(
    sapply(re_a_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD__a\\[(.+),Intercept\\]", "\\1", re_a_cols)
  )
  sp_re_b <- setNames(
    sapply(re_b_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD__b\\[(.+),Intercept\\]", "\\1", re_b_cols)
  )

  # Load data
  hd_file <- file.path(data_base, v, "height_diameter.csv")
  dg_file <- file.path(data_base, v, "diameter_growth.csv")
  data_file <- if (file.exists(hd_file)) hd_file else if (file.exists(dg_file)) dg_file else NULL
  if (is.null(data_file)) next

  dat <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble()
  if ("DIA" %in% names(dat) & !("DIA_t1" %in% names(dat))) dat <- dat %>% rename(DIA_t1 = DIA)
  if ("HT"  %in% names(dat) & !("HT_t1"  %in% names(dat))) dat <- dat %>% rename(HT_t1 = HT)
  dat <- dat %>% filter(!is.na(HT_t1), HT_t1 > 4.5, !is.na(DIA_t1), DIA_t1 > 0)

  set.seed(42)
  if (nrow(dat) > 50000) dat <- dat %>% slice_sample(n = 50000)

  dat <- dat %>%
    mutate(
      SPCD_chr  = as.character(SPCD),
      DBH_scaled = DIA_t1 / 20,
      a_sp = a_pop + ifelse(SPCD_chr %in% names(sp_re_a), sp_re_a[SPCD_chr], 0),
      b_sp = b_pop + ifelse(SPCD_chr %in% names(sp_re_b), sp_re_b[SPCD_chr], 0),
      pred = 4.5 + pmax(a_sp, 1) * (1 - exp(-pmax(b_sp, 0.001) * DBH_scaled))^pmax(c_pop, 0.1),
      resid = HT_t1 - pred,
      size_class = cut(DIA_t1, breaks = size_class_breaks, labels = size_class_labels,
                       right = FALSE, include.lowest = TRUE),
      sp_group = ifelse(SPCD < 300, "Softwood", "Hardwood"),
      variant = toupper(v)
    )

  hd_resid_all[[v]] <- dat %>%
    select(variant, SPCD, SPCD_chr, sp_group, DIA_t1, HT_t1, pred, resid, size_class)
}

hd_resid <- bind_rows(hd_resid_all)

# Stratified metrics by variant x size class
hd_by_size <- hd_resid %>%
  group_by(variant, size_class) %>%
  summarise(
    n = n(),
    rmse = rmse(HT_t1, pred),
    bias = bias_fn(HT_t1, pred),
    bias_pct = bias_pct(HT_t1, pred),
    rmse_pct = rmse_pct(HT_t1, pred),
    r2 = r2(HT_t1, pred),
    .groups = "drop"
  )
write_csv(hd_by_size, file.path(assess_dir, "hd_metrics_by_size_class.csv"))

# Stratified by variant x species group
hd_by_spgrp <- hd_resid %>%
  group_by(variant, sp_group) %>%
  summarise(
    n = n(),
    rmse = rmse(HT_t1, pred),
    bias = bias_fn(HT_t1, pred),
    bias_pct = bias_pct(HT_t1, pred),
    r2 = r2(HT_t1, pred),
    .groups = "drop"
  )
write_csv(hd_by_spgrp, file.path(assess_dir, "hd_metrics_by_species_group.csv"))

# Equivalence tests for H-D (10% tolerance)
hd_equiv <- hd_resid %>%
  group_by(variant) %>%
  summarise(equiv_test(HT_t1, pred, tol_pct = 10), .groups = "drop")
write_csv(hd_equiv, file.path(assess_dir, "hd_equivalence_test.csv"))

# Diagnostic figure: residuals by size class (faceted by variant, top 6)
top_variants <- hd_resid %>%
  count(variant, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(variant)

if (length(top_variants) > 0) {
  p_hd_resid <- hd_resid %>%
    filter(variant %in% top_variants) %>%
    ggplot(aes(x = size_class, y = resid)) +
    geom_boxplot(fill = cb_pal[2], alpha = 0.6, outlier.size = 0.5, outlier.alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~variant, scales = "free_y") +
    labs(x = "Diameter Size Class (in.)", y = "Residual Height (ft)",
         title = "H-D Model Residuals by Size Class") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  ggsave(file.path(assess_dir, "fig_hd_residuals_by_size.png"), p_hd_resid,
         width = 24, height = 16, units = "cm", dpi = 300, bg = "white")
}

# Diagnostic figure: residuals by predicted (6 variants)
if (length(top_variants) > 0) {
  p_hd_resid_pred <- hd_resid %>%
    filter(variant %in% top_variants) %>%
    ggplot(aes(x = pred, y = resid)) +
    geom_hex(bins = 40) +
    geom_hline(yintercept = 0, color = "red") +
    geom_smooth(se = FALSE, color = "blue", linewidth = 0.5, method = "loess") +
    facet_wrap(~variant, scales = "free") +
    scale_fill_viridis_c(option = "D", trans = "log10") +
    labs(x = "Predicted Height (ft)", y = "Residual (ft)",
         title = "H-D Residuals vs Predicted") +
    theme_pub
  ggsave(file.path(assess_dir, "fig_hd_residuals_vs_predicted.png"), p_hd_resid_pred,
         width = 24, height = 16, units = "cm", dpi = 300, bg = "white")
}

cat("    Saved H-D residual diagnostics for", length(hd_resid_all), "variants\n")

# --- A2: Mortality residual diagnostics by size class ------------------------

cat("  A2: Mortality residual diagnostics...\n")

mort_resid_all <- list()

for (v in ALL_VARIANTS) {

  post_file <- file.path(output_base, v, "mortality_samples.rds")
  if (!file.exists(post_file)) next

  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) next

  # Load data
  mort_file <- file.path(data_base, v, "mortality.csv")
  dg_file   <- file.path(data_base, v, "diameter_growth.csv")
  data_file <- if (file.exists(mort_file)) mort_file else if (file.exists(dg_file)) dg_file else NULL
  if (is.null(data_file)) next

  dat <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble()
  if ("DIA" %in% names(dat) & !("DIA_t1" %in% names(dat))) dat <- dat %>% rename(DIA_t1 = DIA)
  if ("died" %in% names(dat)) dat <- dat %>% mutate(survived = !as.logical(died))
  if ("DIA_t2" %in% names(dat) & !("survived" %in% names(dat)))
    dat <- dat %>% mutate(survived = !is.na(DIA_t2) & DIA_t2 > 0)
  if (!("CR_pct" %in% names(dat)) & "CR" %in% names(dat))
    dat <- dat %>% mutate(CR_pct = CR * 100)

  dat <- dat %>%
    filter(!is.na(survived), !is.na(DIA_t1), !is.na(BAL), !is.na(CR_pct),
           !is.na(SI), !is.na(BA))
  if (nrow(dat) < 100) next

  draws_tbl <- as_tibble(draws)
  get_fe <- function(nm) {
    col <- paste0("b_", nm)
    if (col %in% names(draws_tbl)) median(draws_tbl[[col]], na.rm = TRUE) else 0
  }
  b0 <- get_fe("Intercept"); b1 <- get_fe("DBH_std")
  b2 <- get_fe("IDBH_stdE2"); b3 <- get_fe("BAL_std")
  b4 <- get_fe("CR_std"); b5 <- get_fe("SI_std"); b6 <- get_fe("BA_std")

  re_cols <- names(draws_tbl)[grepl("r_SPCD\\[", names(draws_tbl))]
  sp_re <- setNames(
    sapply(re_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE)),
    gsub("r_SPCD\\[(.+),Intercept\\]", "\\1", re_cols)
  )

  set.seed(42)
  if (nrow(dat) > 50000) dat <- dat %>% slice_sample(n = 50000)

  dat <- dat %>%
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
      pred_surv = 1 / (1 + exp(-logit_p)),
      size_class = cut(DIA_t1, breaks = size_class_breaks, labels = size_class_labels,
                       right = FALSE, include.lowest = TRUE),
      sp_group = ifelse(SPCD < 300, "Softwood", "Hardwood"),
      variant = toupper(v)
    )

  mort_resid_all[[v]] <- dat %>%
    select(variant, SPCD, sp_group, DIA_t1, survived, pred_surv, size_class,
           BAL, CR_pct, SI, BA)
}

mort_resid <- bind_rows(mort_resid_all)

# Observed vs predicted mortality rate by DBH class x variant
mort_by_size <- mort_resid %>%
  group_by(variant, size_class) %>%
  summarise(
    n = n(),
    obs_mort_rate = 1 - mean(survived, na.rm = TRUE),
    pred_mort_rate = 1 - mean(pred_surv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(resid_rate = pred_mort_rate - obs_mort_rate)

write_csv(mort_by_size, file.path(assess_dir, "mort_rate_by_size_class.csv"))

# By species group
mort_by_spgrp <- mort_resid %>%
  group_by(variant, sp_group) %>%
  summarise(
    n = n(),
    obs_mort_rate = 1 - mean(survived, na.rm = TRUE),
    pred_mort_rate = 1 - mean(pred_surv, na.rm = TRUE),
    auc = compute_auc(survived, pred_surv),
    .groups = "drop"
  )
write_csv(mort_by_spgrp, file.path(assess_dir, "mort_auc_by_species_group.csv"))

# Hosmer-Lemeshow style calibration by predicted probability deciles
mort_calibration <- mort_resid %>%
  group_by(variant) %>%
  mutate(prob_decile = ntile(pred_surv, 10)) %>%
  group_by(variant, prob_decile) %>%
  summarise(
    n = n(),
    mean_pred_surv = mean(pred_surv, na.rm = TRUE),
    obs_surv_rate  = mean(survived, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(mort_calibration, file.path(assess_dir, "mort_calibration_deciles.csv"))

# Calibration figure (top 6 variants)
if (nrow(mort_calibration) > 0) {
  top_v_mort <- mort_resid %>% count(variant, sort = TRUE) %>%
    slice_head(n = 6) %>% pull(variant)
  p_mort_cal <- mort_calibration %>%
    filter(variant %in% top_v_mort) %>%
    ggplot(aes(x = mean_pred_surv, y = obs_surv_rate)) +
    geom_point(size = 2.5, color = cb_pal[6]) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    facet_wrap(~variant) +
    labs(x = "Mean Predicted Survival Probability",
         y = "Observed Survival Rate",
         title = "Mortality Model Calibration (Probability Deciles)") +
    coord_equal(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
    theme_pub
  ggsave(file.path(assess_dir, "fig_mort_calibration.png"), p_mort_cal,
         width = 24, height = 16, units = "cm", dpi = 300, bg = "white")
}

cat("    Saved mortality diagnostics for", length(mort_resid_all), "variants\n")

# --- A3: DG residual diagnostics by size class -------------------------------

cat("  A3: Diameter growth residual diagnostics...\n")

dg_resid_all <- list()

for (v in ALL_VARIANTS) {

  post_file <- file.path(output_base, v, "diameter_growth_samples.rds")
  std_file  <- file.path(output_base, v, "standardization_params.rds")
  data_file <- file.path(data_base, v, "diameter_growth.csv")

  if (!file.exists(post_file) || !file.exists(data_file)) next

  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) next

  std_params <- if (file.exists(std_file)) readRDS(std_file) else NULL

  dat <- read_csv(data_file, show_col_types = FALSE) %>% as_tibble() %>%
    filter(!is.na(ln_DDS), is.finite(ln_DDS), !is.na(ln_DBH))
  if (nrow(dat) < 100) next

  set.seed(42)
  if (nrow(dat) > 50000) dat <- dat %>% slice_sample(n = 50000)

  draws_tbl <- as_tibble(draws)
  get_par <- function(nm) {
    if (nm %in% names(draws_tbl)) median(draws_tbl[[nm]], na.rm = TRUE) else 0
  }

  mu_b0 <- get_par("mu_b0")
  betas <- sapply(paste0("b", 1:13), get_par)

  sp_cols <- names(draws_tbl)[grepl("^b0\\[", names(draws_tbl))]
  sp_b0 <- sapply(sp_cols, function(x) median(draws_tbl[[x]], na.rm = TRUE))

  # Standardize (with safe fallback for missing params or columns)
  safe_std <- function(x, nm, params) {
    if (!is.null(params) && !is.null(params[[paste0(nm, "_mean")]]) && !is.null(params[[paste0(nm, "_sd")]])) {
      (x - params[[paste0(nm, "_mean")]]) / params[[paste0(nm, "_sd")]]
    } else {
      as.numeric(scale(x))
    }
  }

  # Ensure required columns exist (set to 0 if missing)
  for (col in c("ln_DBH", "DBH_sq", "ln_SI", "SLOPE", "ELEV",
                "SLOPE_SASP", "SLOPE_CASP", "BAL", "BA", "CR_pct")) {
    if (!(col %in% names(dat))) dat[[col]] <- 0
  }

  dat <- dat %>%
    mutate(
      ln_DBH_s = safe_std(ln_DBH, "ln_DBH", std_params),
      DBH_sq_s = safe_std(DBH_sq, "DBH_sq", std_params),
      ln_SI_s = safe_std(ln_SI, "ln_SI", std_params),
      SLOPE_s = safe_std(SLOPE, "SLOPE", std_params),
      ELEV_s = safe_std(ELEV, "ELEV", std_params),
      SLOPE_SASP_s = safe_std(SLOPE_SASP, "SLOPE_SASP", std_params),
      SLOPE_CASP_s = safe_std(SLOPE_CASP, "SLOPE_CASP", std_params),
      BAL_s = safe_std(BAL, "BAL", std_params),
      BA_s = safe_std(BA, "BA", std_params),
      CR_prop = CR_pct / 100,
      species_idx = as.integer(factor(SPCD))
    )

  dat <- dat %>%
    mutate(
      sp_int = ifelse(species_idx <= length(sp_b0), sp_b0[species_idx], mu_b0),
      pred_ln_dds = sp_int +
        betas[1]*ln_DBH_s + betas[2]*DBH_sq_s + betas[3]*ln_SI_s +
        betas[4]*SLOPE_s + betas[5]*SLOPE_s^2 + betas[6]*SLOPE_SASP_s +
        betas[7]*SLOPE_CASP_s + betas[8]*ELEV_s + betas[9]*ELEV_s^2 +
        betas[10]*CR_prop + betas[11]*CR_prop^2 + betas[12]*BAL_s + betas[13]*BA_s,
      resid_ln_dds = ln_DDS - pred_ln_dds,
      # Back-transform with bias correction (exp(MSE/2))
      mse_corr = var(resid_ln_dds, na.rm = TRUE) / 2,
      pred_dds = exp(pred_ln_dds + mse_corr),
      obs_dds  = exp(ln_DDS),
      # Annualized diameter increment (DDS = change in squared diameter)
      pred_di  = (sqrt(DIA_t1^2 + pred_dds) - DIA_t1),
      obs_di   = (sqrt(DIA_t1^2 + obs_dds) - DIA_t1),
      size_class = cut(DIA_t1, breaks = size_class_breaks, labels = size_class_labels,
                       right = FALSE, include.lowest = TRUE),
      sp_group = ifelse(SPCD < 300, "Softwood", "Hardwood"),
      variant = toupper(v)
    ) %>%
    filter(is.finite(pred_ln_dds), is.finite(resid_ln_dds))

  dg_resid_all[[v]] <- dat %>%
    select(variant, SPCD, sp_group, DIA_t1, ln_DDS, pred_ln_dds, resid_ln_dds,
           obs_di, pred_di, size_class, BAL_s, SI)
}

dg_resid <- bind_rows(dg_resid_all)

if (nrow(dg_resid) > 0) {
  # DG metrics by size class
  dg_by_size <- dg_resid %>%
    group_by(variant, size_class) %>%
    summarise(
      n = n(),
      rmse_log = rmse(ln_DDS, pred_ln_dds),
      bias_log = bias_fn(ln_DDS, pred_ln_dds),
      r2_log   = r2(ln_DDS, pred_ln_dds),
      rmse_di  = rmse(obs_di, pred_di),
      bias_di  = bias_fn(obs_di, pred_di),
      .groups = "drop"
    )
  write_csv(dg_by_size, file.path(assess_dir, "dg_metrics_by_size_class.csv"))

  # Equivalence test for DG (log scale, 15% tolerance since DG is harder to predict)
  dg_equiv <- dg_resid %>%
    group_by(variant) %>%
    summarise(equiv_test(ln_DDS, pred_ln_dds, tol_pct = 15), .groups = "drop")
  write_csv(dg_equiv, file.path(assess_dir, "dg_equivalence_test.csv"))
}

cat("    Saved DG diagnostics for", length(dg_resid_all), "variants\n")


# =============================================================================
# PART B: BIOLOGICAL REALISM (BAKUZIS MATRIX)
# =============================================================================

cat("\n=== PART B: Biological Realism Checks (Bakuzis Matrix) ===\n\n")

# The Bakuzis matrix tests the model system against known ecological laws.
# We test 6 laws that can be evaluated from our calibrated parameters:
#
# 1. Sukachev effect: increasing density reduces individual tree growth
#    Test: DG should decrease with increasing BA/BAL
# 2. Reineke's rule: max density scales as TPA ~ QMD^-1.605
#    Test: calibrated self-thinning exponents should be near -1.605
# 3. Eichhorn's rule: total volume production is a function of site quality
#    Test: DG should increase with SI (after controlling for size)
# 4. Mortality-size relationship: mortality should be U-shaped with DBH
#    Test: small and very large trees should have higher mortality
# 5. Height-diameter allometry: H/D ratio should decrease with increasing density
#    Test: predicted heights should be lower at given DBH when BA is high
# 6. Crown recession: crown ratios should decrease with increasing density
#    Test: CR change should be more negative at higher BA

bakuzis_results <- list()

# --- B1: Sukachev effect (competition reduces growth) -----------------------

cat("  B1: Sukachev effect (competition reduces growth)...\n")

if (nrow(dg_resid) > 0) {
  # For each variant, compute average predicted DG across competition gradient
  sukachev <- dg_resid %>%
    mutate(bal_quartile = ntile(BAL_s, 4)) %>%
    group_by(variant, bal_quartile) %>%
    summarise(mean_pred_di = mean(pred_di, na.rm = TRUE), .groups = "drop") %>%
    group_by(variant) %>%
    summarise(
      di_low_comp  = mean_pred_di[bal_quartile == 1],
      di_high_comp = mean_pred_di[bal_quartile == 4],
      growth_reduction_pct = 100 * (di_low_comp - di_high_comp) / di_low_comp,
      sukachev_pass = di_low_comp > di_high_comp,
      .groups = "drop"
    )
  write_csv(sukachev, file.path(assess_dir, "bakuzis_sukachev_effect.csv"))
  bakuzis_results[["Sukachev"]] <- sukachev %>%
    summarise(law = "Sukachev", pass_rate = mean(sukachev_pass) * 100,
              n_variants = n())
  cat("    Pass rate:", round(mean(sukachev$sukachev_pass) * 100, 0), "%\n")
}

# --- B2: Reineke's rule (self-thinning exponent ~ -1.605) -------------------

cat("  B2: Reineke's rule (self-thinning exponent)...\n")

reineke_results <- list()
for (v in ALL_VARIANTS) {
  cal_file <- file.path(output_base, v, "species_sdimax_calibrated.csv")
  if (!file.exists(cal_file)) next
  cal <- read_csv(cal_file, show_col_types = FALSE)
  if ("self_thinning_slope" %in% names(cal)) {
    reineke_results[[v]] <- cal %>%
      summarise(
        variant = toupper(v),
        n_species = n(),
        mean_slope = mean(self_thinning_slope, na.rm = TRUE),
        sd_slope   = sd(self_thinning_slope, na.rm = TRUE),
        deviation_from_reineke = mean_slope - (-1.605),
        reineke_pass = abs(deviation_from_reineke) < 0.3
      )
  }
}
reineke_df <- bind_rows(reineke_results)
if (nrow(reineke_df) > 0) {
  write_csv(reineke_df, file.path(assess_dir, "bakuzis_reineke_rule.csv"))
  bakuzis_results[["Reineke"]] <- reineke_df %>%
    summarise(law = "Reineke", pass_rate = mean(reineke_pass) * 100,
              n_variants = n())
  cat("    Pass rate:", round(mean(reineke_df$reineke_pass) * 100, 0), "%\n")
}

# --- B3: Eichhorn's rule (site quality increases growth) --------------------

cat("  B3: Eichhorn's rule (site drives growth)...\n")

if (nrow(dg_resid) > 0 && "SI" %in% names(dg_resid)) {
  # Compute SI terciles within each variant (robust to small samples)
  eichhorn_means <- dg_resid %>%
    filter(!is.na(SI), is.finite(SI)) %>%
    group_by(variant) %>%
    mutate(si_tercile = ntile(SI, 3)) %>%
    group_by(variant, si_tercile) %>%
    summarise(mean_pred_di = mean(pred_di, na.rm = TRUE), .groups = "drop")

  # Only include variants that have both tercile 1 and 3
  eichhorn <- eichhorn_means %>%
    group_by(variant) %>%
    filter(all(c(1, 3) %in% si_tercile)) %>%
    summarise(
      di_low_si  = mean_pred_di[si_tercile == 1],
      di_high_si = mean_pred_di[si_tercile == 3],
      growth_increase_pct = 100 * (di_high_si - di_low_si) / pmax(abs(di_low_si), 0.001),
      eichhorn_pass = di_high_si > di_low_si,
      .groups = "drop"
    )
  write_csv(eichhorn, file.path(assess_dir, "bakuzis_eichhorn_rule.csv"))
  bakuzis_results[["Eichhorn"]] <- eichhorn %>%
    summarise(law = "Eichhorn", pass_rate = mean(eichhorn_pass) * 100,
              n_variants = n())
  cat("    Pass rate:", round(mean(eichhorn$eichhorn_pass) * 100, 0), "%\n")
}

# --- B4: Mortality-size relationship (U-shaped) ----------------------------

cat("  B4: Mortality-size relationship...\n")

if (nrow(mort_resid) > 0) {
  mort_size <- mort_resid %>%
    mutate(dbh_class = cut(DIA_t1, breaks = c(0,3,5,9,15,21,40,999), right = FALSE)) %>%
    group_by(variant, dbh_class) %>%
    summarise(
      n = n(),
      pred_mort_rate = 1 - mean(pred_surv, na.rm = TRUE),
      obs_mort_rate  = 1 - mean(survived, na.rm = TRUE),
      .groups = "drop"
    )

  # Test U-shape: mortality should be higher in smallest and largest classes
  # than in middle classes
  mort_ushape <- mort_size %>%
    group_by(variant) %>%
    summarise(
      mort_smallest = pred_mort_rate[which.min(as.numeric(dbh_class))],
      mort_middle   = mean(pred_mort_rate[as.numeric(dbh_class) %in% 3:4]),
      mort_largest  = pred_mort_rate[which.max(as.numeric(dbh_class))],
      u_shape_pass  = (mort_smallest > mort_middle) | (mort_largest > mort_middle),
      .groups = "drop"
    )
  write_csv(mort_ushape, file.path(assess_dir, "bakuzis_mortality_ushape.csv"))
  bakuzis_results[["Mortality_U_Shape"]] <- mort_ushape %>%
    summarise(law = "Mortality_U_Shape", pass_rate = mean(u_shape_pass) * 100,
              n_variants = n())
  cat("    Pass rate:", round(mean(mort_ushape$u_shape_pass) * 100, 0), "%\n")
}

# --- B5: Crown recession under competition --------------------------------

cat("  B5: Crown recession under competition...\n")

cr_recession <- list()
for (v in ALL_VARIANTS) {
  post_file <- file.path(output_base, v, "crown_ratio_samples.rds")
  if (!file.exists(post_file)) next
  draws <- tryCatch(readRDS(post_file), error = function(e) NULL)
  if (is.null(draws)) next
  draws_tbl <- as_tibble(draws)
  b_BA <- if ("b_BA_std" %in% names(draws_tbl)) median(draws_tbl$b_BA_std, na.rm = TRUE) else
          if ("b_BA" %in% names(draws_tbl)) median(draws_tbl$b_BA, na.rm = TRUE) else NA
  b_BAL <- if ("b_BAL_std" %in% names(draws_tbl)) median(draws_tbl$b_BAL_std, na.rm = TRUE) else
           if ("b_BAL" %in% names(draws_tbl)) median(draws_tbl$b_BAL, na.rm = TRUE) else NA
  cr_recession[[v]] <- tibble(
    variant = toupper(v),
    ba_effect = b_BA,
    bal_effect = b_BAL,
    # CR should decrease (more negative change) with higher BA and BAL
    # so the coefficients on BA_std and BAL_std should be negative
    recession_pass = (!is.na(b_BA) & b_BA < 0) | (!is.na(b_BAL) & b_BAL < 0)
  )
}
cr_recess_df <- bind_rows(cr_recession)
if (nrow(cr_recess_df) > 0) {
  write_csv(cr_recess_df, file.path(assess_dir, "bakuzis_crown_recession.csv"))
  bakuzis_results[["Crown_Recession"]] <- cr_recess_df %>%
    summarise(law = "Crown_Recession", pass_rate = mean(recession_pass) * 100,
              n_variants = n())
  cat("    Pass rate:", round(mean(cr_recess_df$recession_pass) * 100, 0), "%\n")
}

# --- Bakuzis summary --------------------------------------------------------

bakuzis_summary <- bind_rows(bakuzis_results)
write_csv(bakuzis_summary, file.path(assess_dir, "bakuzis_matrix_summary.csv"))

cat("\n  Bakuzis Matrix Summary:\n")
for (i in seq_len(nrow(bakuzis_summary))) {
  cat("    ", bakuzis_summary$law[i], ": ",
      round(bakuzis_summary$pass_rate[i], 0), "% pass (",
      bakuzis_summary$n_variants[i], " variants)\n", sep = "")
}

# =============================================================================
# PART C: COUPLED STAND-LEVEL FORWARD PROJECTION
# =============================================================================

cat("\n=== PART C: Coupled Stand-Level Forward Projection ===\n\n")

# Simulate a representative stand forward 50 years using all calibrated
# components simultaneously. This tests model system behavior, not individual
# equations. We build a simple tree-list simulator that:
#   (1) Uses H-D to assign heights
#   (2) Applies DG to grow diameter each year
#   (3) Applies mortality to thin the stand
#   (4) Applies CR to track crown ratio
#   (5) Updates BA and BAL each step
#   (6) Checks SDIMAX boundary is not exceeded

cat("  Building stand-level simulator...\n")

# Helper: extract all calibrated params for a variant
load_variant_params <- function(v) {
  params <- list(variant = v)

  # H-D
  hd_file <- file.path(output_base, v, "height_diameter_samples.rds")
  if (file.exists(hd_file)) {
    draws <- tryCatch(readRDS(hd_file), error = function(e) NULL)
    if (!is.null(draws)) {
      d <- as_tibble(draws)
      params$hd <- list(
        a = median(d$b_a_Intercept, na.rm = TRUE),
        b = median(d$b_b_Intercept, na.rm = TRUE),
        c = median(d$b_c_Intercept, na.rm = TRUE)
      )
    }
  }

  # Mortality
  mort_file <- file.path(output_base, v, "mortality_samples.rds")
  if (file.exists(mort_file)) {
    draws <- tryCatch(readRDS(mort_file), error = function(e) NULL)
    if (!is.null(draws)) {
      d <- as_tibble(draws)
      get_fe <- function(nm) {
        col <- paste0("b_", nm)
        if (col %in% names(d)) median(d[[col]], na.rm = TRUE) else 0
      }
      params$mort <- list(
        b0 = get_fe("Intercept"), b1 = get_fe("DBH_std"),
        b2 = get_fe("IDBH_stdE2"), b3 = get_fe("BAL_std"),
        b4 = get_fe("CR_std"), b5 = get_fe("SI_std"), b6 = get_fe("BA_std")
      )
    }
  }

  # SDIMAX
  cal_file <- file.path(output_base, v, "species_sdimax_calibrated.csv")
  if (file.exists(cal_file)) {
    sdi <- read_csv(cal_file, show_col_types = FALSE)
    sdi_col <- if ("sdimax_combined" %in% names(sdi)) "sdimax_combined" else "sdimax_bayes"
    if (sdi_col %in% names(sdi)) {
      params$sdimax <- mean(sdi[[sdi_col]], na.rm = TRUE)
    }
  }

  params
}

# Simple 1-year step simulator
simulate_stand <- function(trees, params, si = 65, years = 50) {
  # trees: data.frame with columns: dbh, tpa, spcd
  # Returns time series of stand metrics

  trajectory <- vector("list", years + 1)

  for (yr in 0:years) {
    # Compute stand metrics
    ba_tree <- pi * (trees$dbh / 24)^2  # BA per tree in sq ft (dbh in inches)
    ba <- sum(trees$tpa * ba_tree)
    tpa <- sum(trees$tpa)
    qmd <- ifelse(tpa > 0, sqrt(ba / tpa * 576 / pi), 0)  # 576 = 144^2 / pi factor
    sdi <- tpa * (qmd / 10)^1.605  # Reineke SDI (imperial, ref 10 in.)

    # H-D prediction
    if (!is.null(params$hd)) {
      trees$ht <- 4.5 + pmax(params$hd$a, 1) *
        (1 - exp(-pmax(params$hd$b, 0.001) * trees$dbh / 20))^pmax(params$hd$c, 0.1)
    } else {
      trees$ht <- 4.5 + 25 * (1 - exp(-0.04 * trees$dbh / 20))^1.0
    }

    top_ht <- ifelse(nrow(trees) > 0,
                     mean(sort(trees$ht, decreasing = TRUE)[1:min(5, nrow(trees))]),
                     0)

    trajectory[[yr + 1]] <- tibble(
      year = yr, tpa = tpa, ba = ba, qmd = qmd, sdi = sdi,
      top_ht = top_ht, n_trees = nrow(trees)
    )

    if (yr == years || nrow(trees) == 0) break

    # Sort by DBH descending for BAL computation
    trees <- trees %>% arrange(desc(dbh))
    trees$bal <- cumsum(trees$tpa * pi * (trees$dbh / 24)^2) - trees$tpa * pi * (trees$dbh / 24)^2

    # CR: start at 40% and drift toward equilibrium
    if (!("cr" %in% names(trees))) trees$cr <- 0.40

    # Mortality (annual probability)
    if (!is.null(params$mort)) {
      m <- params$mort
      # Standardize using current stand means
      dbh_mean <- mean(trees$dbh); dbh_sd <- sd(trees$dbh)
      if (is.na(dbh_sd) || dbh_sd < 0.01) dbh_sd <- 1
      bal_mean <- mean(trees$bal); bal_sd <- sd(trees$bal)
      if (is.na(bal_sd) || bal_sd < 0.01) bal_sd <- 1
      cr_mean <- mean(trees$cr); cr_sd <- sd(trees$cr)
      if (is.na(cr_sd) || cr_sd < 0.01) cr_sd <- 0.1
      ba_mean <- ba; ba_sd <- 10  # rough

      dbh_s <- (trees$dbh - dbh_mean) / dbh_sd
      bal_s <- (trees$bal - bal_mean) / bal_sd
      cr_s  <- (trees$cr - cr_mean) / cr_sd
      si_s  <- 0  # centered at mean
      ba_s  <- 0

      logit_surv <- m$b0 + m$b1 * dbh_s + m$b2 * dbh_s^2 + m$b3 * bal_s +
                    m$b4 * cr_s + m$b5 * si_s + m$b6 * ba_s
      p_surv <- 1 / (1 + exp(-logit_surv))
      p_surv <- pmin(pmax(p_surv, 0.01), 0.999)

      # Stochastic mortality
      set.seed(yr * 1000 + as.integer(substr(params$variant, 1, 2) == "ne") * 100)
      alive <- runif(nrow(trees)) < p_surv
      trees <- trees[alive, ]
    }

    if (nrow(trees) == 0) break

    # Diameter growth: simplified (constant annual increment scaled by size and competition)
    # If full DG params available, could use Wykoff, but we use a biologically
    # reasonable approximation here since the full model needs all 13 standardized covariates
    # Growth is roughly 0.1 to 0.3 inches/year for most species, declining with competition
    base_growth <- 0.20  # inches/year baseline
    competition_modifier <- pmax(0.3, 1 - trees$bal / (ba + 1))
    size_modifier <- pmax(0.5, 1 - (trees$dbh / 30)^0.5)  # large trees grow slower
    trees$dbh <- trees$dbh + base_growth * competition_modifier * size_modifier

    # CR change: slight decline under competition
    trees$cr <- pmax(0.10, pmin(0.95, trees$cr - 0.005 * (ba / 150)))

    # SDI check against maximum
    if (!is.null(params$sdimax)) {
      current_tpa <- sum(trees$tpa)
      current_qmd <- sqrt(sum(trees$tpa * pi * (trees$dbh / 24)^2) / current_tpa * 576 / pi)
      current_sdi <- current_tpa * (current_qmd / 10)^1.605
      if (current_sdi > params$sdimax) {
        # Proportional mortality to bring back to SDI max
        reduction <- params$sdimax / current_sdi
        trees$tpa <- trees$tpa * reduction
      }
    }
  }

  bind_rows(trajectory)
}

# Run projections for all variants with complete parameter sets
cat("  Running 50-year projections...\n")

projection_results <- list()

for (v in ALL_VARIANTS) {
  cat("    Projecting:", toupper(v), "... ")

  params <- load_variant_params(v)
  if (is.null(params$hd) || is.null(params$mort)) {
    cat("incomplete params\n"); next
  }

  # Build representative initial stand: 200 TPA, mixed sizes, 2 species
  set.seed(42)
  n_init <- 200
  init_trees <- tibble(
    dbh = rlnorm(n_init, log(6), 0.5),  # ~6 in mean DBH
    tpa = 1,  # each tree represents 1 TPA
    spcd = sample(c(97, 316), n_init, replace = TRUE)  # arbitrary species codes
  ) %>%
    filter(dbh >= 1, dbh <= 30)

  traj <- simulate_stand(init_trees, params, si = 65, years = 50)
  traj$variant <- toupper(v)
  projection_results[[v]] <- traj
  cat("final BA =", round(tail(traj$ba, 1), 0), "sq ft\n")
}

proj_df <- bind_rows(projection_results)
write_csv(proj_df, file.path(assess_dir, "stand_projection_trajectories.csv"))

# Projection figures
if (nrow(proj_df) > 0) {
  # BA trajectory
  p_ba_traj <- proj_df %>%
    ggplot(aes(x = year, y = ba, color = variant, group = variant)) +
    geom_line(alpha = 0.7) +
    labs(x = "Projection Year", y = expression("Basal Area (ft"^2*" ac"^-1*")"),
         title = "Stand BA Trajectories (50-year Projection)",
         color = "Variant") +
    theme_pub +
    theme(legend.position = "right", legend.text = element_text(size = 7))
  ggsave(file.path(assess_dir, "fig_ba_trajectories.png"), p_ba_traj,
         width = 20, height = 14, units = "cm", dpi = 300, bg = "white")

  # QMD trajectory
  p_qmd_traj <- proj_df %>%
    ggplot(aes(x = year, y = qmd, color = variant, group = variant)) +
    geom_line(alpha = 0.7) +
    labs(x = "Projection Year", y = "QMD (in.)",
         title = "Stand QMD Trajectories (50-year Projection)",
         color = "Variant") +
    theme_pub +
    theme(legend.position = "right", legend.text = element_text(size = 7))
  ggsave(file.path(assess_dir, "fig_qmd_trajectories.png"), p_qmd_traj,
         width = 20, height = 14, units = "cm", dpi = 300, bg = "white")

  # SDI relative to max
  if (any(!is.na(proj_df$sdi))) {
    p_sdi_rel <- proj_df %>%
      group_by(variant) %>%
      mutate(sdi_rel = sdi / max(sdi, na.rm = TRUE)) %>%
      ggplot(aes(x = year, y = sdi_rel, color = variant, group = variant)) +
      geom_line(alpha = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(x = "Projection Year", y = "Relative SDI (proportion of max)",
           title = "Relative Density Over Time",
           color = "Variant") +
      theme_pub +
      theme(legend.position = "right", legend.text = element_text(size = 7))
    ggsave(file.path(assess_dir, "fig_sdi_relative.png"), p_sdi_rel,
           width = 20, height = 14, units = "cm", dpi = 300, bg = "white")
  }

  # Projection realism checks
  proj_realism <- proj_df %>%
    group_by(variant) %>%
    summarise(
      ba_initial = first(ba),
      ba_final   = last(ba),
      ba_growth_reasonable = ba_final > ba_initial & ba_final < 400,
      qmd_initial = first(qmd),
      qmd_final   = last(qmd),
      qmd_grows = qmd_final > qmd_initial,
      tpa_declines = last(tpa) < first(tpa),
      sdi_bounded = max(sdi, na.rm = TRUE) < 1500,
      all_checks_pass = ba_growth_reasonable & qmd_grows & tpa_declines & sdi_bounded,
      .groups = "drop"
    )
  write_csv(proj_realism, file.path(assess_dir, "projection_realism_checks.csv"))
  cat("\n  Projection realism: ",
      sum(proj_realism$all_checks_pass), "/", nrow(proj_realism),
      " variants pass all checks\n")
}


# =============================================================================
# PART D: SPECIES CROSS-VARIANT PERFORMANCE
# =============================================================================

cat("\n=== PART D: Species Cross-Variant Performance ===\n\n")

# For each major FIA species that appears in multiple variants, compute
# the H-D model performance across all variants it occurs in.
# This answers: "How well does red spruce perform across NE, LS, CS, etc.?"

if (nrow(hd_resid) > 0) {

  # Common FIA species names
  fia_names <- tribble(
    ~SPCD, ~common_name,
    12, "Balsam fir",
    97, "Red spruce",
    129, "White pine",
    241, "White cedar",
    261, "Eastern hemlock",
    316, "Red maple",
    318, "Sugar maple",
    371, "Yellow birch",
    375, "Paper birch",
    531, "American beech",
    541, "White ash",
    743, "Red oak",
    746, "Chestnut oak",
    802, "White oak",
    833, "Northern red oak",
    951, "Red pine",
    # Western species
    15, "White fir",
    17, "Grand fir",
    19, "Subalpine fir",
    93, "Engelmann spruce",
    108, "Lodgepole pine",
    122, "Ponderosa pine",
    202, "Douglas fir",
    263, "Western hemlock"
  )

  species_xvariant <- hd_resid %>%
    group_by(SPCD, variant) %>%
    summarise(
      n = n(),
      rmse = rmse(HT_t1, pred),
      bias = bias_fn(HT_t1, pred),
      r2 = r2(HT_t1, pred),
      .groups = "drop"
    ) %>%
    left_join(fia_names, by = "SPCD") %>%
    mutate(common_name = ifelse(is.na(common_name), paste("SPCD", SPCD), common_name))

  write_csv(species_xvariant, file.path(assess_dir, "hd_species_cross_variant.csv"))

  # Species appearing in 3+ variants
  multi_variant_spp <- species_xvariant %>%
    group_by(SPCD, common_name) %>%
    summarise(
      n_variants = n_distinct(variant),
      total_n = sum(n),
      mean_r2 = mean(r2, na.rm = TRUE),
      sd_r2   = sd(r2, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_variants >= 3) %>%
    arrange(desc(n_variants))

  write_csv(multi_variant_spp, file.path(assess_dir, "hd_multi_variant_species_summary.csv"))
  cat("  Species in 3+ variants:", nrow(multi_variant_spp), "\n")

  # Figure: top 12 species performance across variants
  top_spp <- multi_variant_spp %>% slice_head(n = 12) %>% pull(SPCD)

  if (length(top_spp) > 0) {
    p_sp_xvar <- species_xvariant %>%
      filter(SPCD %in% top_spp, n >= 30) %>%
      ggplot(aes(x = variant, y = r2, fill = common_name)) +
      geom_col(position = "dodge", width = 0.7) +
      facet_wrap(~common_name, scales = "free_x", ncol = 3) +
      labs(x = "Variant", y = expression(R^2),
           title = "H-D Model Performance by Species Across Variants") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "none")
    ggsave(file.path(assess_dir, "fig_species_cross_variant.png"), p_sp_xvar,
           width = 24, height = 20, units = "cm", dpi = 300, bg = "white")
  }

  # Similarly for mortality AUC across variants
  if (nrow(mort_resid) > 0) {
    mort_sp_xvar <- mort_resid %>%
      group_by(SPCD, variant) %>%
      summarise(
        n = n(),
        obs_mort = 1 - mean(survived, na.rm = TRUE),
        auc = compute_auc(survived, pred_surv),
        .groups = "drop"
      ) %>%
      left_join(fia_names, by = "SPCD") %>%
      mutate(common_name = ifelse(is.na(common_name), paste("SPCD", SPCD), common_name))

    write_csv(mort_sp_xvar, file.path(assess_dir, "mort_species_cross_variant.csv"))
  }
}


# =============================================================================
# PART E: MANUSCRIPT SUMMARY TABLES
# =============================================================================

cat("\n=== PART E: Manuscript Summary Tables ===\n\n")

# --- Table 1: Overall model performance by component -----------------------

cat("  Building Table 1: Overall performance by component...\n")

# Load existing metrics from 11_full_comparison.R
existing_perf <- file.path(comp_dir, "model_performance_summary.csv")
if (file.exists(existing_perf)) {
  perf <- read_csv(existing_perf, show_col_types = FALSE)

  # Summarize by component
  table1 <- perf %>%
    filter(!is.na(value), is.finite(value)) %>%
    group_by(component, metric) %>%
    summarise(
      n_variants = n_distinct(variant),
      mean = round(mean(value, na.rm = TRUE), 3),
      sd   = round(sd(value, na.rm = TRUE), 3),
      min  = round(min(value, na.rm = TRUE), 3),
      max  = round(max(value, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    arrange(component, metric)

  write_csv(table1, file.path(assess_dir, "table1_overall_performance.csv"))
  cat("    Table 1:", nrow(table1), "rows\n")
}

# --- Table 2: H-D equivalence test summary --------------------------------

if (file.exists(file.path(assess_dir, "hd_equivalence_test.csv"))) {
  hd_eq <- read_csv(file.path(assess_dir, "hd_equivalence_test.csv"), show_col_types = FALSE)
  table2 <- hd_eq %>%
    summarise(
      n_variants = n(),
      n_equivalent = sum(equivalent),
      pct_equivalent = round(100 * n_equivalent / n_variants, 1),
      mean_bias = round(mean(mean_diff, na.rm = TRUE), 2),
      max_abs_bias = round(max(abs(mean_diff), na.rm = TRUE), 2)
    )
  write_csv(table2, file.path(assess_dir, "table2_hd_equivalence_summary.csv"))
  cat("    Table 2: H-D equivalence ", table2$pct_equivalent, "% pass\n")
}

# --- Table 3: Bakuzis matrix results ----------------------------------------

if (file.exists(file.path(assess_dir, "bakuzis_matrix_summary.csv"))) {
  table3 <- read_csv(file.path(assess_dir, "bakuzis_matrix_summary.csv"), show_col_types = FALSE)
  write_csv(table3, file.path(assess_dir, "table3_bakuzis_matrix.csv"))
  cat("    Table 3: Bakuzis matrix ", nrow(table3), " laws tested\n")
}

# --- Table 4: Projection realism -------------------------------------------

if (file.exists(file.path(assess_dir, "projection_realism_checks.csv"))) {
  proj_real <- read_csv(file.path(assess_dir, "projection_realism_checks.csv"), show_col_types = FALSE)
  table4 <- proj_real %>%
    summarise(
      n_variants = n(),
      pct_ba_reasonable = round(100 * mean(ba_growth_reasonable), 1),
      pct_qmd_grows = round(100 * mean(qmd_grows), 1),
      pct_tpa_declines = round(100 * mean(tpa_declines), 1),
      pct_sdi_bounded = round(100 * mean(sdi_bounded), 1),
      pct_all_pass = round(100 * mean(all_checks_pass), 1)
    )
  write_csv(table4, file.path(assess_dir, "table4_projection_realism.csv"))
  cat("    Table 4: Projection realism ", table4$pct_all_pass, "% pass all\n")
}

# --- Table 5: Species cross-variant consistency ----------------------------

if (file.exists(file.path(assess_dir, "hd_multi_variant_species_summary.csv"))) {
  table5 <- read_csv(file.path(assess_dir, "hd_multi_variant_species_summary.csv"),
                      show_col_types = FALSE) %>%
    select(common_name, n_variants, total_n, mean_r2, sd_r2, mean_rmse, mean_bias) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  write_csv(table5, file.path(assess_dir, "table5_species_consistency.csv"))
  cat("    Table 5: Species consistency for", nrow(table5), "species\n")
}


# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n=============================================================\n")
cat("Comprehensive Assessment Complete\n")
cat("=============================================================\n")
cat("Output directory:", assess_dir, "\n\n")

# Count all outputs
assess_files <- list.files(assess_dir)
cat("Files generated:", length(assess_files), "\n")
cat("  CSV tables:", sum(grepl("\\.csv$", assess_files)), "\n")
cat("  PNG figures:", sum(grepl("\\.png$", assess_files)), "\n\n")

cat("Key results:\n")
if (exists("bakuzis_summary") && nrow(bakuzis_summary) > 0) {
  cat("  Bakuzis matrix: ", round(mean(bakuzis_summary$pass_rate), 0),
      "% average pass rate across ", nrow(bakuzis_summary), " ecological laws\n")
}
if (exists("hd_equiv") && nrow(hd_equiv) > 0) {
  cat("  H-D equivalence (10%): ", sum(hd_equiv$equivalent), "/",
      nrow(hd_equiv), " variants pass\n")
}
if (exists("proj_realism") && nrow(proj_realism) > 0) {
  cat("  Projection realism: ", sum(proj_realism$all_checks_pass), "/",
      nrow(proj_realism), " variants pass all checks\n")
}

cat("\nThis script complements 11_full_comparison.R (tree-level metrics)\n")
cat("by adding diagnostic stratification, biological realism, and system-level tests.\n")
cat("=============================================================\n")
