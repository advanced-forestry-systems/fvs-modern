#!/usr/bin/env Rscript

################################################################################
# FIA Benchmark Engine v2: Tree-Level Bayesian Projection Benchmark
################################################################################
# Purpose: Benchmark calibrated FVS parameters against actual FIA remeasurement
#          data using tree-level projection with Bayesian posteriors.
#
# Key improvements over v1:
#   - Tree-level projection using actual Bayesian posterior medians
#   - DG: Wykoff ln(DDS) model with 13 covariates (calibrated)
#   - Mortality: Logistic survival model with species random effects
#   - H-D: Chapman-Richards with species random effects on a and b
#   - Proper plot pair linkage via PLT_CN / PREV_PLT_CN
#   - Expected-value mortality (no stochastic noise)
#   - Interval scaling: models trained on mean interval, applied to actual
#
# Author: A. Weiskittel
# Date: 2026-04-04
################################################################################

library(data.table)
library(ggplot2)

# Setup ========================================================================
cat("FIA Benchmark Engine v2: Tree-Level Bayesian Projection\n")
cat(strrep("=", 80), "\n\n")

project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (project_root == "") project_root <- getwd()

fia_root    <- "/users/PUOM0008/crsfaaron/FIA"
calib_root  <- file.path(project_root, "calibration/output/variants")
data_root   <- file.path(project_root, "calibration/data/processed")
output_root <- file.path(project_root, "calibration/output/comparisons")

dir.create(file.path(output_root, "manuscript_tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "manuscript_figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "intermediate"), showWarnings = FALSE, recursive = TRUE)

cat("Project root:", project_root, "\n")
cat("FIA data path:", fia_root, "\n\n")

ALL_VARIANTS <- c("acd", "ak", "bc", "bm", "ca", "ci", "cr", "cs", "ec",
                   "em", "ie", "kt", "ls", "nc", "ne", "oc", "on", "op",
                   "pn", "sn", "so", "tt", "ut", "wc", "ws")

# Mortality strategy:
#   "full"       = all fixed+random effects with training-data standardization
#   "re_only"    = species RE + intercept only (no covariates)
#   "within"     = within-condition standardization for DBH/BAL/CR, SI/BA zeroed
#   "no_dbh_si"  = training-data std for BAL/CR/BA, DBH and SI zeroed
#   "no_si"      = training-data std for DBH/BAL/CR/BA, only SI zeroed (captures U-shape)
MORT_STRATEGY <- "no_si"
cat("Mortality strategy:", MORT_STRATEGY, "\n")

# Environment variable overrides for ablation study
env_bool <- function(nm, default) {
  v <- Sys.getenv(nm)
  if (v == "") return(default)
  toupper(v) %in% c("TRUE", "1", "YES")
}
OUTPUT_TAG <- Sys.getenv("FVS_OUTPUT_TAG")  # "" = default run, otherwise ablation tag

# Ingrowth toggle
INGROWTH_ENABLED <- env_bool("FVS_INGROWTH", TRUE)
cat("Ingrowth model:", ifelse(INGROWTH_ENABLED, "ENABLED", "DISABLED"), "\n")

# ClimateSI + Emmerson SDIMAX raster lookup
CLIMATE_SI_ENABLED <- env_bool("FVS_CLIMATE_SI", TRUE)
SDIMAX_RASTER_ENABLED <- env_bool("FVS_SDIMAX", TRUE)
RASTER_LOOKUP_CSV <- file.path(project_root, "calibration/data/plot_raster_lookup.csv")
cat("ClimateSI raster:", ifelse(CLIMATE_SI_ENABLED, "ENABLED", "DISABLED"), "\n")
cat("SDIMAX raster:", ifelse(SDIMAX_RASTER_ENABLED, "ENABLED", "DISABLED"), "\n")

# Uncertainty estimation
CI_BRACKETS_ENABLED <- env_bool("FVS_CI_BRACKETS", TRUE)
BOOTSTRAP_CI_ENABLED <- env_bool("FVS_BOOTSTRAP", TRUE)
BOOTSTRAP_N <- 1000
cat("CI bracket projections:", ifelse(CI_BRACKETS_ENABLED, "ENABLED (80% CI)", "DISABLED"), "\n")
cat("Bootstrap CI on stats:", ifelse(BOOTSTRAP_CI_ENABLED, paste0("ENABLED (", BOOTSTRAP_N, " reps)"), "DISABLED"), "\n")
if (OUTPUT_TAG != "") cat("OUTPUT TAG:", OUTPUT_TAG, "\n")

# ==============================================================================
# Helper Functions
# ==============================================================================

calc_ba_vec <- function(dia, tpa) {
  sum(tpa * pi * (dia / 24)^2, na.rm = TRUE)
}

calc_qmd_scalar <- function(ba, tpa) {
  if (is.na(tpa) | is.na(ba) | tpa <= 0 | ba <= 0) return(NA_real_)
  sqrt(ba / tpa * 576 / pi)
}

calc_sdi_scalar <- function(tpa, qmd) {
  if (is.na(tpa) | is.na(qmd) | tpa <= 0 | qmd <= 0) return(NA_real_)
  tpa * (qmd / 10)^1.605
}

safe_scale <- function(x) {
  if (length(x) < 2) return(rep(0, length(x)))
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) | s < 1e-10) return(rep(0, length(x)))
  (x - m) / s
}

# ==============================================================================
# STEP 1: Load Calibrated Parameters (Bayesian posteriors)
# ==============================================================================

cat("STEP 1: Load Calibrated Parameters\n")
cat(strrep("-", 80), "\n\n")

load_variant_params <- function(v) {
  params <- list(variant = v)
  vdir <- file.path(calib_root, v)

  # Quantile extraction helpers (for CI bracket projections)
  CI_LO <- 0.025  # 95% CI on parameters -> ~80% CI on predictions (joint shrinkage)
  CI_HI <- 0.975

  # --- Measurement interval ---
  dg_data_file <- file.path(data_root, v, "diameter_growth.csv")
  if (file.exists(dg_data_file)) {
    dg_data <- fread(dg_data_file, select = "years_interval", nrows = 50000)
    if ("years_interval" %in% names(dg_data)) {
      params$meas_interval <- mean(dg_data$years_interval, na.rm = TRUE)
    }
  }
  if (is.null(params$meas_interval)) params$meas_interval <- 7

  # --- H-D (Chapman-Richards) ---
  hd_file <- file.path(vdir, "height_diameter_samples.rds")
  if (file.exists(hd_file)) {
    draws <- tryCatch(readRDS(hd_file), error = function(e) NULL)
    if (!is.null(draws)) {
      d <- as.data.frame(draws)
      get_med <- function(nm) if (nm %in% names(d)) median(d[[nm]], na.rm = TRUE) else NA_real_
      get_lo  <- function(nm) if (nm %in% names(d)) quantile(d[[nm]], CI_LO, na.rm = TRUE) else NA_real_
      get_hi  <- function(nm) if (nm %in% names(d)) quantile(d[[nm]], CI_HI, na.rm = TRUE) else NA_real_
      params$hd <- list(
        a = get_med("b_a_Intercept"),
        b = get_med("b_b_Intercept"),
        c = get_med("b_c_Intercept")
      )
      params$hd_lo <- list(a = get_lo("b_a_Intercept"), b = get_lo("b_b_Intercept"), c = get_lo("b_c_Intercept"))
      params$hd_hi <- list(a = get_hi("b_a_Intercept"), b = get_hi("b_b_Intercept"), c = get_hi("b_c_Intercept"))
      # Species random effects on a (asymptote) keyed by SPCD
      a_cols <- names(d)[grepl("^r_SPCD__a\\[", names(d))]
      if (length(a_cols) > 0) {
        spcd_a <- as.integer(gsub("r_SPCD__a\\[(\\d+),Intercept\\]", "\\1", a_cols))
        params$hd$sp_a <- setNames(
          sapply(a_cols, function(x) median(d[[x]], na.rm = TRUE)),
          as.character(spcd_a)
        )
        # CI bounds for species RE: use median (RE uncertainty is secondary)
        params$hd_lo$sp_a <- params$hd$sp_a
        params$hd_hi$sp_a <- params$hd$sp_a
      }
    }
  }

  # --- DG (Wykoff ln(DDS)) ---
  dg_file <- file.path(vdir, "diameter_growth_samples.rds")
  std_file <- file.path(vdir, "standardization_params.rds")
  if (file.exists(dg_file)) {
    draws <- tryCatch(readRDS(dg_file), error = function(e) NULL)
    if (!is.null(draws)) {
      d <- as.data.frame(draws)
      get_med <- function(nm) if (nm %in% names(d)) median(d[[nm]], na.rm = TRUE) else 0
      get_lo  <- function(nm) if (nm %in% names(d)) quantile(d[[nm]], CI_LO, na.rm = TRUE) else 0
      get_hi  <- function(nm) if (nm %in% names(d)) quantile(d[[nm]], CI_HI, na.rm = TRUE) else 0
      params$dg <- list(
        mu_b0 = get_med("mu_b0"),
        sigma_b0 = get_med("sigma_b0"),
        betas = sapply(paste0("b", 1:13), get_med),
        sigma = get_med("sigma")
      )
      # CI bracket parameters for DG
      # For CI brackets: use Q2.5/Q97.5 on parameters for systematic uncertainty
      # Also store sigma_shift: a stand-level residual term scaled by 1/sqrt(N_trees)
      # This additive shift on ln(DDS) propagates residual variance through the nonlinear model
      dg_sigma_med <- get_med("sigma")
      params$dg_lo <- list(
        mu_b0 = get_lo("mu_b0"), sigma_b0 = get_med("sigma_b0"),
        betas = sapply(paste0("b", 1:13), get_lo), sigma = dg_sigma_med,
        sigma_shift = -1.28 * dg_sigma_med  # 80% CI lower tail on residuals
      )
      params$dg_hi <- list(
        mu_b0 = get_hi("mu_b0"), sigma_b0 = get_med("sigma_b0"),
        betas = sapply(paste0("b", 1:13), get_hi), sigma = dg_sigma_med,
        sigma_shift = +1.28 * dg_sigma_med  # 80% CI upper tail on residuals
      )
      # Reconstructed species intercepts: b0[i] already in posterior
      b0_cols <- names(d)[grepl("^b0\\[\\d+\\]$", names(d))]
      if (length(b0_cols) > 0) {
        b0_medians <- sapply(b0_cols, function(x) median(d[[x]], na.rm = TRUE))
        b0_lo <- sapply(b0_cols, function(x) quantile(d[[x]], CI_LO, na.rm = TRUE))
        b0_hi <- sapply(b0_cols, function(x) quantile(d[[x]], CI_HI, na.rm = TRUE))
        params$dg$sp_b0_vals <- b0_medians
        params$dg_lo$sp_b0_vals <- b0_lo
        params$dg_hi$sp_b0_vals <- b0_hi
      }

      # Build SPCD -> species_index mapping from training data
      if (file.exists(dg_data_file)) {
        spcd_map_data <- tryCatch({
          fread(dg_data_file, select = c("SPCD", "fvs_sp_idx"), nrows = 500000)
        }, error = function(e) NULL)
        if (!is.null(spcd_map_data) && "fvs_sp_idx" %in% names(spcd_map_data)) {
          sp_map <- unique(spcd_map_data[, .(SPCD, fvs_sp_idx)])
          spcd_idx <- setNames(sp_map$fvs_sp_idx, as.character(sp_map$SPCD))
          params$dg$spcd_to_idx <- spcd_idx
          params$dg_lo$spcd_to_idx <- spcd_idx
          params$dg_hi$spcd_to_idx <- spcd_idx
          cat(sprintf("[%d spp mapped] ", nrow(sp_map)))
        } else {
          spcd_fallback <- tryCatch({
            sdat <- fread(dg_data_file, select = "SPCD", nrows = 500000)
            sort(unique(sdat$SPCD))
          }, error = function(e) NULL)
          if (!is.null(spcd_fallback)) {
            spcd_idx <- setNames(seq_along(spcd_fallback), as.character(spcd_fallback))
            params$dg$spcd_to_idx <- spcd_idx
            params$dg_lo$spcd_to_idx <- spcd_idx
            params$dg_hi$spcd_to_idx <- spcd_idx
            cat(sprintf("[%d spp fallback] ", length(spcd_fallback)))
          }
        }
      }
    }
  }
  if (file.exists(std_file)) {
    params$dg_std <- tryCatch(readRDS(std_file), error = function(e) NULL)
  }

  # --- Mortality (logistic survival) ---
  mort_file <- file.path(vdir, "mortality_samples.rds")
  if (file.exists(mort_file)) {
    draws <- tryCatch(readRDS(mort_file), error = function(e) NULL)
    if (!is.null(draws)) {
      d <- as.data.frame(draws)
      get_fe <- function(nm) {
        col <- paste0("b_", nm)
        if (col %in% names(d)) median(d[[col]], na.rm = TRUE) else 0
      }
      get_fe_lo <- function(nm) {
        col <- paste0("b_", nm)
        if (col %in% names(d)) quantile(d[[col]], CI_LO, na.rm = TRUE) else 0
      }
      get_fe_hi <- function(nm) {
        col <- paste0("b_", nm)
        if (col %in% names(d)) quantile(d[[col]], CI_HI, na.rm = TRUE) else 0
      }
      fe_names <- c("Intercept", "DBH_std", "IDBH_stdE2", "BAL_std", "CR_std", "SI_std", "BA_std")
      params$mort <- setNames(lapply(fe_names, get_fe), paste0("b", 0:6))
      params$mort_lo <- setNames(lapply(fe_names, get_fe_lo), paste0("b", 0:6))
      params$mort_hi <- setNames(lapply(fe_names, get_fe_hi), paste0("b", 0:6))
      # Species random intercepts keyed by SPCD
      sp_cols <- names(d)[grepl("^r_SPCD\\[", names(d))]
      if (length(sp_cols) > 0) {
        spcd_vals <- as.integer(gsub("r_SPCD\\[(\\d+),Intercept\\]", "\\1", sp_cols))
        sp_re_med <- setNames(
          sapply(sp_cols, function(x) median(d[[x]], na.rm = TRUE)),
          as.character(spcd_vals)
        )
        params$mort$sp_re <- sp_re_med
        # Use median RE for brackets (RE uncertainty is secondary to fixed effects)
        params$mort_lo$sp_re <- sp_re_med
        params$mort_hi$sp_re <- sp_re_med
      }
    }
  }

  # --- Standardization params for mortality (from training data) ---
  # IMPORTANT: Training script (04_fit_mortality.R) uses:
  #   DBH_std = scale(DIA_t1)    -- column is DIA in the CSV
  #   CR_std  = scale(CR_pct/100) -- CR_pct/100 = CR on 0-99 scale
  #   SI_std  = scale(log(SI))    -- log-transformed SI
  #   BAL_std = scale(BAL)
  #   BA_std  = scale(BA)
  mort_data_file <- file.path(data_root, v, "mortality.csv")
  if (file.exists(mort_data_file)) {
    mort_data <- tryCatch({
      md <- fread(mort_data_file, nrows = 100000,
                  select = c("DIA", "BAL", "CR", "CR_pct", "SI", "BA"))
      # Use DIA column (training script called it DIA_t1 but CSV stores as DIA)
      dia_vals <- if ("DIA" %in% names(md)) md$DIA else NULL
      # CR: training uses CR_pct/100 which equals CR on 0-99 scale
      cr_vals <- if ("CR" %in% names(md)) md$CR else if ("CR_pct" %in% names(md)) md$CR_pct / 100 else NULL
      # SI: training standardizes log(SI), not raw SI
      si_log_vals <- if ("SI" %in% names(md)) log(pmax(md$SI, 1)) else NULL
      list(
        DBH_mean = mean(dia_vals, na.rm = TRUE), DBH_sd = sd(dia_vals, na.rm = TRUE),
        BAL_mean = mean(md$BAL, na.rm = TRUE),   BAL_sd = sd(md$BAL, na.rm = TRUE),
        CR_mean  = mean(cr_vals, na.rm = TRUE),   CR_sd  = sd(cr_vals, na.rm = TRUE),
        SI_log_mean = mean(si_log_vals, na.rm = TRUE), SI_log_sd = sd(si_log_vals, na.rm = TRUE),
        BA_mean  = mean(md$BA, na.rm = TRUE),     BA_sd  = sd(md$BA, na.rm = TRUE)
      )
    }, error = function(e) NULL)
    params$mort_std <- mort_data
  }

  # --- SDIMAX calibrated values ---
  sdimax_file <- file.path(vdir, "species_sdimax_calibrated.csv")
  if (file.exists(sdimax_file)) {
    sdimax_data <- tryCatch({
      sd <- fread(sdimax_file)
      if ("sdimax_combined" %in% names(sd) & "SPCD" %in% names(sd)) {
        setNames(sd$sdimax_combined, as.character(sd$SPCD))
      } else NULL
    }, error = function(e) NULL)
    params$sdimax <- sdimax_data
  }

  # --- Height Increment (HTG) model ---
  # ln(HTG) = b0[sp] + b1[sp]*ln(DBH) + (b2 + gamma_shape[spgrp])*ln(HT)
  #           + b3*CR + b4*SI + b5*BAL + b6*BA + b7*SLOPE + b8*cos(ASP)
  htg_file <- file.path(vdir, "height_increment_summary.csv")
  htg_method_file <- file.path(vdir, "height_increment_method.txt")
  htg_samples_file <- file.path(vdir, "height_increment_samples.rds")
  if (file.exists(htg_samples_file)) {
    htg_draws <- tryCatch(readRDS(htg_samples_file), error = function(e) NULL)
    if (!is.null(htg_draws)) {
      hd_d <- as.data.frame(htg_draws)
      get_med <- function(nm) if (nm %in% names(hd_d)) median(hd_d[[nm]], na.rm = TRUE) else 0
      get_lo  <- function(nm) if (nm %in% names(hd_d)) quantile(hd_d[[nm]], CI_LO, na.rm = TRUE) else 0
      get_hi  <- function(nm) if (nm %in% names(hd_d)) quantile(hd_d[[nm]], CI_HI, na.rm = TRUE) else 0

      # Fixed effects
      htg_betas <- c(b2 = get_med("b2"), b3 = get_med("b3"), b4 = get_med("b4"),
                     b5 = get_med("b5"), b6 = get_med("b6"), b7 = get_med("b7"),
                     b8 = get_med("b8"))
      htg_sigma <- get_med("sigma")

      # Species intercepts and ln(DBH) coefficients
      b0_cols <- names(hd_d)[grepl("^b0\\[\\d+\\]$", names(hd_d))]
      b1_cols <- names(hd_d)[grepl("^b1\\[\\d+\\]$", names(hd_d))]
      gs_cols <- names(hd_d)[grepl("^gamma_shape\\[\\d+\\]$", names(hd_d))]

      htg_b0 <- sapply(b0_cols, function(x) median(hd_d[[x]], na.rm = TRUE))
      htg_b1 <- sapply(b1_cols, function(x) median(hd_d[[x]], na.rm = TRUE))
      htg_gs <- sapply(gs_cols, function(x) median(hd_d[[x]], na.rm = TRUE))

      params$htg <- list(
        b0 = htg_b0, b1 = htg_b1, betas = htg_betas,
        gamma_shape = htg_gs, sigma = htg_sigma
      )

      # Build SPCD -> species_index mapping from training data
      htg_data_file <- file.path(data_root, v, "height_growth.csv")
      if (file.exists(htg_data_file)) {
        spcd_map_htg <- tryCatch({
          hd_raw <- fread(htg_data_file, select = "SPCD", nrows = 500000)
          sp_list <- sort(unique(hd_raw$SPCD))
          setNames(seq_along(sp_list), as.character(sp_list))
        }, error = function(e) NULL)
        if (!is.null(spcd_map_htg)) {
          params$htg$spcd_to_idx <- spcd_map_htg
          # Species group mapping: cycle through groups
          n_spgrp_htg <- length(htg_gs)
          if (n_spgrp_htg > 0) {
            params$htg$spcd_to_spgrp <- setNames(
              ((seq_along(spcd_map_htg) - 1) %% n_spgrp_htg) + 1,
              names(spcd_map_htg)
            )
          }
        }
      }
      cat("[HTG] ")
    }
  } else if (file.exists(htg_file) && file.info(htg_file)$size > 20) {
    # Try loading from summary CSV (MAP estimates)
    htg_summ <- tryCatch(fread(htg_file), error = function(e) NULL)
    if (!is.null(htg_summ) && nrow(htg_summ) > 0 && "median" %in% names(htg_summ)) {
      get_val <- function(nm) {
        row <- htg_summ[variable == nm]
        if (nrow(row) > 0) row$median[1] else 0
      }
      htg_betas <- c(b2 = get_val("b2"), b3 = get_val("b3"), b4 = get_val("b4"),
                     b5 = get_val("b5"), b6 = get_val("b6"), b7 = get_val("b7"),
                     b8 = get_val("b8"))

      b0_rows <- htg_summ[grepl("^b0\\[", variable)]
      b1_rows <- htg_summ[grepl("^b1\\[", variable)]
      gs_rows <- htg_summ[grepl("^gamma_shape\\[", variable)]

      params$htg <- list(
        b0 = setNames(b0_rows$median, b0_rows$variable),
        b1 = setNames(b1_rows$median, b1_rows$variable),
        betas = htg_betas,
        gamma_shape = setNames(gs_rows$median, gs_rows$variable),
        sigma = get_val("sigma")
      )
      cat("[HTG-MAP] ")
    }
  }

  params
}

# Load all variant parameters
variant_params <- list()
for (v in ALL_VARIANTS) {
  cat("  Loading", toupper(v), "... ")
  p <- load_variant_params(v)
  has_components <- c(
    DG = !is.null(p$dg),
    Mort = !is.null(p$mort),
    HD = !is.null(p$hd),
    HTG = !is.null(p$htg)
  )
  cat(paste(names(has_components)[has_components], collapse = "+"), "\n")
  variant_params[[toupper(v)]] <- p
}
cat("Loaded parameters for", length(variant_params), "variants\n\n")

# ==============================================================================
# STEP 2: Tree-Level Projection Functions
# ==============================================================================

#' Project a condition's tree list forward using calibrated Bayesian parameters
#' @param trees data.table with DIA, TPA_UNADJ, SPCD, CR, HT, VOLCFGRS columns
#' @param params Output of load_variant_params()
#' @param interval_years Years between measurements
#' @param si Site index
#' @return list(TPA_pred, BA_pred, QMD_pred, SDI_pred, VOL_CFGRS_pred)
project_condition_calibrated <- function(trees, params, interval_years, si = 65,
                                        slope_cond = NA, aspect_cond = NA, elev_cond = NA,
                                        sdimax_raster_val = NA) {
  null_result <- list(TPA_pred = NA_real_, BA_pred = NA_real_,
                      QMD_pred = NA_real_, SDI_pred = NA_real_,
                      VOL_CFGRS_pred = NA_real_,
                      VOL_CFNET_pred = NA_real_,
                      VOL_BFNET_pred = NA_real_,
                      VOL_CFGRS_gross = NA_real_,
                      HT_top = NA_real_)
  if (nrow(trees) == 0) return(null_result)

  dt <- copy(trees)
  dt <- dt[!is.na(DIA) & DIA >= 1.0 & !is.na(TPA_UNADJ) & TPA_UNADJ > 0]
  if (nrow(dt) == 0) return(null_result)

  n <- nrow(dt)
  meas_int <- params$meas_interval

  # --- Stand-level predictors ---
  ba_tree <- pi * (dt$DIA / 24)^2
  ba_total <- sum(dt$TPA_UNADJ * ba_tree, na.rm = TRUE)
  tpa_total <- sum(dt$TPA_UNADJ, na.rm = TRUE)

  # BAL: sort descending by diameter, cumulative BA above
  dt <- dt[order(-DIA)]
  ba_tree <- pi * (dt$DIA / 24)^2
  dt[, BAL := cumsum(TPA_UNADJ * ba_tree) - TPA_UNADJ * ba_tree]

  # Crown ratio: use FIA value if available (FIA stores 0-99), else 0.40
  dt[, CR_prop := fifelse(!is.na(CR) & CR > 0 & CR <= 99, CR / 100, 0.40)]

  # ---- MORTALITY: expected TPA reduction ----
  if (!is.null(params$mort)) {
    m <- params$mort

    # Species random intercepts
    sp_re <- rep(0, n)
    if (!is.null(m$sp_re)) {
      spcd_char <- as.character(dt$SPCD)
      matched <- spcd_char %in% names(m$sp_re)
      if (any(matched)) {
        sp_re[matched] <- m$sp_re[spcd_char[matched]]
      }
    }

    # Compute logit of survival based on strategy
    if (MORT_STRATEGY == "re_only") {
      # Species RE + intercept only: no covariate fixed effects
      logit_surv_period <- m$b0 + sp_re

    } else if (MORT_STRATEGY == "within") {
      # Within-condition standardization: captures relative competitive position
      # SI and BA are constant within condition => zeroed after centering
      dbh_s <- safe_scale(dt$DIA)
      bal_s <- safe_scale(dt$BAL)
      cr_s  <- safe_scale(dt$CR_prop * 100)  # 0-99 scale
      si_s  <- rep(0, n)
      ba_s  <- rep(0, n)

      logit_surv_period <- (m$b0 + sp_re) +
        m$b1 * dbh_s + m$b2 * dbh_s^2 +
        m$b3 * bal_s + m$b4 * cr_s +
        m$b5 * si_s + m$b6 * ba_s

    } else if (MORT_STRATEGY == "no_dbh_si") {
      # Training-data standardization for BAL/CR/BA; DBH and SI zeroed
      # Removes size-dependent and site effects from mortality
      if (!is.null(params$mort_std)) {
        ms <- params$mort_std
        std_safe <- function(x, mn, sd_val) {
          if (is.null(mn) | is.null(sd_val) | is.na(sd_val) | sd_val < 1e-10) return(rep(0, length(x)))
          (x - mn) / sd_val
        }
        bal_s <- std_safe(dt$BAL, ms$BAL_mean, ms$BAL_sd)
        cr_s  <- std_safe(dt$CR_prop * 100, ms$CR_mean, ms$CR_sd)
        ba_s  <- std_safe(rep(ba_total, n), ms$BA_mean, ms$BA_sd)
      } else {
        bal_s <- safe_scale(dt$BAL)
        cr_s  <- safe_scale(dt$CR_prop * 100)
        ba_s  <- rep(0, n)
      }
      dbh_s <- rep(0, n)
      si_s  <- rep(0, n)

      logit_surv_period <- (m$b0 + sp_re) +
        m$b1 * dbh_s + m$b2 * dbh_s^2 +
        m$b3 * bal_s + m$b4 * cr_s +
        m$b5 * si_s + m$b6 * ba_s

    } else if (MORT_STRATEGY == "no_si") {
      # Training-data standardization for DBH/BAL/CR/BA; only SI zeroed
      # Retains size-dependent mortality (U-shaped via DBH + DBH^2)
      # This captures senescence mortality in large trees
      if (!is.null(params$mort_std)) {
        ms <- params$mort_std
        std_safe <- function(x, mn, sd_val) {
          if (is.null(mn) | is.null(sd_val) | is.na(sd_val) | sd_val < 1e-10) return(rep(0, length(x)))
          (x - mn) / sd_val
        }
        dbh_s <- std_safe(dt$DIA, ms$DBH_mean, ms$DBH_sd)
        bal_s <- std_safe(dt$BAL, ms$BAL_mean, ms$BAL_sd)
        cr_s  <- std_safe(dt$CR_prop * 100, ms$CR_mean, ms$CR_sd)
        ba_s  <- std_safe(rep(ba_total, n), ms$BA_mean, ms$BA_sd)
      } else {
        dbh_s <- safe_scale(dt$DIA)
        bal_s <- safe_scale(dt$BAL)
        cr_s  <- safe_scale(dt$CR_prop * 100)
        ba_s  <- rep(0, n)
      }
      si_s  <- rep(0, n)

      logit_surv_period <- (m$b0 + sp_re) +
        m$b1 * dbh_s + m$b2 * dbh_s^2 +
        m$b3 * bal_s + m$b4 * cr_s +
        m$b5 * si_s + m$b6 * ba_s

    } else {
      # "full": training-data standardization with all covariates
      if (!is.null(params$mort_std)) {
        ms <- params$mort_std
        std_safe <- function(x, mn, sd_val) {
          if (is.null(mn) | is.null(sd_val) | is.na(sd_val) | sd_val < 1e-10) return(rep(0, length(x)))
          (x - mn) / sd_val
        }
        dbh_s   <- std_safe(dt$DIA, ms$DBH_mean, ms$DBH_sd)
        bal_s   <- std_safe(dt$BAL, ms$BAL_mean, ms$BAL_sd)
        cr_s    <- std_safe(dt$CR_prop * 100, ms$CR_mean, ms$CR_sd)
        si_log  <- log(pmax(si, 1))
        si_s    <- std_safe(rep(si_log, n), ms$SI_log_mean, ms$SI_log_sd)
        ba_s    <- std_safe(rep(ba_total, n), ms$BA_mean, ms$BA_sd)
      } else {
        # Fallback to within-condition if training stats unavailable
        dbh_s <- safe_scale(dt$DIA)
        bal_s <- safe_scale(dt$BAL)
        cr_s  <- safe_scale(dt$CR_prop * 100)
        si_s  <- rep(0, n)
        ba_s  <- rep(0, n)
      }

      logit_surv_period <- (m$b0 + sp_re) +
        m$b1 * dbh_s + m$b2 * dbh_s^2 +
        m$b3 * bal_s + m$b4 * cr_s +
        m$b5 * si_s + m$b6 * ba_s
    }

    logit_surv_period[!is.finite(logit_surv_period)] <- 2

    # ---- SDIMAX RELATIVE DENSITY ADJUSTMENT ----
    # If Emmerson SDIMAX is available, apply a carrying capacity modifier:
    # When SDI/SDIMAX > 0.55 (Reineke's zone of imminent mortality), increase mort
    # When SDI/SDIMAX < 0.35 (understocked), decrease mort
    # This is a post-hoc density dependent modifier, not a regression covariate
    if (SDIMAX_RASTER_ENABLED && !is.null(params$sdimax)) {
      # Compute current SDI: Reineke's SDI = sum(TPA * (DIA/10)^1.605)
      sdi_current <- sum(dt$TPA_UNADJ * (dt$DIA / 10)^1.605, na.rm = TRUE)

      # Get species-weighted SDIMAX from Emmerson raster or calibrated values
      sdimax_val <- NA_real_
      # Prefer plot-level Emmerson raster value if available
      if (!is.null(sdimax_raster_val) && !is.na(sdimax_raster_val) && sdimax_raster_val > 0) {
        sdimax_val <- sdimax_raster_val
      } else {
        # Fall back to species-specific calibrated SDIMAX (TPA-weighted mean)
        spcd_char <- as.character(dt$SPCD)
        sp_sdi <- params$sdimax[spcd_char]
        sp_sdi <- sp_sdi[!is.na(sp_sdi) & sp_sdi > 0]
        if (length(sp_sdi) > 0) sdimax_val <- mean(sp_sdi)
      }

      if (!is.na(sdimax_val) && sdimax_val > 50) {
        rdi <- sdi_current / sdimax_val
        rdi <- pmin(rdi, 1.5)  # cap at 150% to avoid extreme adjustments

        # Logistic modifier centered at RDI=0.55 (onset of competition mortality)
        # At RDI=0.55: modifier=0 (neutral)
        # At RDI=0.80: modifier increases mortality
        # At RDI=0.30: modifier decreases mortality
        rdi_modifier <- -2.0 * (rdi - 0.55)  # subtract from logit: higher RDI -> lower survival
        logit_surv_period <- logit_surv_period + rdi_modifier
      }
    }

    p_surv_period <- plogis(logit_surv_period)
    p_surv_period <- pmin(pmax(p_surv_period, 0.01), 0.999)

    # Scale to actual interval:
    # Model was trained on data with mean interval = meas_int
    # For interval_years: p_surv_actual = p_surv_period^(interval_years / meas_int)
    p_surv_actual <- p_surv_period^(interval_years / meas_int)

    p_surv_actual <- pmin(pmax(p_surv_actual, 0.01), 0.999)
    dt[, TPA_pred := TPA_UNADJ * p_surv_actual]
  } else {
    # No calibrated mortality: use 1.5% annual default
    dt[, TPA_pred := TPA_UNADJ * (1 - 0.015)^interval_years]
  }

  # ---- DIAMETER GROWTH: Wykoff ln(DDS) ----
  if (!is.null(params$dg) && !is.null(params$dg_std)) {
    dg  <- params$dg
    std <- params$dg_std

    safe_std <- function(x, nm) {
      mn <- std[[paste0(nm, "_mean")]]
      sd_val <- std[[paste0(nm, "_sd")]]
      if (!is.null(mn) && !is.null(sd_val) && !is.na(sd_val) && sd_val > 1e-10) {
        return((x - mn) / sd_val)
      }
      rep(0, length(x))
    }

    ln_dbh <- log(pmax(dt$DIA, 0.1))
    dbh_sq <- dt$DIA^2
    ln_si  <- log(pmax(si, 1))

    # Use actual condition-level site/topo from FIA, fallback to means
    slope_val <- fifelse(!is.na(slope_cond) & slope_cond >= 0, slope_cond,
                         ifelse(!is.null(std$SLOPE_mean), std$SLOPE_mean, 15))
    elev_val  <- fifelse(!is.na(elev_cond) & elev_cond > 0, elev_cond,
                         ifelse(!is.null(std$ELEV_mean), std$ELEV_mean, 1500))
    asp_rad   <- fifelse(!is.na(aspect_cond) & aspect_cond >= 0,
                         aspect_cond * pi / 180, 0)
    slope_sasp <- slope_val * sin(asp_rad)
    slope_casp <- slope_val * cos(asp_rad)

    ln_dbh_s    <- safe_std(ln_dbh, "ln_DBH")
    dbh_sq_s    <- safe_std(dbh_sq, "DBH_sq")
    ln_si_s     <- safe_std(rep(ln_si, n), "ln_SI")
    slope_s     <- safe_std(rep(slope_val, n), "SLOPE")
    slope_sasp_s <- safe_std(rep(slope_sasp, n), "SLOPE_SASP")
    slope_casp_s <- safe_std(rep(slope_casp, n), "SLOPE_CASP")
    elev_s      <- safe_std(rep(elev_val, n), "ELEV")
    bal_s_dg    <- safe_std(dt$BAL, "BAL")
    ba_s_dg     <- safe_std(rep(ba_total, n), "BA")
    cr_dg       <- dt$CR_prop   # CR on 0-1 scale (not standardized per Stan model)

    # Guard NaN
    for (vec_name in c("ln_dbh_s", "dbh_sq_s", "ln_si_s", "slope_s",
                       "slope_sasp_s", "slope_casp_s", "elev_s", "bal_s_dg", "ba_s_dg")) {
      vv <- get(vec_name)
      vv[!is.finite(vv)] <- 0
      assign(vec_name, vv)
    }
    cr_dg[!is.finite(cr_dg)] <- 0.40

    b <- dg$betas
    b[!is.finite(b)] <- 0

    # Species-specific intercepts: look up b0[species_index] for each tree
    sp_int <- rep(dg$mu_b0, n)  # default to grand mean
    if (!is.null(dg$sp_b0_vals) && !is.null(dg$spcd_to_idx)) {
      spcd_char <- as.character(dt$SPCD)
      matched_sp <- spcd_char %in% names(dg$spcd_to_idx)
      if (any(matched_sp)) {
        sp_indices <- dg$spcd_to_idx[spcd_char[matched_sp]]
        # b0 columns are named "b0[1]", "b0[2]", etc.
        valid_idx <- sp_indices >= 1 & sp_indices <= length(dg$sp_b0_vals)
        sp_int[matched_sp][valid_idx] <- dg$sp_b0_vals[sp_indices[valid_idx]]
      }
    }

    # Stan model: b0[sp] + b1*ln_DBH + b2*DBH_sq + b3*ln_SI +
    #             b4*SLOPE + b5*SLOPE_sq + b6*SLOPE_SASP + b7*SLOPE_CASP +
    #             b8*ELEV + b9*ELEV_sq + b10*CR + b11*CR_sq + b12*BAL + b13*BA
    pred_ln_dds <- sp_int +
      b[1]*ln_dbh_s + b[2]*dbh_sq_s + b[3]*ln_si_s +
      b[4]*slope_s + b[5]*slope_s^2 + b[6]*slope_sasp_s + b[7]*slope_casp_s +
      b[8]*elev_s + b[9]*elev_s^2 +
      b[10]*cr_dg + b[11]*cr_dg^2 + b[12]*bal_s_dg + b[13]*ba_s_dg

    pred_ln_dds[!is.finite(pred_ln_dds)] <- log(0.5)

    # Apply residual sigma shift for CI bracket projections
    # sigma_shift is ±1.28*sigma (80% CI), scaled by 1/sqrt(N) at tree level
    # This captures within-model prediction uncertainty that parameter swaps miss
    if (!is.null(dg$sigma_shift) && is.finite(dg$sigma_shift)) {
      n_eff <- max(sqrt(n), 1)  # effective sample size for stand-level averaging
      pred_ln_dds <- pred_ln_dds + dg$sigma_shift / n_eff
    }

    # Back-transform with bias correction
    sigma_sq <- ifelse(is.finite(dg$sigma), dg$sigma^2, 0)
    dds_period <- exp(pred_ln_dds + sigma_sq / 2)

    # Scale DDS to actual interval
    dds_actual <- dds_period * (interval_years / meas_int)
    dds_actual <- pmax(dds_actual, 0)
    dds_actual[!is.finite(dds_actual)] <- 0

    # D2 = sqrt(D1^2 + DDS)
    dt[, DIA_pred := pmax(sqrt(DIA^2 + dds_actual), DIA)]
  } else {
    # Default: 0.10 in/yr
    dt[, DIA_pred := DIA + 0.10 * interval_years]
  }

  # ---- HEIGHT PREDICTION ----
  # Strategy (following FVS workflow):
  #   1. Impute missing heights at t1 using calibrated H-D equation
  #   2. Predict height increment using calibrated HTG model (if available)
  #      OR use H-D ratio approach as fallback
  #   3. HT_pred = HT_t1 + HTG (increment model) or HT_t1 * ratio (H-D fallback)
  # The H-D model is always needed for imputation and volume scaling

  has_obs_ht <- !is.na(dt$HT) & dt$HT > 4.5

  # Step 1: Compute H-D predictions at t1 and t2 diameters (needed for imputation + volume)
  if (!is.null(params$hd)) {
    hd <- params$hd
    a_base <- hd$a; b_val <- hd$b; c_val <- hd$c
    a_vec <- rep(a_base, n)
    if (!is.null(hd$sp_a)) {
      spcd_char <- as.character(dt$SPCD)
      sp_match <- spcd_char %in% names(hd$sp_a)
      if (any(sp_match)) a_vec[sp_match] <- a_base + hd$sp_a[spcd_char[sp_match]]
    }
    ht_pred_t1 <- pmax(4.5 + a_vec * (1 - exp(-b_val * pmax(dt$DIA, 0.1)))^c_val, 5)
    ht_pred_t2 <- pmax(4.5 + a_vec * (1 - exp(-b_val * pmax(dt$DIA_pred, 0.1)))^c_val, 5)
  } else {
    # No H-D model: use observed HT or NA
    ht_pred_t1 <- dt$HT
    ht_pred_t2 <- fifelse(has_obs_ht, dt$HT * (dt$DIA_pred / pmax(dt$DIA, 0.1))^0.5, NA_real_)
  }

  # Step 1b: Establish initial height (impute missing with calibrated H-D)
  ht_t1 <- fifelse(has_obs_ht, dt$HT, ht_pred_t1)

  # Step 2: Predict height at t2
  # Primary method: Calibrated H-D ratio anchored to observed/imputed height
  # This approach leverages the calibrated H-D equation while preserving the
  # observed height anchor. The H-D model captures species-specific allometry,
  # and the ratio HD(DIA_pred)/HD(DIA_t1) transfers that growth signal.
  #
  # The standalone HTG increment model (b0+b1*ln(DBH)+(b2+gamma)*ln(HT)+...)
  # is available but currently has unconstrained gamma_shape parameters that
  # make the effective ln(HT) coefficient positive for many variants, causing
  # biologically unreasonable predictions (taller trees grow faster -> runaway).
  # HTG can be enabled for testing via FVS_USE_HTG=true environment variable.
  use_htg <- tolower(Sys.getenv("FVS_USE_HTG", "false")) == "true"

  if (use_htg && !is.null(params$htg)) {
    # --- Experimental: Calibrated Height Increment Model ---
    htg <- params$htg
    betas <- htg$betas

    spcd_char <- as.character(dt$SPCD)
    htg_intercept <- rep(0, n)
    htg_b1_val    <- rep(0, n)
    htg_gamma_val <- rep(0, n)

    if (!is.null(htg$spcd_to_idx)) {
      sp_matched <- spcd_char %in% names(htg$spcd_to_idx)
      if (any(sp_matched)) {
        sp_idx <- htg$spcd_to_idx[spcd_char[sp_matched]]
        b0_keys <- paste0("b0[", sp_idx, "]")
        b1_keys <- paste0("b1[", sp_idx, "]")
        valid_b0 <- b0_keys %in% names(htg$b0)
        valid_b1 <- b1_keys %in% names(htg$b1)
        if (any(valid_b0)) htg_intercept[sp_matched][valid_b0] <- htg$b0[b0_keys[valid_b0]]
        if (any(valid_b1)) htg_b1_val[sp_matched][valid_b1] <- htg$b1[b1_keys[valid_b1]]
      }
      if (!is.null(htg$spcd_to_spgrp)) {
        spgrp_matched <- spcd_char %in% names(htg$spcd_to_spgrp)
        if (any(spgrp_matched)) {
          grp_idx <- htg$spcd_to_spgrp[spcd_char[spgrp_matched]]
          gs_keys <- paste0("gamma_shape[", grp_idx, "]")
          valid_gs <- gs_keys %in% names(htg$gamma_shape)
          if (any(valid_gs)) htg_gamma_val[spgrp_matched][valid_gs] <- htg$gamma_shape[gs_keys[valid_gs]]
        }
      }
    }

    # Constrain effective ln(HT) coefficient to be non-positive
    # Positive values create runaway: taller trees -> more height growth -> taller
    eff_lnht_coef <- pmin(betas["b2"] + htg_gamma_val, 0)

    ln_dbh_htg <- log(pmax(dt$DIA, 0.1))
    ln_ht_htg  <- log(pmax(ht_t1, 5))
    cr_htg     <- pmin(pmax(fifelse(!is.na(dt$CR_prop), dt$CR_prop, 0.40), 0.01), 1.0)
    si_htg     <- si / 100
    bal_htg    <- dt$BAL / 100
    ba_htg     <- ba_total / 100
    slope_htg  <- fifelse(!is.na(slope_cond), slope_cond, 15) / 100
    casp_htg   <- fifelse(!is.na(aspect_cond) & aspect_cond >= 0,
                           cos(aspect_cond * pi / 180), 0)

    pred_ln_htg <- htg_intercept +
      htg_b1_val * ln_dbh_htg +
      eff_lnht_coef * ln_ht_htg +
      betas["b3"] * cr_htg +
      betas["b4"] * si_htg +
      betas["b5"] * bal_htg +
      betas["b6"] * ba_htg +
      betas["b7"] * slope_htg +
      betas["b8"] * casp_htg

    pred_ln_htg[!is.finite(pred_ln_htg)] <- log(0.5)
    htg_annual <- exp(pred_ln_htg)
    htg_annual <- pmin(pmax(htg_annual, 0.01), 5)
    htg_total  <- htg_annual * interval_years
    dt[, HT_pred := pmax(ht_t1 + htg_total, 5)]

  } else if (!is.null(params$hd)) {
    # --- Primary: Calibrated H-D ratio anchored to initial height ---
    # HT_pred = HT_t1 * [HD(DIA_pred) / HD(DIA_t1)]
    # This preserves the observed height anchor while using the calibrated
    # H-D equation to predict the proportional change in height.
    ht_ratio_hd <- ht_pred_t2 / pmax(ht_pred_t1, 5)
    ht_ratio_hd <- pmin(pmax(ht_ratio_hd, 0.95), 2.0)
    dt[, HT_pred := pmax(ht_t1 * ht_ratio_hd, 5)]

  } else {
    # --- Last resort: simple allometric scaling ---
    dt[, HT_pred := fifelse(has_obs_ht, HT * (DIA_pred / pmax(DIA, 0.1))^0.5, NA_real_)]
    ht_pred_t1 <- dt$HT
    ht_pred_t2 <- dt$HT_pred
  }

  # ---- VOLUME PREDICTION: combined variable ratio ----
  # VOL_pred_i = VOL_t1_i * (DIA_pred^2 * HT_pred_t2) / (DIA^2 * HT_t1)
  has_ht <- !is.na(ht_pred_t1) & ht_pred_t1 > 0 & !is.na(ht_pred_t2) & dt$DIA > 0
  cv_ratio_vec <- rep(1, n)
  if (any(has_ht)) {
    cv_ratio_vec[has_ht] <- (dt$DIA_pred[has_ht]^2 * ht_pred_t2[has_ht]) /
                             (dt$DIA[has_ht]^2 * ht_pred_t1[has_ht])
    cv_ratio_vec <- pmin(pmax(cv_ratio_vec, 0.5), 5.0)
  }

  # Gross cubic foot
  has_cfgrs <- !is.na(dt$VOLCFGRS) & dt$VOLCFGRS > 0 & has_ht
  dt[, VOL_pred := 0]
  if (any(has_cfgrs)) dt$VOL_pred[has_cfgrs] <- dt$VOLCFGRS[has_cfgrs] * cv_ratio_vec[has_cfgrs]

  # Net cubic foot
  has_cfnet <- !is.na(dt$VOLCFNET) & dt$VOLCFNET > 0 & has_ht
  dt[, VOL_CFNET_pred := 0]
  if (any(has_cfnet)) dt$VOL_CFNET_pred[has_cfnet] <- dt$VOLCFNET[has_cfnet] * cv_ratio_vec[has_cfnet]

  # Net board foot
  has_bfnet <- !is.na(dt$VOLBFNET) & dt$VOLBFNET > 0 & has_ht
  dt[, VOL_BFNET_pred := 0]
  if (any(has_bfnet)) dt$VOL_BFNET_pred[has_bfnet] <- dt$VOLBFNET[has_bfnet] * cv_ratio_vec[has_bfnet]

  # ---- AGGREGATE to stand-level predicted metrics ----
  ba_pred  <- sum(dt$TPA_pred * pi * (dt$DIA_pred / 24)^2, na.rm = TRUE)
  tpa_pred <- sum(dt$TPA_pred, na.rm = TRUE)
  vol_cfgrs_pred <- sum(dt$TPA_pred * dt$VOL_pred, na.rm = TRUE)
  vol_cfnet_pred <- sum(dt$TPA_pred * dt$VOL_CFNET_pred, na.rm = TRUE)
  vol_bfnet_pred <- sum(dt$TPA_pred * dt$VOL_BFNET_pred, na.rm = TRUE)

  # ---- GROSS VOLUME: what survivors would produce if no mortality ----
  # Gross = growth on all trees (using original TPA_UNADJ), ignoring mortality
  vol_cfgrs_gross <- sum(dt$TPA_UNADJ * dt$VOL_pred, na.rm = TRUE)

  # ---- TOP HEIGHT: TPA-weighted top 40 trees/ac (Lord method) ----
  ht_top <- NA_real_
  if (!is.null(dt$HT_pred) && any(!is.na(dt$HT_pred))) {
    dt_ht <- dt[!is.na(HT_pred) & HT_pred > 0][order(-HT_pred)]
    if (nrow(dt_ht) > 0) {
      cum_tpa <- cumsum(dt_ht$TPA_pred)
      top_idx <- cum_tpa <= 40  # top 40 TPA (≈ 100 TPH)
      if (sum(top_idx) == 0) top_idx[1] <- TRUE
      ht_top <- weighted.mean(dt_ht$HT_pred[top_idx], dt_ht$TPA_pred[top_idx])
    }
  }

  # ---- INGROWTH: add empirical ingrowth scaled to interval ----
  if (INGROWTH_ENABLED && !is.null(ingrowth_lookup)) {
    var_upper <- toupper(params$variant)
    ig <- ingrowth_lookup[[var_upper]]
    if (is.null(ig)) ig <- ingrowth_lookup[["OVERALL"]]
    if (!is.null(ig)) {
      tpa_pred <- tpa_pred + ig$med_ann_TPA * interval_years
      ba_pred  <- ba_pred  + ig$med_ann_BA  * interval_years
      vol_cfgrs_pred <- vol_cfgrs_pred + ig$med_ann_VOL_CFGRS * interval_years
      vol_cfnet_pred <- vol_cfnet_pred + ig$med_ann_VOL_CFNET * interval_years
      vol_bfnet_pred <- vol_bfnet_pred + ig$med_ann_VOL_BFNET * interval_years
    }
  }

  qmd_pred <- calc_qmd_scalar(ba_pred, tpa_pred)
  sdi_pred <- calc_sdi_scalar(tpa_pred, qmd_pred)

  list(TPA_pred = tpa_pred, BA_pred = ba_pred,
       QMD_pred = qmd_pred, SDI_pred = sdi_pred,
       VOL_CFGRS_pred = vol_cfgrs_pred,
       VOL_CFNET_pred = vol_cfnet_pred,
       VOL_BFNET_pred = vol_bfnet_pred,
       VOL_CFGRS_gross = vol_cfgrs_gross,
       HT_top = ht_top)
}


#' Project a condition's tree list forward using default FVS-like parameters
project_condition_default <- function(trees, interval_years, variant_code = NULL) {
  null_result <- list(TPA_pred = NA_real_, BA_pred = NA_real_,
                      QMD_pred = NA_real_, SDI_pred = NA_real_,
                      VOL_CFGRS_pred = NA_real_,
                      VOL_CFNET_pred = NA_real_,
                      VOL_BFNET_pred = NA_real_,
                      VOL_CFGRS_gross = NA_real_,
                      HT_top = NA_real_)
  if (nrow(trees) == 0) return(null_result)

  dt <- copy(trees)
  dt <- dt[!is.na(DIA) & DIA >= 1.0 & !is.na(TPA_UNADJ) & TPA_UNADJ > 0]
  if (nrow(dt) == 0) return(null_result)

  # Stand metrics
  ba_tree <- pi * (dt$DIA / 24)^2
  ba_total <- sum(dt$TPA_UNADJ * ba_tree, na.rm = TRUE)

  # BAL
  dt <- dt[order(-DIA)]
  ba_tree <- pi * (dt$DIA / 24)^2
  dt[, BAL := cumsum(TPA_UNADJ * ba_tree) - TPA_UNADJ * ba_tree]

  # Default mortality: 1.5% annual, modified by size and competition
  base_mort <- 0.015
  size_mod <- fifelse(dt$DIA < 5, 1.5, fifelse(dt$DIA > 20, 1.2, 1.0))
  comp_mod <- 1 + 0.5 * pmin(ba_total / 200, 1)
  annual_surv <- pmax(1 - base_mort * size_mod * comp_mod, 0.50)
  p_surv_period <- annual_surv^interval_years
  dt[, TPA_pred := TPA_UNADJ * p_surv_period]

  # Default DG: 0.10 in/yr with competition and size modifiers
  competition_mod <- pmax(0.3, 1 - dt$BAL / pmax(ba_total + 1, 10))
  dg_size_mod <- pmax(0.5, 1 - (dt$DIA / 30)^0.5)
  dt[, DIA_pred := DIA + 0.10 * interval_years * competition_mod * dg_size_mod]

  # Default volume: scale by BA ratio (DIA_pred^2 / DIA^2)
  ba_ratio_vec <- rep(1, nrow(dt))
  has_dia <- dt$DIA > 0
  if (any(has_dia)) {
    ba_ratio_vec[has_dia] <- (dt$DIA_pred[has_dia] / dt$DIA[has_dia])^2
    ba_ratio_vec <- pmin(pmax(ba_ratio_vec, 0.5), 5.0)
  }

  has_cfgrs <- !is.na(dt$VOLCFGRS) & dt$VOLCFGRS > 0 & has_dia
  dt[, VOL_pred := 0]
  if (any(has_cfgrs)) dt$VOL_pred[has_cfgrs] <- dt$VOLCFGRS[has_cfgrs] * ba_ratio_vec[has_cfgrs]

  has_cfnet <- !is.na(dt$VOLCFNET) & dt$VOLCFNET > 0 & has_dia
  dt[, VOL_CFNET_pred := 0]
  if (any(has_cfnet)) dt$VOL_CFNET_pred[has_cfnet] <- dt$VOLCFNET[has_cfnet] * ba_ratio_vec[has_cfnet]

  has_bfnet <- !is.na(dt$VOLBFNET) & dt$VOLBFNET > 0 & has_dia
  dt[, VOL_BFNET_pred := 0]
  if (any(has_bfnet)) dt$VOL_BFNET_pred[has_bfnet] <- dt$VOLBFNET[has_bfnet] * ba_ratio_vec[has_bfnet]

  # Aggregate
  ba_pred  <- sum(dt$TPA_pred * pi * (dt$DIA_pred / 24)^2, na.rm = TRUE)
  tpa_pred <- sum(dt$TPA_pred, na.rm = TRUE)
  vol_cfgrs_pred <- sum(dt$TPA_pred * dt$VOL_pred, na.rm = TRUE)
  vol_cfnet_pred <- sum(dt$TPA_pred * dt$VOL_CFNET_pred, na.rm = TRUE)
  vol_bfnet_pred <- sum(dt$TPA_pred * dt$VOL_BFNET_pred, na.rm = TRUE)

  # Gross volume: growth on all trees as if no mortality
  vol_cfgrs_gross <- sum(dt$TPA_UNADJ * dt$VOL_pred, na.rm = TRUE)

  # Top height: TPA-weighted top 40 trees/ac
  ht_top <- NA_real_
  ht_avail <- !is.na(dt$HT) & dt$HT > 0
  if (any(ht_avail)) {
    ht_ratio <- fifelse(dt$DIA > 0, (dt$DIA_pred / pmax(dt$DIA, 0.1))^0.5, 1)
    dt[, HT_pred_def := fifelse(ht_avail, HT * ht_ratio, NA_real_)]
    dt_ht <- dt[!is.na(HT_pred_def) & HT_pred_def > 0][order(-HT_pred_def)]
    if (nrow(dt_ht) > 0) {
      cum_tpa <- cumsum(dt_ht$TPA_pred)
      top_idx <- cum_tpa <= 40
      if (sum(top_idx) == 0) top_idx[1] <- TRUE
      ht_top <- weighted.mean(dt_ht$HT_pred_def[top_idx], dt_ht$TPA_pred[top_idx])
    }
  }

  # ---- INGROWTH: add empirical ingrowth scaled to interval ----
  # Default projection also gets ingrowth for fair comparison
  if (INGROWTH_ENABLED && !is.null(ingrowth_lookup)) {
    ig <- if (!is.null(variant_code)) ingrowth_lookup[[variant_code]] else NULL
    if (is.null(ig)) ig <- ingrowth_lookup[["OVERALL"]]
    if (!is.null(ig)) {
      tpa_pred <- tpa_pred + ig$med_ann_TPA * interval_years
      ba_pred  <- ba_pred  + ig$med_ann_BA  * interval_years
      vol_cfgrs_pred <- vol_cfgrs_pred + ig$med_ann_VOL_CFGRS * interval_years
      vol_cfnet_pred <- vol_cfnet_pred + ig$med_ann_VOL_CFNET * interval_years
      vol_bfnet_pred <- vol_bfnet_pred + ig$med_ann_VOL_BFNET * interval_years
    }
  }

  qmd_pred <- calc_qmd_scalar(ba_pred, tpa_pred)
  sdi_pred <- calc_sdi_scalar(tpa_pred, qmd_pred)

  list(TPA_pred = tpa_pred, BA_pred = ba_pred,
       QMD_pred = qmd_pred, SDI_pred = sdi_pred,
       VOL_CFGRS_pred = vol_cfgrs_pred,
       VOL_CFNET_pred = vol_cfnet_pred,
       VOL_BFNET_pred = vol_bfnet_pred,
       VOL_CFGRS_gross = vol_cfgrs_gross,
       HT_top = ht_top)
}


# ==============================================================================
# STEP 2: Data Assembly
# ==============================================================================

cat("STEP 2: Data Assembly\n")
cat(strrep("-", 80), "\n\n")

# Read plot data (remeasured plots only)
cat("Reading ENTIRE_PLOT.csv ...\n")
plots <- fread(
  file.path(fia_root, "ENTIRE_PLOT.csv"),
  select = c("CN", "PREV_PLT_CN", "INVYR", "STATECD", "LAT", "LON", "ELEV"),
  colClasses = c(CN = "character", PREV_PLT_CN = "character", INVYR = "integer")
)
plots <- plots[!is.na(PREV_PLT_CN) & PREV_PLT_CN != ""]
cat("  Remeasured plots:", nrow(plots), "\n")

# Build plot pairs: time2 CN -> time1 CN
plot_pairs <- unique(plots[, .(PLT_CN_t2 = CN, PLT_CN_t1 = PREV_PLT_CN,
                                INVYR_t2 = INVYR, ELEV_t2 = ELEV)])

# We need INVYR for time1 plots. Look it up from plots table (CN -> INVYR)
all_plots_invyr <- fread(
  file.path(fia_root, "ENTIRE_PLOT.csv"),
  select = c("CN", "INVYR"),
  colClasses = c(CN = "character", INVYR = "integer")
)
setnames(all_plots_invyr, c("PLT_CN_t1", "INVYR_t1"))
plot_pairs <- merge(plot_pairs, all_plots_invyr, by = "PLT_CN_t1", all.x = TRUE)
plot_pairs[, interval_years := INVYR_t2 - INVYR_t1]
plot_pairs <- plot_pairs[!is.na(interval_years) & interval_years > 0 & interval_years <= 20]
cat("  Valid plot pairs:", nrow(plot_pairs), "\n")
cat("  Mean interval:", round(mean(plot_pairs$interval_years), 2), "years\n")

# FVS variant mapping
cat("Reading FVS variant mapping...\n")
variant_map <- fread(
  file.path(fia_root, "ENTIRE_FVS_STANDINIT_COND.csv"),
  select = c("STAND_CN", "VARIANT"),
  colClasses = c(STAND_CN = "character", VARIANT = "character")
)
variant_map <- unique(variant_map)

# FVS STANDINIT_PLOT: pull FVS-assigned SITE_INDEX and MAX_SDI for enrichment
cat("Reading FVS STANDINIT_PLOT for SITE_INDEX and MAX_SDI...\n")
fvs_plot_init_file <- file.path(fia_root, "ENTIRE_FVS_STANDINIT_PLOT.csv")
if (file.exists(fvs_plot_init_file)) {
  fvs_plot_init <- tryCatch({
    fread(fvs_plot_init_file,
          select = c("STAND_CN", "SITE_INDEX", "MAX_SDI", "MAX_SDI_FIA",
                      "SITE_SPECIES", "ELEVATION", "FOREST_TYPE"),
          colClasses = c(STAND_CN = "character"))
  }, error = function(e) {
    cat("  WARNING: Could not read FVS STANDINIT_PLOT:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(fvs_plot_init)) {
    fvs_plot_init <- unique(fvs_plot_init)
    cat("  FVS plot init records:", nrow(fvs_plot_init), "\n")
  }
} else {
  cat("  FVS STANDINIT_PLOT not found, skipping\n")
  fvs_plot_init <- NULL
}

# Collect all plot CNs involved
all_plt_cns <- unique(c(plot_pairs$PLT_CN_t1, plot_pairs$PLT_CN_t2))
cat("  Total unique plot CNs:", length(all_plt_cns), "\n")

# Read TREE table (only needed columns)
cat("Reading ENTIRE_TREE.csv ...\n")
tree_cols <- c("CN", "PLT_CN", "INVYR", "STATECD", "CONDID", "STATUSCD",
               "SPCD", "DIA", "HT", "CR", "TPA_UNADJ",
               "VOLCFGRS", "VOLCFNET", "VOLBFNET", "VOLBFGRS")
if (INGROWTH_ENABLED) tree_cols <- c(tree_cols, "PREV_TRE_CN")
fia_trees <- fread(
  file.path(fia_root, "ENTIRE_TREE.csv"),
  select = tree_cols,
  colClasses = c(PLT_CN = "character", CN = "character",
                 if (INGROWTH_ENABLED) c(PREV_TRE_CN = "character") else NULL)
)
cat("  Total tree records:", nrow(fia_trees), "\n")

# Filter to remeasurement plots
fia_trees <- fia_trees[PLT_CN %in% all_plt_cns]
cat("  After filtering:", nrow(fia_trees), "trees\n")

# Read COND table
cat("Reading ENTIRE_COND.csv ...\n")
fia_conds <- fread(
  file.path(fia_root, "ENTIRE_COND.csv"),
  select = c("CN", "PLT_CN", "CONDID", "COND_STATUS_CD", "FORTYPCD",
             "STDAGE", "SITECLCD", "SICOND", "SISP", "BALIVE",
             "SLOPE", "ASPECT"),
  colClasses = c(CN = "character", PLT_CN = "character")
)
fia_conds <- fia_conds[PLT_CN %in% all_plt_cns & COND_STATUS_CD == 1]
setnames(fia_conds, "CN", "COND_CN")
cat("  Forested conditions:", nrow(fia_conds), "\n")

# Merge trees with conditions
fia_trees <- merge(fia_trees, fia_conds,
                   by.x = c("PLT_CN", "CONDID"),
                   by.y = c("PLT_CN", "CONDID"),
                   all.x = FALSE)
cat("  Trees with condition info:", nrow(fia_trees), "\n")

# Assign FVS variants
fia_trees <- merge(fia_trees, variant_map,
                   by.x = "COND_CN", by.y = "STAND_CN", all.x = TRUE)

# Merge FVS STANDINIT_PLOT data (SITE_INDEX, MAX_SDI) if available
# Note: STANDINIT_PLOT.STAND_CN = PLOT.CN (plot-level), not COND.CN
if (!is.null(fvs_plot_init)) {
  fia_trees <- merge(fia_trees, fvs_plot_init,
                     by.x = "PLT_CN", by.y = "STAND_CN", all.x = TRUE)
  cat("  Trees with FVS plot init info:", sum(!is.na(fia_trees$SITE_INDEX)), "\n")
  # Use FVS SITE_INDEX as fallback when SICOND is missing
  fia_trees[is.na(SICOND) & !is.na(SITE_INDEX), SICOND := SITE_INDEX]
}

# Fill missing variants using most common for that state
for (st in unique(fia_trees[is.na(VARIANT), STATECD])) {
  state_vars <- fia_trees[STATECD == st & !is.na(VARIANT)]
  if (nrow(state_vars) > 0) {
    default_var <- names(sort(table(state_vars$VARIANT), decreasing = TRUE))[1]
    fia_trees[STATECD == st & is.na(VARIANT), VARIANT := default_var]
  }
}

# Filter to live trees with diameter
fia_trees <- fia_trees[STATUSCD == 1 & !is.na(DIA) & DIA >= 1.0 & !is.na(TPA_UNADJ)]
cat("  Live trees >= 1 in.:", nrow(fia_trees), "\n\n")

# ==============================================================================
# STEP 2b: Compute Empirical Ingrowth Rates from FIA Data
# ==============================================================================

if (INGROWTH_ENABLED && "PREV_TRE_CN" %in% names(fia_trees)) {
  cat("STEP 2b: Compute Empirical Ingrowth Rates\n")
  cat(strrep("-", 80), "\n\n")

  # Ingrowth trees: live at t2, on remeasurement plots, with no previous tree CN
  # These are trees that grew past the 1-inch threshold between measurements
  ingrowth_trees <- fia_trees[PLT_CN %in% plot_pairs$PLT_CN_t2 &
                               (is.na(PREV_TRE_CN) | PREV_TRE_CN == "" |
                                PREV_TRE_CN == "0" | PREV_TRE_CN == "NA")]
  cat("  Ingrowth trees identified:", nrow(ingrowth_trees), "\n")

  # Compute ingrowth BA and TPA per condition (PLT_CN + CONDID)
  ingrowth_by_cond <- ingrowth_trees[, .(
    ingrowth_TPA = sum(TPA_UNADJ, na.rm = TRUE),
    ingrowth_BA  = sum(TPA_UNADJ * pi * (DIA / 24)^2, na.rm = TRUE),
    ingrowth_VOL_CFGRS = sum(TPA_UNADJ * fifelse(!is.na(VOLCFGRS) & VOLCFGRS > 0, VOLCFGRS, 0), na.rm = TRUE),
    ingrowth_VOL_CFNET = sum(TPA_UNADJ * fifelse(!is.na(VOLCFNET) & VOLCFNET > 0, VOLCFNET, 0), na.rm = TRUE),
    ingrowth_VOL_BFNET = sum(TPA_UNADJ * fifelse(!is.na(VOLBFNET) & VOLBFNET > 0, VOLBFNET, 0), na.rm = TRUE),
    ingrowth_n = .N
  ), by = .(PLT_CN, CONDID)]

  # Merge with variant info to get variant-level rates
  # Use cond variant mapping from fia_trees
  cond_variant <- unique(fia_trees[, .(PLT_CN, CONDID, VARIANT)])
  ingrowth_by_cond <- merge(ingrowth_by_cond, cond_variant,
                             by = c("PLT_CN", "CONDID"), all.x = TRUE)

  # Also get the t1 BA for each condition to compute density-dependent rates
  # t1 BA comes from the t1 trees on the paired plots
  ingrowth_by_cond <- merge(
    ingrowth_by_cond,
    plot_pairs[, .(PLT_CN_t2 = PLT_CN_t2, interval_years)],
    by.x = "PLT_CN", by.y = "PLT_CN_t2", all.x = TRUE, allow.cartesian = TRUE
  )

  # Compute annual ingrowth rate per condition
  ingrowth_by_cond[, `:=`(
    ann_ingrowth_TPA = ingrowth_TPA / pmax(interval_years, 1),
    ann_ingrowth_BA  = ingrowth_BA / pmax(interval_years, 1),
    ann_ingrowth_VOL_CFGRS = ingrowth_VOL_CFGRS / pmax(interval_years, 1),
    ann_ingrowth_VOL_CFNET = ingrowth_VOL_CFNET / pmax(interval_years, 1),
    ann_ingrowth_VOL_BFNET = ingrowth_VOL_BFNET / pmax(interval_years, 1)
  )]

  # Compute median ingrowth rates by variant
  ingrowth_rates <- ingrowth_by_cond[!is.na(VARIANT), .(
    med_ann_TPA   = median(ann_ingrowth_TPA, na.rm = TRUE),
    med_ann_BA    = median(ann_ingrowth_BA, na.rm = TRUE),
    med_ann_VOL_CFGRS = median(ann_ingrowth_VOL_CFGRS, na.rm = TRUE),
    med_ann_VOL_CFNET = median(ann_ingrowth_VOL_CFNET, na.rm = TRUE),
    med_ann_VOL_BFNET = median(ann_ingrowth_VOL_BFNET, na.rm = TRUE),
    n_conds = .N
  ), by = VARIANT]

  # Also compute an overall rate as fallback
  overall_ingrowth <- ingrowth_by_cond[!is.na(VARIANT), .(
    med_ann_TPA   = median(ann_ingrowth_TPA, na.rm = TRUE),
    med_ann_BA    = median(ann_ingrowth_BA, na.rm = TRUE),
    med_ann_VOL_CFGRS = median(ann_ingrowth_VOL_CFGRS, na.rm = TRUE),
    med_ann_VOL_CFNET = median(ann_ingrowth_VOL_CFNET, na.rm = TRUE),
    med_ann_VOL_BFNET = median(ann_ingrowth_VOL_BFNET, na.rm = TRUE)
  )]

  cat("  Ingrowth rates by variant:\n")
  for (i in 1:nrow(ingrowth_rates)) {
    cat(sprintf("    %s: TPA=%.1f/yr, BA=%.2f/yr, VOL=%.1f/yr (n=%d)\n",
                ingrowth_rates$VARIANT[i],
                ingrowth_rates$med_ann_TPA[i],
                ingrowth_rates$med_ann_BA[i],
                ingrowth_rates$med_ann_VOL_CFGRS[i],
                ingrowth_rates$n_conds[i]))
  }
  cat(sprintf("  Overall: TPA=%.1f/yr, BA=%.2f/yr, VOL=%.1f/yr\n\n",
              overall_ingrowth$med_ann_TPA,
              overall_ingrowth$med_ann_BA,
              overall_ingrowth$med_ann_VOL_CFGRS))

  # Store as named lists for fast lookup in projection
  ingrowth_lookup <- list()
  for (i in 1:nrow(ingrowth_rates)) {
    v <- ingrowth_rates$VARIANT[i]
    ingrowth_lookup[[v]] <- as.list(ingrowth_rates[i])
  }
  ingrowth_lookup[["OVERALL"]] <- as.list(overall_ingrowth)

  # Clean up
  rm(ingrowth_trees, ingrowth_by_cond, cond_variant)
} else {
  INGROWTH_ENABLED <- FALSE
  ingrowth_lookup <- NULL
  cat("Ingrowth model: SKIPPED (PREV_TRE_CN not available)\n\n")
}

# ==============================================================================
# STEP 2c: Load Raster Lookup (ClimateSI + Emmerson SDIMAX)
# ==============================================================================

raster_lookup <- NULL
if ((CLIMATE_SI_ENABLED || SDIMAX_RASTER_ENABLED) && file.exists(RASTER_LOOKUP_CSV)) {
  cat("STEP 2c: Load Raster Lookup (ClimateSI + SDIMAX)\n")
  cat(strrep("-", 80), "\n\n")

  raster_lookup <- fread(RASTER_LOOKUP_CSV,
                          colClasses = c(PLT_CN = "character"))
  cat("  Raster lookup records:", nrow(raster_lookup), "\n")

  if (CLIMATE_SI_ENABLED && "ClimateSI_ft" %in% names(raster_lookup)) {
    valid_si <- raster_lookup[!is.na(ClimateSI_ft) & ClimateSI_ft > 0]
    cat("  ClimateSI values:", nrow(valid_si),
        "| mean:", round(mean(valid_si$ClimateSI_ft), 1), "ft",
        "| range:", round(min(valid_si$ClimateSI_ft), 1), "-",
        round(max(valid_si$ClimateSI_ft), 1), "\n")
  } else {
    CLIMATE_SI_ENABLED <- FALSE
    cat("  ClimateSI: column not found, DISABLED\n")
  }

  if (SDIMAX_RASTER_ENABLED && "SDIMAX_imperial" %in% names(raster_lookup)) {
    valid_sdi <- raster_lookup[!is.na(SDIMAX_imperial) & SDIMAX_imperial > 0]
    cat("  SDIMAX (imperial) values:", nrow(valid_sdi),
        "| mean:", round(mean(valid_sdi$SDIMAX_imperial), 1),
        "| range:", round(min(valid_sdi$SDIMAX_imperial), 1), "-",
        round(max(valid_sdi$SDIMAX_imperial), 1), "\n")
  } else {
    SDIMAX_RASTER_ENABLED <- FALSE
    cat("  SDIMAX: column not found, DISABLED\n")
  }

  setkey(raster_lookup, PLT_CN)
  cat("\n")
} else {
  if (CLIMATE_SI_ENABLED || SDIMAX_RASTER_ENABLED) {
    cat("  Raster lookup not found:", RASTER_LOOKUP_CSV, "\n")
    cat("  Falling back to FIA SICOND and calibrated SDIMAX\n\n")
  }
  CLIMATE_SI_ENABLED <- FALSE
  SDIMAX_RASTER_ENABLED <- FALSE
}

# ==============================================================================
# STEP 3: Compute Observed Stand Metrics by Condition
# ==============================================================================

cat("STEP 3: Compute Observed Metrics\n")
cat(strrep("-", 80), "\n\n")

# Compute metrics by PLT_CN + CONDID (this covers both t1 and t2 plots)
cat("Computing stand metrics for all conditions...\n")
cond_metrics <- fia_trees[, .(
  TPA = sum(TPA_UNADJ, na.rm = TRUE),
  BA  = sum(TPA_UNADJ * pi * (DIA / 24)^2, na.rm = TRUE),
  VOL_CFGRS = sum(TPA_UNADJ * fifelse(!is.na(VOLCFGRS) & VOLCFGRS > 0, VOLCFGRS, 0), na.rm = TRUE),
  VOL_CFNET = sum(TPA_UNADJ * fifelse(!is.na(VOLCFNET) & VOLCFNET > 0, VOLCFNET, 0), na.rm = TRUE),
  VOL_BFNET = sum(TPA_UNADJ * fifelse(!is.na(VOLBFNET) & VOLBFNET > 0, VOLBFNET, 0), na.rm = TRUE),
  N_trees = .N,
  VARIANT = VARIANT[1],
  SICOND  = SICOND[1],
  FORTYPCD = FORTYPCD[1],
  SLOPE_cond = SLOPE[1],
  ASPECT_cond = ASPECT[1],
  FVS_SITE_INDEX = if ("SITE_INDEX" %in% names(fia_trees)) SITE_INDEX[1] else NA_real_,
  FVS_MAX_SDI    = if ("MAX_SDI" %in% names(fia_trees)) MAX_SDI[1] else NA_real_,
  FVS_MAX_SDI_FIA = if ("MAX_SDI_FIA" %in% names(fia_trees)) MAX_SDI_FIA[1] else NA_real_
), by = .(PLT_CN, CONDID)]

cond_metrics[, QMD := fifelse(TPA > 0 & BA > 0, sqrt(BA / TPA * 576 / pi), NA_real_)]
cond_metrics[, SDI := fifelse(!is.na(QMD) & QMD > 0, TPA * (QMD / 10)^1.605, NA_real_)]

# Top height: TPA-weighted mean height of top 40 TPA (Lord method)
cat("Computing observed top height per condition...\n")
ht_top_obs <- fia_trees[!is.na(HT) & HT > 0 & !is.na(TPA_UNADJ) & TPA_UNADJ > 0,
  {
    dt_ht <- .SD[order(-HT)]
    cum_tpa <- cumsum(dt_ht$TPA_UNADJ)
    top_idx <- cum_tpa <= 40
    if (sum(top_idx) == 0) top_idx[1] <- TRUE
    list(HT_top = weighted.mean(dt_ht$HT[top_idx], dt_ht$TPA_UNADJ[top_idx]))
  }, by = .(PLT_CN, CONDID)]
cond_metrics <- merge(cond_metrics, ht_top_obs, by = c("PLT_CN", "CONDID"), all.x = TRUE)
cat("  Top height computed for", sum(!is.na(cond_metrics$HT_top)), "/",
    nrow(cond_metrics), "conditions\n")

cat("  Conditions computed:", nrow(cond_metrics), "\n")

# Build matched condition pairs using plot_pairs
cat("Matching t1/t2 condition pairs...\n")

# Time1 metrics (keyed by PLT_CN_t1)
t1_metrics <- cond_metrics[PLT_CN %in% plot_pairs$PLT_CN_t1]
setnames(t1_metrics, c("PLT_CN", "TPA", "BA", "VOL_CFGRS", "VOL_CFNET", "VOL_BFNET",
                        "N_trees", "QMD", "SDI", "HT_top"),
         c("PLT_CN_t1", "TPA_t1", "BA_t1", "VOL_CFGRS_t1", "VOL_CFNET_t1", "VOL_BFNET_t1",
           "N_trees_t1", "QMD_t1", "SDI_t1", "HT_top_t1"))

# Time2 metrics (keyed by PLT_CN_t2)
t2_metrics <- cond_metrics[PLT_CN %in% plot_pairs$PLT_CN_t2]
setnames(t2_metrics, c("PLT_CN", "TPA", "BA", "VOL_CFGRS", "VOL_CFNET", "VOL_BFNET",
                        "N_trees", "QMD", "SDI", "HT_top"),
         c("PLT_CN_t2", "TPA_t2", "BA_t2", "VOL_CFGRS_t2", "VOL_CFNET_t2", "VOL_BFNET_t2",
           "N_trees_t2", "QMD_t2", "SDI_t2", "HT_top_t2"))

# Join through plot_pairs to get matched t1/t2 conditions
matched <- merge(
  plot_pairs[, .(PLT_CN_t1, PLT_CN_t2, interval_years, ELEV_t2)],
  t1_metrics[, .(PLT_CN_t1, CONDID, TPA_t1, BA_t1, QMD_t1, SDI_t1,
                  VOL_CFGRS_t1, VOL_CFNET_t1, VOL_BFNET_t1, HT_top_t1,
                  VARIANT, SICOND, FORTYPCD, SLOPE_cond, ASPECT_cond,
                  FVS_SITE_INDEX, FVS_MAX_SDI, FVS_MAX_SDI_FIA)],
  by = "PLT_CN_t1", allow.cartesian = TRUE
)
matched <- merge(
  matched,
  t2_metrics[, .(PLT_CN_t2, CONDID, TPA_t2, BA_t2, QMD_t2, SDI_t2,
                  VOL_CFGRS_t2, VOL_CFNET_t2, VOL_BFNET_t2, HT_top_t2)],
  by = c("PLT_CN_t2", "CONDID")
)

# Remove conditions with zero or missing metrics
matched <- matched[TPA_t1 > 0 & BA_t1 > 0 & TPA_t2 > 0 & BA_t2 > 0 &
                   !is.na(QMD_t1) & !is.na(QMD_t2)]

cat("  Matched condition pairs:", nrow(matched), "\n")
cat("  Variants present:", paste(sort(unique(matched$VARIANT)), collapse = ", "), "\n")
cat("  Mean interval:", round(mean(matched$interval_years), 2), "years\n\n")

# Join raster lookup values (ClimateSI, SDIMAX) to matched conditions
if (!is.null(raster_lookup)) {
  cat("Joining raster lookup to matched conditions...\n")
  matched <- merge(matched, raster_lookup[, .(PLT_CN, ClimateSI_ft, SDIMAX_imperial)],
                   by.x = "PLT_CN_t1", by.y = "PLT_CN", all.x = TRUE)

  if (CLIMATE_SI_ENABLED) {
    n_si <- sum(!is.na(matched$ClimateSI_ft) & matched$ClimateSI_ft > 0)
    cat("  ClimateSI matched:", n_si, "/", nrow(matched),
        sprintf("(%.1f%%)\n", 100 * n_si / nrow(matched)))
  }
  if (SDIMAX_RASTER_ENABLED) {
    n_sdi <- sum(!is.na(matched$SDIMAX_imperial) & matched$SDIMAX_imperial > 0)
    cat("  SDIMAX matched:", n_sdi, "/", nrow(matched),
        sprintf("(%.1f%%)\n", 100 * n_sdi / nrow(matched)))
  }
  cat("\n")
}

# Save observed changes
fwrite(matched, file.path(output_root, "intermediate/observed_changes.csv"))

# ==============================================================================
# STEP 4: Run Tree-Level Projections
# ==============================================================================

cat("STEP 4: Run Tree-Level Projections\n")
cat(strrep("-", 80), "\n\n")

# Index trees by PLT_CN for fast lookup
setkey(fia_trees, PLT_CN, CONDID)

# Process by variant for efficiency
variants_in_data <- sort(unique(matched$VARIANT))
projection_results <- vector("list", nrow(matched))
processed <- 0
failed <- 0

for (var in variants_in_data) {
  var_rows <- which(matched$VARIANT == var)
  n_var <- length(var_rows)
  cat("  Projecting", toupper(var), ":", n_var, "conditions... ")

  params <- variant_params[[var]]
  if (is.null(params)) {
    cat("SKIP (no params)\n")
    failed <- failed + n_var
    next
  }

  t_start <- proc.time()

  for (idx in var_rows) {
    row <- matched[idx]

    # Get time1 trees for this condition
    t1_trees <- fia_trees[PLT_CN == row$PLT_CN_t1 & CONDID == row$CONDID]

    if (nrow(t1_trees) == 0) {
      projection_results[[idx]] <- data.table(
        PLT_CN_t1 = row$PLT_CN_t1, PLT_CN_t2 = row$PLT_CN_t2,
        CONDID = row$CONDID, VARIANT = var,
        interval_years = row$interval_years,
        TPA_pred_calib = NA_real_, BA_pred_calib = NA_real_,
        QMD_pred_calib = NA_real_, SDI_pred_calib = NA_real_,
        VOL_CFGRS_pred_calib = NA_real_,
        VOL_CFNET_pred_calib = NA_real_,
        VOL_BFNET_pred_calib = NA_real_,
        VOL_CFGRS_gross_calib = NA_real_,
        HT_top_calib = NA_real_,
        BA_pred_calib_lo = NA_real_, BA_pred_calib_hi = NA_real_,
        VOL_CFGRS_pred_calib_lo = NA_real_, VOL_CFGRS_pred_calib_hi = NA_real_,
        HT_top_calib_lo = NA_real_, HT_top_calib_hi = NA_real_,
        TPA_pred_default = NA_real_, BA_pred_default = NA_real_,
        QMD_pred_default = NA_real_, SDI_pred_default = NA_real_,
        VOL_CFGRS_pred_default = NA_real_,
        VOL_CFNET_pred_default = NA_real_,
        VOL_BFNET_pred_default = NA_real_,
        VOL_CFGRS_gross_default = NA_real_,
        HT_top_default = NA_real_
      )
      failed <- failed + 1
      next
    }

    # Site index: prefer ClimateSI (raster), fall back to SICOND, default 65
    si <- NA_real_
    if (CLIMATE_SI_ENABLED && "ClimateSI_ft" %in% names(row) &&
        !is.na(row$ClimateSI_ft) && row$ClimateSI_ft > 0) {
      si <- row$ClimateSI_ft
    } else {
      si <- row$SICOND
    }
    if (is.na(si) | si <= 0 | si > 200) si <- 65

    # SDIMAX: prefer Emmerson raster, fall back to NA
    sdimax_val <- NA_real_
    if (SDIMAX_RASTER_ENABLED && "SDIMAX_imperial" %in% names(row) &&
        !is.na(row$SDIMAX_imperial) && row$SDIMAX_imperial > 0) {
      sdimax_val <- row$SDIMAX_imperial
    }

    # Calibrated projection (median posterior)
    proj_args <- list(trees = t1_trees, params = params,
                      interval_years = row$interval_years, si = si,
                      slope_cond = row$SLOPE_cond,
                      aspect_cond = row$ASPECT_cond,
                      elev_cond = row$ELEV_t2,
                      sdimax_raster_val = sdimax_val)
    null_err <- list(TPA_pred = NA_real_, BA_pred = NA_real_,
                      QMD_pred = NA_real_, SDI_pred = NA_real_,
                      VOL_CFGRS_pred = NA_real_, VOL_CFNET_pred = NA_real_,
                      VOL_BFNET_pred = NA_real_,
                      VOL_CFGRS_gross = NA_real_, HT_top = NA_real_)
    calib <- tryCatch(do.call(project_condition_calibrated, proj_args),
                       error = function(e) null_err)

    # Default projection
    default <- tryCatch(
      project_condition_default(t1_trees, row$interval_years, variant_code = var),
      error = function(e) null_err)

    # CI bracket projections (Q10/Q90 parameter swap)
    ba_lo <- NA_real_; ba_hi <- NA_real_
    vol_lo <- NA_real_; vol_hi <- NA_real_
    ht_lo <- NA_real_; ht_hi <- NA_real_
    if (CI_BRACKETS_ENABLED && !is.null(params$dg_lo)) {
      # Build "low growth" params: Q10 DG betas, Q90 mortality intercept (more death)
      params_lo <- params
      params_lo$dg <- params$dg_lo
      params_lo$hd <- params$hd_lo
      if (!is.null(params$mort_lo)) {
        params_lo$mort <- params$mort_hi  # higher mort intercept = more mortality = lower stand
      }
      calib_lo <- tryCatch(do.call(project_condition_calibrated,
                                    c(list(trees = t1_trees, params = params_lo,
                                           interval_years = row$interval_years, si = si,
                                           slope_cond = row$SLOPE_cond,
                                           aspect_cond = row$ASPECT_cond,
                                           elev_cond = row$ELEV_t2,
                                           sdimax_raster_val = sdimax_val))),
                            error = function(e) null_err)
      # Build "high growth" params: Q90 DG betas, Q10 mortality intercept (less death)
      params_hi <- params
      params_hi$dg <- params$dg_hi
      params_hi$hd <- params$hd_hi
      if (!is.null(params$mort_hi)) {
        params_hi$mort <- params$mort_lo  # lower mort intercept = less mortality = higher stand
      }
      calib_hi <- tryCatch(do.call(project_condition_calibrated,
                                    c(list(trees = t1_trees, params = params_hi,
                                           interval_years = row$interval_years, si = si,
                                           slope_cond = row$SLOPE_cond,
                                           aspect_cond = row$ASPECT_cond,
                                           elev_cond = row$ELEV_t2,
                                           sdimax_raster_val = sdimax_val))),
                            error = function(e) null_err)
      # Ensure lo < hi by taking min/max
      ba_lo  <- min(calib_lo$BA_pred, calib_hi$BA_pred, na.rm = TRUE)
      ba_hi  <- max(calib_lo$BA_pred, calib_hi$BA_pred, na.rm = TRUE)
      vol_lo <- min(calib_lo$VOL_CFGRS_pred, calib_hi$VOL_CFGRS_pred, na.rm = TRUE)
      vol_hi <- max(calib_lo$VOL_CFGRS_pred, calib_hi$VOL_CFGRS_pred, na.rm = TRUE)
      ht_lo  <- min(calib_lo$HT_top, calib_hi$HT_top, na.rm = TRUE)
      ht_hi  <- max(calib_lo$HT_top, calib_hi$HT_top, na.rm = TRUE)
      # Guard against Inf from all-NA inputs
      if (!is.finite(ba_lo))  ba_lo  <- NA_real_
      if (!is.finite(ba_hi))  ba_hi  <- NA_real_
      if (!is.finite(vol_lo)) vol_lo <- NA_real_
      if (!is.finite(vol_hi)) vol_hi <- NA_real_
      if (!is.finite(ht_lo))  ht_lo  <- NA_real_
      if (!is.finite(ht_hi))  ht_hi  <- NA_real_
    }

    projection_results[[idx]] <- data.table(
      PLT_CN_t1 = row$PLT_CN_t1, PLT_CN_t2 = row$PLT_CN_t2,
      CONDID = row$CONDID, VARIANT = var,
      interval_years = row$interval_years,
      TPA_pred_calib = calib$TPA_pred, BA_pred_calib = calib$BA_pred,
      QMD_pred_calib = calib$QMD_pred, SDI_pred_calib = calib$SDI_pred,
      VOL_CFGRS_pred_calib = calib$VOL_CFGRS_pred,
      VOL_CFNET_pred_calib = calib$VOL_CFNET_pred,
      VOL_BFNET_pred_calib = calib$VOL_BFNET_pred,
      VOL_CFGRS_gross_calib = calib$VOL_CFGRS_gross,
      HT_top_calib = calib$HT_top,
      BA_pred_calib_lo = ba_lo, BA_pred_calib_hi = ba_hi,
      VOL_CFGRS_pred_calib_lo = vol_lo, VOL_CFGRS_pred_calib_hi = vol_hi,
      HT_top_calib_lo = ht_lo, HT_top_calib_hi = ht_hi,
      TPA_pred_default = default$TPA_pred, BA_pred_default = default$BA_pred,
      QMD_pred_default = default$QMD_pred, SDI_pred_default = default$SDI_pred,
      VOL_CFGRS_pred_default = default$VOL_CFGRS_pred,
      VOL_CFNET_pred_default = default$VOL_CFNET_pred,
      VOL_BFNET_pred_default = default$VOL_BFNET_pred,
      VOL_CFGRS_gross_default = default$VOL_CFGRS_gross,
      HT_top_default = default$HT_top
    )
    processed <- processed + 1
  }

  elapsed <- (proc.time() - t_start)[3]
  cat(sprintf("%.1f sec (%.1f cond/sec)\n", elapsed, n_var / max(elapsed, 0.001)))
}

cat("\nProjections completed:", processed, "successful,", failed, "failed\n\n")

# Combine results
projections_dt <- rbindlist(projection_results[!sapply(projection_results, is.null)],
                            fill = TRUE)

# ==============================================================================
# STEP 5: Compute Validation Statistics
# ==============================================================================

cat("STEP 5: Compute Validation Statistics\n")
cat(strrep("-", 80), "\n\n")

# Merge projections with observed time2 metrics
# Select only columns that exist in matched (FVS columns may be missing)
merge_cols <- c("PLT_CN_t1", "PLT_CN_t2", "CONDID", "VARIANT", "interval_years",
                "TPA_t1", "BA_t1", "QMD_t1", "SDI_t1",
                "VOL_CFGRS_t1", "VOL_CFNET_t1", "VOL_BFNET_t1",
                "TPA_t2", "BA_t2", "QMD_t2", "SDI_t2",
                "VOL_CFGRS_t2", "VOL_CFNET_t2", "VOL_BFNET_t2",
                "HT_top_t1", "HT_top_t2")
opt_cols <- c("FVS_SITE_INDEX", "FVS_MAX_SDI", "FVS_MAX_SDI_FIA",
              "ClimateSI_ft", "SDIMAX_imperial")
merge_cols <- c(merge_cols, intersect(opt_cols, names(matched)))

validation_data <- merge(
  matched[, ..merge_cols],
  projections_dt,
  by = c("PLT_CN_t1", "PLT_CN_t2", "CONDID", "VARIANT", "interval_years")
)

# Compute NET PAI: (VOL_t2 - VOL_t1) / interval = observed stand change per year
# Net PAI confounds growth + mortality + ingrowth (what actually happened)
validation_data[, PAI_net_obs := fifelse(interval_years > 0,
  (VOL_CFGRS_t2 - VOL_CFGRS_t1) / interval_years, NA_real_)]
validation_data[, PAI_net_pred_calib := fifelse(interval_years > 0,
  (VOL_CFGRS_pred_calib - VOL_CFGRS_t1) / interval_years, NA_real_)]
validation_data[, PAI_net_pred_default := fifelse(interval_years > 0,
  (VOL_CFGRS_pred_default - VOL_CFGRS_t1) / interval_years, NA_real_)]

# Compute GROSS PAI: growth on all trees ignoring mortality
# Gross PAI uses VOL_CFGRS_gross (sum of TPA_UNADJ * VOL_pred, no mortality applied)
validation_data[, PAI_gross_pred_calib := fifelse(interval_years > 0,
  (VOL_CFGRS_gross_calib - VOL_CFGRS_t1) / interval_years, NA_real_)]
validation_data[, PAI_gross_pred_default := fifelse(interval_years > 0,
  (VOL_CFGRS_gross_default - VOL_CFGRS_t1) / interval_years, NA_real_)]

# Legacy PAI columns (kept for backward compatibility)
validation_data[, PAI_CFGRS_obs := PAI_net_obs]
validation_data[, PAI_CFGRS_pred_calib := PAI_net_pred_calib]
validation_data[, PAI_CFGRS_pred_default := PAI_net_pred_default]
validation_data[, PAI_CFNET_obs := fifelse(interval_years > 0,
  (VOL_CFNET_t2 - VOL_CFNET_t1) / interval_years, NA_real_)]
validation_data[, PAI_CFNET_pred_calib := fifelse(interval_years > 0,
  (VOL_CFNET_pred_calib - VOL_CFNET_t1) / interval_years, NA_real_)]
validation_data[, PAI_CFNET_pred_default := fifelse(interval_years > 0,
  (VOL_CFNET_pred_default - VOL_CFNET_t1) / interval_years, NA_real_)]
validation_data[, PAI_BFNET_obs := fifelse(interval_years > 0,
  (VOL_BFNET_t2 - VOL_BFNET_t1) / interval_years, NA_real_)]
validation_data[, PAI_BFNET_pred_calib := fifelse(interval_years > 0,
  (VOL_BFNET_pred_calib - VOL_BFNET_t1) / interval_years, NA_real_)]
validation_data[, PAI_BFNET_pred_default := fifelse(interval_years > 0,
  (VOL_BFNET_pred_default - VOL_BFNET_t1) / interval_years, NA_real_)]

# Remove NAs
validation_data <- validation_data[!is.na(BA_pred_calib) & !is.na(BA_pred_default)]
cat("Validation pairs:", nrow(validation_data), "\n")

# Save intermediate
tag_suffix <- if (OUTPUT_TAG != "") paste0("_", OUTPUT_TAG) else ""
fwrite(validation_data, file.path(output_root, paste0("intermediate/validation_data", tag_suffix, ".csv")))

# --- Validation metrics ---
calc_rmse <- function(pred, obs) sqrt(mean((pred - obs)^2, na.rm = TRUE))
calc_bias <- function(pred, obs) mean(pred - obs, na.rm = TRUE)
calc_mae <- function(pred, obs) mean(abs(pred - obs), na.rm = TRUE)
calc_bias_pct <- function(pred, obs) {
  valid <- !is.na(pred) & !is.na(obs) & obs != 0
  if (sum(valid) == 0) return(NA_real_)
  100 * mean((pred[valid] - obs[valid]) / obs[valid])
}
calc_rmse_pct <- function(pred, obs) {
  obs_mean <- mean(obs, na.rm = TRUE)
  if (is.na(obs_mean) || abs(obs_mean) < 1e-10) return(NA_real_)
  100 * calc_rmse(pred, obs) / abs(obs_mean)
}
calc_mae_pct <- function(pred, obs) {
  obs_mean <- mean(obs, na.rm = TRUE)
  if (is.na(obs_mean) || abs(obs_mean) < 1e-10) return(NA_real_)
  100 * calc_mae(pred, obs) / abs(obs_mean)
}
calc_r2 <- function(pred, obs) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  max(1 - ss_res / ss_tot, 0)
}
calc_equiv <- function(pred, obs, tol = 0.20) {
  valid <- !is.na(pred) & !is.na(obs) & obs > 0
  if (sum(valid) == 0) return(NA_real_)
  100 * mean(abs((pred[valid] - obs[valid]) / obs[valid]) <= tol)
}
# Willmott's agreement index (d)
calc_willmott_d <- function(pred, obs) {
  valid <- !is.na(pred) & !is.na(obs)
  if (sum(valid) < 3) return(NA_real_)
  p <- pred[valid]; o <- obs[valid]
  obs_mean <- mean(o)
  ss_res <- sum((o - p)^2)
  ss_pot <- sum((abs(p - obs_mean) + abs(o - obs_mean))^2)
  if (ss_pot == 0) return(NA_real_)
  1 - ss_res / ss_pot
}
# Helper: compute full metric suite for a pred/obs pair
compute_metrics <- function(pred, obs, prefix, suffix) {
  setNames(
    list(
      calc_rmse(pred, obs),
      calc_mae(pred, obs),
      calc_bias(pred, obs),
      calc_bias_pct(pred, obs),
      calc_rmse_pct(pred, obs),
      calc_mae_pct(pred, obs),
      calc_r2(pred, obs),
      calc_equiv(pred, obs),
      calc_willmott_d(pred, obs),
      mean(obs, na.rm = TRUE),
      mean(pred, na.rm = TRUE)
    ),
    paste0(prefix, c("_RMSE_", "_MAE_", "_bias_", "_bias_pct_", "_RMSE_pct_",
                      "_MAE_pct_", "_r2_", "_equiv_", "_d_",
                      "_obs_mean_", "_pred_mean_"), suffix)
  )
}

# --- Bootstrap CI for validation statistics ---
# Returns 95% CI on RMSE, R2, bias, MAE via nonparametric bootstrap
bootstrap_stat_ci <- function(pred, obs, n_boot = BOOTSTRAP_N, seed = 42) {
  n <- length(pred)
  if (n < 20) return(NULL)
  set.seed(seed)

  boot_rmse <- boot_r2 <- boot_bias <- boot_mae <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    p <- pred[idx]; o <- obs[idx]
    resid <- p - o
    boot_rmse[b] <- sqrt(mean(resid^2))
    boot_bias[b] <- mean(resid)
    boot_mae[b]  <- mean(abs(resid))
    ss_res <- sum((o - p)^2)
    ss_tot <- sum((o - mean(o))^2)
    boot_r2[b] <- if (ss_tot > 0) max(1 - ss_res / ss_tot, 0) else NA_real_
  }
  list(
    RMSE_lo = quantile(boot_rmse, 0.025, na.rm = TRUE),
    RMSE_hi = quantile(boot_rmse, 0.975, na.rm = TRUE),
    r2_lo   = quantile(boot_r2,   0.025, na.rm = TRUE),
    r2_hi   = quantile(boot_r2,   0.975, na.rm = TRUE),
    bias_lo = quantile(boot_bias, 0.025, na.rm = TRUE),
    bias_hi = quantile(boot_bias, 0.975, na.rm = TRUE),
    MAE_lo  = quantile(boot_mae,  0.025, na.rm = TRUE),
    MAE_hi  = quantile(boot_mae,  0.975, na.rm = TRUE)
  )
}

# Wrapper: compute bootstrap CI and return named list with prefix/suffix
compute_boot_ci <- function(pred, obs, prefix, suffix) {
  if (!BOOTSTRAP_CI_ENABLED) return(list())
  ci <- bootstrap_stat_ci(pred, obs)
  if (is.null(ci)) return(list())
  setNames(
    list(ci$RMSE_lo, ci$RMSE_hi, ci$r2_lo, ci$r2_hi,
         ci$bias_lo, ci$bias_hi, ci$MAE_lo, ci$MAE_hi),
    paste0(prefix, c("_RMSE_lo_", "_RMSE_hi_", "_r2_lo_", "_r2_hi_",
                      "_bias_lo_", "_bias_hi_", "_MAE_lo_", "_MAE_hi_"), suffix)
  )
}

# --- Helper: compute all metrics for a given data subset ---
compute_variant_stats <- function(vd, label) {
  stats <- list(VARIANT = label, n_conditions = nrow(vd))

  # Core stand attributes: BA, TPA, QMD, SDI (full metric suite)
  core_attrs <- list(
    BA  = list(pred_c = "BA_pred_calib",  pred_d = "BA_pred_default",  obs = "BA_t2"),
    TPA = list(pred_c = "TPA_pred_calib", pred_d = "TPA_pred_default", obs = "TPA_t2"),
    QMD = list(pred_c = "QMD_pred_calib", pred_d = "QMD_pred_default", obs = "QMD_t2"),
    SDI = list(pred_c = "SDI_pred_calib", pred_d = "SDI_pred_default", obs = "SDI_t2")
  )
  for (attr_nm in names(core_attrs)) {
    a <- core_attrs[[attr_nm]]
    stats <- c(stats, compute_metrics(vd[[a$pred_c]], vd[[a$obs]], attr_nm, "calib"))
    stats <- c(stats, compute_metrics(vd[[a$pred_d]], vd[[a$obs]], attr_nm, "default"))
    # Bootstrap CI on calibrated stats
    stats <- c(stats, compute_boot_ci(vd[[a$pred_c]], vd[[a$obs]], attr_nm, "calib"))
  }

  # Volume attributes: CFGRS, CFNET, BFNET (subset to positive observed volume)
  vol_attrs <- list(
    VOL_CFGRS = list(pred_c = "VOL_CFGRS_pred_calib", pred_d = "VOL_CFGRS_pred_default",
                     obs = "VOL_CFGRS_t2"),
    VOL_CFNET = list(pred_c = "VOL_CFNET_pred_calib", pred_d = "VOL_CFNET_pred_default",
                     obs = "VOL_CFNET_t2"),
    VOL_BFNET = list(pred_c = "VOL_BFNET_pred_calib", pred_d = "VOL_BFNET_pred_default",
                     obs = "VOL_BFNET_t2")
  )
  for (vol_nm in names(vol_attrs)) {
    va <- vol_attrs[[vol_nm]]
    if (all(c(va$pred_c, va$pred_d, va$obs) %in% names(vd))) {
      vd_v <- vd[get(va$obs) > 0 & !is.na(get(va$pred_c)) & !is.na(get(va$pred_d))]
      if (nrow(vd_v) >= 10) {
        stats[[paste0("n_", vol_nm, "_conditions")]] <- nrow(vd_v)
        stats <- c(stats, compute_metrics(vd_v[[va$pred_c]], vd_v[[va$obs]], vol_nm, "calib"))
        stats <- c(stats, compute_metrics(vd_v[[va$pred_d]], vd_v[[va$obs]], vol_nm, "default"))
        stats <- c(stats, compute_boot_ci(vd_v[[va$pred_c]], vd_v[[va$obs]], vol_nm, "calib"))
      }
    }
  }

  # PAI attributes: CFGRS, CFNET, BFNET (legacy net PAI)
  pai_attrs <- list(
    PAI_CFGRS = list(pred_c = "PAI_CFGRS_pred_calib", pred_d = "PAI_CFGRS_pred_default",
                     obs = "PAI_CFGRS_obs", vol_t1 = "VOL_CFGRS_t1", vol_t2 = "VOL_CFGRS_t2"),
    PAI_CFNET = list(pred_c = "PAI_CFNET_pred_calib", pred_d = "PAI_CFNET_pred_default",
                     obs = "PAI_CFNET_obs", vol_t1 = "VOL_CFNET_t1", vol_t2 = "VOL_CFNET_t2"),
    PAI_BFNET = list(pred_c = "PAI_BFNET_pred_calib", pred_d = "PAI_BFNET_pred_default",
                     obs = "PAI_BFNET_obs", vol_t1 = "VOL_BFNET_t1", vol_t2 = "VOL_BFNET_t2")
  )
  for (pai_nm in names(pai_attrs)) {
    pa <- pai_attrs[[pai_nm]]
    if (all(c(pa$pred_c, pa$pred_d, pa$obs, pa$vol_t1, pa$vol_t2) %in% names(vd))) {
      vd_p <- vd[!is.na(get(pa$obs)) & !is.na(get(pa$pred_c)) & !is.na(get(pa$pred_d)) &
                 get(pa$vol_t1) > 0 & get(pa$vol_t2) > 0]
      if (nrow(vd_p) >= 10) {
        stats[[paste0("n_", pai_nm, "_conditions")]] <- nrow(vd_p)
        stats <- c(stats, compute_metrics(vd_p[[pa$pred_c]], vd_p[[pa$obs]], pai_nm, "calib"))
        stats <- c(stats, compute_metrics(vd_p[[pa$pred_d]], vd_p[[pa$obs]], pai_nm, "default"))
      }
    }
  }

  # Gross vs Net PAI decomposition (CFGRS only)
  # Net PAI = (VOL_t2 - VOL_t1) / interval (observed stand change)
  # Gross PAI = growth on all trees ignoring mortality
  if (all(c("PAI_net_obs", "PAI_net_pred_calib", "PAI_net_pred_default",
            "PAI_gross_pred_calib", "PAI_gross_pred_default",
            "VOL_CFGRS_t1", "VOL_CFGRS_t2") %in% names(vd))) {
    vd_pai <- vd[!is.na(PAI_net_obs) & !is.na(PAI_net_pred_calib) &
                  !is.na(PAI_gross_pred_calib) & VOL_CFGRS_t1 > 0 & VOL_CFGRS_t2 > 0]
    if (nrow(vd_pai) >= 10) {
      stats[["n_PAI_net_conditions"]] <- nrow(vd_pai)
      stats <- c(stats, compute_metrics(vd_pai$PAI_net_pred_calib, vd_pai$PAI_net_obs,
                                         "PAI_net", "calib"))
      stats <- c(stats, compute_metrics(vd_pai$PAI_net_pred_default, vd_pai$PAI_net_obs,
                                         "PAI_net", "default"))
      stats[["n_PAI_gross_conditions"]] <- nrow(vd_pai)
      # For gross PAI, observed = net obs (we cannot observe gross directly from FIA)
      # So we compare gross predicted means to give context
      stats[["PAI_gross_pred_mean_calib"]]   <- mean(vd_pai$PAI_gross_pred_calib, na.rm = TRUE)
      stats[["PAI_gross_pred_mean_default"]]  <- mean(vd_pai$PAI_gross_pred_default, na.rm = TRUE)
      stats[["PAI_net_obs_mean"]]             <- mean(vd_pai$PAI_net_obs, na.rm = TRUE)
      stats[["PAI_net_pred_mean_calib"]]      <- mean(vd_pai$PAI_net_pred_calib, na.rm = TRUE)
      stats[["PAI_net_pred_mean_default"]]    <- mean(vd_pai$PAI_net_pred_default, na.rm = TRUE)
    }
  }

  # Top height (Lord method, top 40 TPA)
  if (all(c("HT_top_t2", "HT_top_calib", "HT_top_default") %in% names(vd))) {
    vd_ht <- vd[!is.na(HT_top_t2) & HT_top_t2 > 0 &
                 !is.na(HT_top_calib) & !is.na(HT_top_default)]
    if (nrow(vd_ht) >= 10) {
      stats[["n_HT_top_conditions"]] <- nrow(vd_ht)
      stats <- c(stats, compute_metrics(vd_ht$HT_top_calib, vd_ht$HT_top_t2, "HT_top", "calib"))
      stats <- c(stats, compute_metrics(vd_ht$HT_top_default, vd_ht$HT_top_t2, "HT_top", "default"))
      stats <- c(stats, compute_boot_ci(vd_ht$HT_top_calib, vd_ht$HT_top_t2, "HT_top", "calib"))
    }
  }

  stats
}

# Compute by variant
cat("Computing statistics by variant...\n")
validation_stats <- list()

for (var in sort(unique(validation_data$VARIANT))) {
  vd <- validation_data[VARIANT == var]
  if (nrow(vd) < 10) next

  stats <- compute_variant_stats(vd, var)
  validation_stats[[var]] <- stats

  ht_str <- ""
  if (!is.null(stats$HT_top_r2_calib)) {
    ht_str <- sprintf(" | HT R2 c/d=%.3f/%.3f", stats$HT_top_r2_calib, stats$HT_top_r2_default)
  }
  cat(sprintf("  %s: n=%d | BA RMSE c/d=%.1f/%.1f | BA R2 c/d=%.3f/%.3f | VOL R2 c/d=%.3f/%.3f%s\n",
              var, nrow(vd),
              stats$BA_RMSE_calib, stats$BA_RMSE_default,
              stats$BA_r2_calib, stats$BA_r2_default,
              ifelse(is.null(stats$VOL_CFGRS_r2_calib), NA, stats$VOL_CFGRS_r2_calib),
              ifelse(is.null(stats$VOL_CFGRS_r2_default), NA, stats$VOL_CFGRS_r2_default),
              ht_str))
}

# Overall
cat("\nComputing overall statistics...\n")
overall <- compute_variant_stats(validation_data, "OVERALL")
validation_stats[["OVERALL"]] <- overall

# Save results
validation_stats_dt <- rbindlist(validation_stats, fill = TRUE)
fwrite(validation_stats_dt, file.path(output_root, paste0("manuscript_tables/fia_benchmark_results", tag_suffix, ".csv")))

cat("\n")
cat(sprintf("OVERALL: BA RMSE calib=%.1f default=%.1f | R2 calib=%.3f default=%.3f\n",
            overall$BA_RMSE_calib, overall$BA_RMSE_default,
            overall$BA_r2_calib, overall$BA_r2_default))
cat(sprintf("         BA MAE calib=%.1f default=%.1f | bias%% calib=%.1f%% default=%.1f%%\n",
            overall$BA_MAE_calib, overall$BA_MAE_default,
            overall$BA_bias_pct_calib, overall$BA_bias_pct_default))
cat(sprintf("         BA equiv calib=%.1f%% default=%.1f%% | d calib=%.3f default=%.3f\n",
            overall$BA_equiv_calib, overall$BA_equiv_default,
            overall$BA_d_calib, overall$BA_d_default))
if (!is.null(overall$VOL_CFGRS_r2_calib)) {
  cat(sprintf("         VOL_CFGRS RMSE calib=%.1f default=%.1f | R2 calib=%.3f default=%.3f\n",
              overall$VOL_CFGRS_RMSE_calib, overall$VOL_CFGRS_RMSE_default,
              overall$VOL_CFGRS_r2_calib, overall$VOL_CFGRS_r2_default))
  cat(sprintf("         VOL_CFGRS bias%% calib=%.1f%% default=%.1f%% | d calib=%.3f default=%.3f\n",
              overall$VOL_CFGRS_bias_pct_calib, overall$VOL_CFGRS_bias_pct_default,
              overall$VOL_CFGRS_d_calib, overall$VOL_CFGRS_d_default))
}
if (!is.null(overall$VOL_CFNET_r2_calib)) {
  cat(sprintf("         VOL_CFNET R2 calib=%.3f default=%.3f | d calib=%.3f default=%.3f\n",
              overall$VOL_CFNET_r2_calib, overall$VOL_CFNET_r2_default,
              overall$VOL_CFNET_d_calib, overall$VOL_CFNET_d_default))
}
if (!is.null(overall$VOL_BFNET_r2_calib)) {
  cat(sprintf("         VOL_BFNET R2 calib=%.3f default=%.3f | d calib=%.3f default=%.3f\n",
              overall$VOL_BFNET_r2_calib, overall$VOL_BFNET_r2_default,
              overall$VOL_BFNET_d_calib, overall$VOL_BFNET_d_default))
}
if (!is.null(overall$PAI_CFGRS_r2_calib)) {
  cat(sprintf("         PAI_CFGRS obs=%.1f calib=%.1f default=%.1f | R2 c/d=%.3f/%.3f\n",
              overall$PAI_CFGRS_obs_mean_calib, overall$PAI_CFGRS_pred_mean_calib,
              overall$PAI_CFGRS_pred_mean_default,
              overall$PAI_CFGRS_r2_calib, overall$PAI_CFGRS_r2_default))
}
# Gross vs Net PAI summary
if (!is.null(overall$PAI_net_obs_mean)) {
  cat(sprintf("         PAI_net  obs=%.1f calib=%.1f default=%.1f | R2 c/d=%.3f/%.3f\n",
              overall$PAI_net_obs_mean, overall$PAI_net_pred_mean_calib,
              overall$PAI_net_pred_mean_default,
              overall$PAI_net_r2_calib, overall$PAI_net_r2_default))
  cat(sprintf("         PAI_gross      calib=%.1f default=%.1f (no obs equivalent)\n",
              overall$PAI_gross_pred_mean_calib, overall$PAI_gross_pred_mean_default))
}
# Top height summary
if (!is.null(overall$HT_top_r2_calib)) {
  cat(sprintf("         HT_top RMSE c/d=%.1f/%.1f | R2 c/d=%.3f/%.3f | bias c/d=%.1f/%.1f ft\n",
              overall$HT_top_RMSE_calib, overall$HT_top_RMSE_default,
              overall$HT_top_r2_calib, overall$HT_top_r2_default,
              overall$HT_top_bias_calib, overall$HT_top_bias_default))
}

# ==============================================================================
# STEP 6: Generate Figures
# ==============================================================================

cat("\nSTEP 6: Generate Figures\n")
cat(strrep("-", 80), "\n\n")

fig_dir <- file.path(output_root, "manuscript_figures")

# Use complete cases for plotting
plot_data <- validation_data[!is.na(BA_pred_calib) & !is.na(BA_pred_default) &
                              !is.na(BA_t2) & BA_t2 > 0]

# --- Figure 1: Predicted vs Observed BA (Calibrated) with 80% CI ribbons ---
cat("Figure 1: BA predicted vs observed (calibrated) with CI ribbons...\n")
has_ba_ci <- all(c("BA_pred_calib_lo", "BA_pred_calib_hi") %in% names(plot_data))
p1 <- ggplot(plot_data, aes(x = BA_t2, y = BA_pred_calib))
if (has_ba_ci && CI_BRACKETS_ENABLED) {
  p1 <- p1 + geom_linerange(
    aes(ymin = BA_pred_calib_lo, ymax = BA_pred_calib_hi),
    alpha = 0.05, color = "steelblue", linewidth = 0.3,
    data = plot_data[!is.na(BA_pred_calib_lo) & !is.na(BA_pred_calib_hi)])
}
p1 <- p1 +
  geom_point(alpha = 0.15, size = 0.5, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~VARIANT, scales = "free") +
  labs(title = "Basal Area: Observed vs Calibrated Prediction (80% CI brackets)",
       x = expression("Observed BA (ft"^2~"ac"^{-1}*")"),
       y = expression("Predicted BA (ft"^2~"ac"^{-1}*")")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "01_ba_pred_vs_obs_calibrated.png"), p1,
       width = 14, height = 10, dpi = 300)

# --- Figure 2: Predicted vs Observed BA (Default) ---
cat("Figure 2: BA predicted vs observed (default)...\n")
p2 <- ggplot(plot_data, aes(x = BA_t2, y = BA_pred_default)) +
  geom_point(alpha = 0.15, size = 0.5, color = "darkorange") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~VARIANT, scales = "free") +
  labs(title = "Basal Area: Observed vs Default Prediction",
       x = expression("Observed BA (ft"^2~"ac"^{-1}*")"),
       y = expression("Predicted BA (ft"^2~"ac"^{-1}*")")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "02_ba_pred_vs_obs_default.png"), p2,
       width = 14, height = 10, dpi = 300)

# --- Figure 3: RMSE comparison bar chart ---
cat("Figure 3: RMSE comparison...\n")
rmse_long <- rbindlist(list(
  validation_stats_dt[VARIANT != "OVERALL",
    .(VARIANT, metric = "BA", Calibrated = BA_RMSE_calib, Default = BA_RMSE_default)],
  validation_stats_dt[VARIANT != "OVERALL",
    .(VARIANT, metric = "QMD", Calibrated = QMD_RMSE_calib, Default = QMD_RMSE_default)]
))
rmse_melt <- melt(rmse_long, id.vars = c("VARIANT", "metric"),
                   variable.name = "approach", value.name = "RMSE")

p3 <- ggplot(rmse_melt[metric == "BA"],
             aes(x = reorder(VARIANT, -RMSE), y = RMSE, fill = approach)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = c(Calibrated = "steelblue", Default = "darkorange")) +
  labs(title = "Basal Area RMSE: Calibrated vs Default",
       x = "Variant", y = expression("RMSE (ft"^2~"ac"^{-1}*")"), fill = "Approach") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "03_ba_rmse_comparison.png"), p3,
       width = 12, height = 6, dpi = 300)

# --- Figure 4: Bias dumbbell ---
cat("Figure 4: BA bias dumbbell...\n")
bias_data <- validation_stats_dt[VARIANT != "OVERALL",
  .(VARIANT, Calibrated = BA_bias_calib, Default = BA_bias_default)]
bias_melt <- melt(bias_data, id.vars = "VARIANT",
                   variable.name = "approach", value.name = "bias")

p4 <- ggplot(bias_melt, aes(x = reorder(VARIANT, bias), y = bias, color = approach)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(Calibrated = "steelblue", Default = "darkorange")) +
  coord_flip() +
  labs(title = "Basal Area Bias by Variant",
       x = "Variant", y = expression("Bias (ft"^2~"ac"^{-1}*")"), color = "Approach") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "04_ba_bias_dumbbell.png"), p4,
       width = 10, height = 8, dpi = 300)

# --- Figure 5: R-squared comparison ---
cat("Figure 5: R-squared comparison...\n")
r2_data <- validation_stats_dt[VARIANT != "OVERALL",
  .(VARIANT, Calibrated = BA_r2_calib, Default = BA_r2_default)]

p5 <- ggplot(r2_data, aes(x = Default, y = Calibrated, label = VARIANT)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "steelblue") +
  geom_text(size = 2.5, nudge_y = 0.01, check_overlap = TRUE) +
  labs(title = expression("BA R"^2*": Calibrated vs Default"),
       x = expression("Default R"^2), y = expression("Calibrated R"^2)) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "05_ba_r2_scatter.png"), p5,
       width = 8, height = 7, dpi = 300)

# --- Figure 6: Equivalence comparison ---
cat("Figure 6: Equivalence comparison...\n")
equiv_data <- validation_stats_dt[VARIANT != "OVERALL",
  .(VARIANT, Calibrated = BA_equiv_calib, Default = BA_equiv_default)]
equiv_melt <- melt(equiv_data, id.vars = "VARIANT",
                    variable.name = "approach", value.name = "equiv_pct")

p6 <- ggplot(equiv_melt, aes(x = reorder(VARIANT, -equiv_pct), y = equiv_pct, fill = approach)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = c(Calibrated = "steelblue", Default = "darkorange")) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  labs(title = "BA Equivalence (within 20% of observed)",
       x = "Variant", y = "% within 20%", fill = "Approach") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "06_ba_equivalence.png"), p6,
       width = 12, height = 6, dpi = 300)

# --- Figure 7: QMD predicted vs observed ---
cat("Figure 7: QMD predicted vs observed...\n")
p7 <- ggplot(plot_data[!is.na(QMD_t2) & !is.na(QMD_pred_calib)],
             aes(x = QMD_t2, y = QMD_pred_calib)) +
  geom_point(alpha = 0.15, size = 0.5, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~VARIANT, scales = "free") +
  labs(title = "QMD: Observed vs Calibrated Prediction",
       x = "Observed QMD (in.)", y = "Predicted QMD (in.)") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "07_qmd_pred_vs_obs_calibrated.png"), p7,
       width = 14, height = 10, dpi = 300)

# --- Figure 8: Volume predicted vs observed (calibrated) with 80% CI ---
cat("Figure 8: Volume CFGRS predicted vs observed with CI...\n")
vol_plot <- plot_data[VOL_CFGRS_t2 > 0 & !is.na(VOL_CFGRS_pred_calib)]
if (nrow(vol_plot) > 0) {
  has_vol_ci <- all(c("VOL_CFGRS_pred_calib_lo", "VOL_CFGRS_pred_calib_hi") %in% names(vol_plot))
  p8 <- ggplot(vol_plot, aes(x = VOL_CFGRS_t2, y = VOL_CFGRS_pred_calib))
  if (has_vol_ci && CI_BRACKETS_ENABLED) {
    p8 <- p8 + geom_linerange(
      aes(ymin = VOL_CFGRS_pred_calib_lo, ymax = VOL_CFGRS_pred_calib_hi),
      alpha = 0.05, color = "steelblue", linewidth = 0.3,
      data = vol_plot[!is.na(VOL_CFGRS_pred_calib_lo) & !is.na(VOL_CFGRS_pred_calib_hi)])
  }
  p8 <- p8 +
    geom_point(alpha = 0.15, size = 0.5, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    facet_wrap(~VARIANT, scales = "free") +
    labs(title = "Gross Cubic Foot Volume: Observed vs Calibrated (80% CI brackets)",
         x = expression("Observed Volume (ft"^3~"ac"^{-1}*")"),
         y = expression("Predicted Volume (ft"^3~"ac"^{-1}*")")) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "08_vol_cfgrs_pred_vs_obs_calibrated.png"), p8,
         width = 14, height = 10, dpi = 300)
}

# --- Figure 9: Volume RMSE comparison ---
cat("Figure 9: Volume RMSE comparison...\n")
vol_rmse <- validation_stats_dt[VARIANT != "OVERALL" & !is.na(VOL_CFGRS_RMSE_calib),
  .(VARIANT, Calibrated = VOL_CFGRS_RMSE_calib, Default = VOL_CFGRS_RMSE_default)]
if (nrow(vol_rmse) > 0) {
  vol_rmse_melt <- melt(vol_rmse, id.vars = "VARIANT",
                         variable.name = "approach", value.name = "RMSE")
  p9 <- ggplot(vol_rmse_melt, aes(x = reorder(VARIANT, -RMSE), y = RMSE, fill = approach)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    scale_fill_manual(values = c(Calibrated = "steelblue", Default = "darkorange")) +
    labs(title = "Volume (CFGRS) RMSE: Calibrated vs Default",
         x = "Variant", y = expression("RMSE (ft"^3~"ac"^{-1}*")"), fill = "Approach") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "09_vol_rmse_comparison.png"), p9,
         width = 12, height = 6, dpi = 300)
}

# --- Figure 10: PAI predicted vs observed ---
cat("Figure 10: PAI predicted vs observed...\n")
pai_plot <- plot_data[!is.na(PAI_CFGRS_obs) & !is.na(PAI_CFGRS_pred_calib) &
                      VOL_CFGRS_t1 > 0 & VOL_CFGRS_t2 > 0]
if (nrow(pai_plot) > 0) {
  pai_long <- rbindlist(list(
    pai_plot[, .(VARIANT, PAI_obs = PAI_CFGRS_obs, PAI_pred = PAI_CFGRS_pred_calib,
                 approach = "Calibrated")],
    pai_plot[, .(VARIANT, PAI_obs = PAI_CFGRS_obs, PAI_pred = PAI_CFGRS_pred_default,
                 approach = "Default")]
  ))
  p10 <- ggplot(pai_long, aes(x = PAI_obs, y = PAI_pred, color = approach)) +
    geom_point(alpha = 0.10, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c(Calibrated = "steelblue", Default = "darkorange")) +
    facet_wrap(~VARIANT, scales = "free") +
    labs(title = "Periodic Annual Increment (Gross CF Volume): Observed vs Predicted",
         x = expression("Observed PAI (ft"^3~"ac"^{-1}~"yr"^{-1}*")"),
         y = expression("Predicted PAI (ft"^3~"ac"^{-1}~"yr"^{-1}*")"),
         color = "Approach") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "10_pai_pred_vs_obs.png"), p10,
         width = 14, height = 10, dpi = 300)
}

# --- Figure 11: Top height predicted vs observed with CI ---
cat("Figure 11: Top height predicted vs observed with CI...\n")
ht_plot <- plot_data[!is.na(HT_top_t2) & HT_top_t2 > 0 &
                      !is.na(HT_top_calib) & !is.na(HT_top_default)]
if (nrow(ht_plot) > 0) {
  has_ht_ci <- all(c("HT_top_calib_lo", "HT_top_calib_hi") %in% names(ht_plot))
  ht_long <- rbindlist(list(
    ht_plot[, .(VARIANT, HT_obs = HT_top_t2, HT_pred = HT_top_calib,
                 approach = "Calibrated")],
    ht_plot[, .(VARIANT, HT_obs = HT_top_t2, HT_pred = HT_top_default,
                 approach = "Default")]
  ))
  p11_ht <- ggplot(ht_long, aes(x = HT_obs, y = HT_pred, color = approach))
  if (has_ht_ci && CI_BRACKETS_ENABLED) {
    p11_ht <- p11_ht + geom_linerange(
      aes(x = HT_top_t2, ymin = HT_top_calib_lo, ymax = HT_top_calib_hi),
      alpha = 0.05, color = "steelblue", linewidth = 0.3, inherit.aes = FALSE,
      data = ht_plot[!is.na(HT_top_calib_lo) & !is.na(HT_top_calib_hi)])
  }
  p11_ht <- p11_ht +
    geom_point(alpha = 0.10, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c(Calibrated = "steelblue", Default = "darkorange")) +
    facet_wrap(~VARIANT, scales = "free") +
    labs(title = "Top Height (Lord, 40 TPA): Observed vs Predicted (80% CI brackets)",
         x = "Observed Top Height (ft)", y = "Predicted Top Height (ft)",
         color = "Approach") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "11_ht_top_pred_vs_obs.png"), p11_ht,
         width = 14, height = 10, dpi = 300)
}

# --- Figure 12: Gross vs Net PAI comparison ---
cat("Figure 12: Gross vs Net PAI comparison...\n")
pai_gn_plot <- plot_data[!is.na(PAI_net_obs) & !is.na(PAI_gross_pred_calib) &
                          VOL_CFGRS_t1 > 0 & VOL_CFGRS_t2 > 0]
if (nrow(pai_gn_plot) > 0) {
  pai_gn_long <- rbindlist(list(
    pai_gn_plot[, .(VARIANT, PAI_obs = PAI_net_obs, PAI_pred = PAI_net_pred_calib,
                     type = "Net (Calibrated)")],
    pai_gn_plot[, .(VARIANT, PAI_obs = PAI_net_obs, PAI_pred = PAI_gross_pred_calib,
                     type = "Gross (Calibrated)")],
    pai_gn_plot[, .(VARIANT, PAI_obs = PAI_net_obs, PAI_pred = PAI_net_pred_default,
                     type = "Net (Default)")]
  ))
  p12 <- ggplot(pai_gn_long, aes(x = PAI_obs, y = PAI_pred, color = type)) +
    geom_point(alpha = 0.08, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c("Net (Calibrated)" = "steelblue",
                                   "Gross (Calibrated)" = "forestgreen",
                                   "Net (Default)" = "darkorange")) +
    facet_wrap(~VARIANT, scales = "free") +
    labs(title = "PAI Decomposition: Gross vs Net (CF Volume)",
         x = expression("Observed Net PAI (ft"^3~"ac"^{-1}~"yr"^{-1}*")"),
         y = expression("Predicted PAI (ft"^3~"ac"^{-1}~"yr"^{-1}*")"),
         color = "Type") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "12_pai_gross_vs_net.png"), p12,
         width = 14, height = 10, dpi = 300)
}

# --- Figure 13: Comprehensive multi-metric R2 heatmap ---
cat("Figure 13: Multi-metric R2 heatmap...\n")
r2_cols <- c("BA_r2_calib", "BA_r2_default", "TPA_r2_calib", "TPA_r2_default",
             "QMD_r2_calib", "QMD_r2_default", "SDI_r2_calib", "SDI_r2_default",
             "VOL_CFGRS_r2_calib", "VOL_CFGRS_r2_default",
             "VOL_CFNET_r2_calib", "VOL_CFNET_r2_default",
             "VOL_BFNET_r2_calib", "VOL_BFNET_r2_default",
             "HT_top_r2_calib", "HT_top_r2_default")
r2_present <- r2_cols[r2_cols %in% names(validation_stats_dt)]
if (length(r2_present) > 0) {
  r2_wide <- validation_stats_dt[VARIANT != "OVERALL", c("VARIANT", r2_present), with = FALSE]
  r2_melt <- melt(r2_wide, id.vars = "VARIANT", variable.name = "metric_approach",
                   value.name = "R2")
  r2_melt[, metric := gsub("_r2_(calib|default)", "", metric_approach)]
  r2_melt[, approach := fifelse(grepl("calib", metric_approach), "Calibrated", "Default")]

  p13 <- ggplot(r2_melt, aes(x = metric, y = reorder(VARIANT, R2), fill = R2)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", mid = "white", high = "steelblue",
                          midpoint = 0.5, limits = c(0, 1)) +
    facet_wrap(~approach) +
    labs(title = expression("R"^2~"by Variant and Metric"),
         x = "Metric", y = "Variant", fill = expression("R"^2)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "13_r2_heatmap.png"), p13,
         width = 12, height = 10, dpi = 300)
}

# --- Figure 14: Bootstrap CI on BA and VOL RMSE by variant ---
if (BOOTSTRAP_CI_ENABLED &&
    all(c("BA_RMSE_lo_calib", "BA_RMSE_hi_calib") %in% names(validation_stats_dt))) {
  cat("Figure 14: Bootstrap CI on calibrated RMSE...\n")
  ci_plot <- validation_stats_dt[VARIANT != "OVERALL" & !is.na(BA_RMSE_lo_calib)]

  if (nrow(ci_plot) > 0) {
    p14 <- ggplot(ci_plot, aes(x = reorder(VARIANT, -BA_RMSE_calib), y = BA_RMSE_calib)) +
      geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
      geom_errorbar(aes(ymin = BA_RMSE_lo_calib, ymax = BA_RMSE_hi_calib),
                    width = 0.25, color = "gray30", linewidth = 0.5) +
      geom_point(aes(y = BA_RMSE_default), shape = 4, size = 2.5, color = "darkorange",
                 stroke = 1) +
      labs(title = "Calibrated BA RMSE with 95% Bootstrap CI (x = default)",
           x = "Variant",
           y = expression("RMSE (ft"^2~"ac"^{-1}*")")) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(fig_dir, "14_ba_rmse_bootstrap_ci.png"), p14,
           width = 12, height = 6, dpi = 300)
  }

  # Figure 14b: Bootstrap CI on R2
  if ("BA_r2_lo_calib" %in% names(validation_stats_dt)) {
    ci_r2 <- validation_stats_dt[VARIANT != "OVERALL" & !is.na(BA_r2_lo_calib)]
    if (nrow(ci_r2) > 0) {
      p14b <- ggplot(ci_r2, aes(x = reorder(VARIANT, BA_r2_calib), y = BA_r2_calib)) +
        geom_point(size = 3, color = "steelblue") +
        geom_errorbar(aes(ymin = BA_r2_lo_calib, ymax = BA_r2_hi_calib),
                      width = 0.25, color = "steelblue", linewidth = 0.5) +
        geom_point(aes(y = BA_r2_default), shape = 4, size = 2.5, color = "darkorange",
                   stroke = 1) +
        geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50") +
        coord_flip() +
        labs(title = expression("Calibrated BA R"^2~"with 95% Bootstrap CI (x = default)"),
             x = "Variant", y = expression("R"^2)) +
        theme_minimal(base_size = 11)
      ggsave(file.path(fig_dir, "14b_ba_r2_bootstrap_ci.png"), p14b,
             width = 10, height = 8, dpi = 300)
    }
  }
}

# --- Figure 15: CI bracket coverage analysis ---
if (CI_BRACKETS_ENABLED &&
    all(c("BA_pred_calib_lo", "BA_pred_calib_hi") %in% names(validation_data))) {
  cat("Figure 15: CI bracket coverage...\n")
  cov_data <- validation_data[!is.na(BA_pred_calib_lo) & !is.na(BA_pred_calib_hi)]
  if (nrow(cov_data) > 0) {
    cov_data[, BA_covered := BA_t2 >= BA_pred_calib_lo & BA_t2 <= BA_pred_calib_hi]
    cov_by_var <- cov_data[, .(coverage = 100 * mean(BA_covered, na.rm = TRUE),
                                 n = .N), by = VARIANT]
    p15 <- ggplot(cov_by_var, aes(x = reorder(VARIANT, -coverage), y = coverage)) +
      geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
      geom_text(aes(label = sprintf("%.0f%%", coverage)), vjust = -0.3, size = 3) +
      labs(title = "80% CI Bracket Coverage: Observed BA within Predicted Interval",
           x = "Variant", y = "Coverage (%)") +
      ylim(0, 105) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(fig_dir, "15_ci_bracket_coverage.png"), p15,
           width = 12, height = 6, dpi = 300)

    # Print coverage summary
    overall_cov <- cov_data[, 100 * mean(BA_covered, na.rm = TRUE)]
    cat(sprintf("  Overall BA 80%% CI coverage: %.1f%% (n=%d)\n",
                overall_cov, nrow(cov_data)))
  }
}

cat("All figures saved to:", fig_dir, "\n\n")

# ==============================================================================
# STEP 7: Summary Report
# ==============================================================================

cat("STEP 7: Summary Report\n")
cat(strrep("=", 80), "\n\n")

# Count improvements
n_improved_rmse <- sum(validation_stats_dt[VARIANT != "OVERALL"]$BA_RMSE_calib <
                       validation_stats_dt[VARIANT != "OVERALL"]$BA_RMSE_default, na.rm = TRUE)
n_total <- nrow(validation_stats_dt[VARIANT != "OVERALL"])

cat(sprintf("BA RMSE improved: %d / %d variants\n", n_improved_rmse, n_total))
cat(sprintf("Overall BA RMSE reduction: %.1f%%\n",
            100 * (1 - overall$BA_RMSE_calib / overall$BA_RMSE_default)))
cat(sprintf("Overall BA R2 improvement: %.3f -> %.3f\n",
            overall$BA_r2_default, overall$BA_r2_calib))
cat(sprintf("Overall BA equivalence: %.1f%% -> %.1f%%\n",
            overall$BA_equiv_default, overall$BA_equiv_calib))
# Bootstrap CI summary
if (BOOTSTRAP_CI_ENABLED && !is.null(overall$BA_RMSE_lo_calib)) {
  cat(sprintf("  BA RMSE 95%% CI: [%.1f, %.1f] | R2 95%% CI: [%.3f, %.3f]\n",
              overall$BA_RMSE_lo_calib, overall$BA_RMSE_hi_calib,
              overall$BA_r2_lo_calib, overall$BA_r2_hi_calib))
  cat(sprintf("  BA bias 95%% CI: [%.1f, %.1f]\n",
              overall$BA_bias_lo_calib, overall$BA_bias_hi_calib))
}
if (!is.null(overall$VOL_CFGRS_RMSE_lo_calib)) {
  cat(sprintf("  VOL_CFGRS RMSE 95%% CI: [%.1f, %.1f] | R2 95%% CI: [%.3f, %.3f]\n",
              overall$VOL_CFGRS_RMSE_lo_calib, overall$VOL_CFGRS_RMSE_hi_calib,
              overall$VOL_CFGRS_r2_lo_calib, overall$VOL_CFGRS_r2_hi_calib))
}
if (!is.null(overall$HT_top_RMSE_lo_calib)) {
  cat(sprintf("  HT_top RMSE 95%% CI: [%.1f, %.1f] | R2 95%% CI: [%.3f, %.3f]\n",
              overall$HT_top_RMSE_lo_calib, overall$HT_top_RMSE_hi_calib,
              overall$HT_top_r2_lo_calib, overall$HT_top_r2_hi_calib))
}

# Volume summary (all types)
for (vol_type in c("VOL_CFGRS", "VOL_CFNET", "VOL_BFNET")) {
  rmse_c <- overall[[paste0(vol_type, "_RMSE_calib")]]
  rmse_d <- overall[[paste0(vol_type, "_RMSE_default")]]
  r2_c   <- overall[[paste0(vol_type, "_r2_calib")]]
  r2_d   <- overall[[paste0(vol_type, "_r2_default")]]
  d_c    <- overall[[paste0(vol_type, "_d_calib")]]
  d_d    <- overall[[paste0(vol_type, "_d_default")]]
  if (!is.null(r2_c)) {
    rmse_col_c <- paste0(vol_type, "_RMSE_calib")
    rmse_col_d <- paste0(vol_type, "_RMSE_default")
    vs_dt <- validation_stats_dt[VARIANT != "OVERALL" & !is.na(get(rmse_col_c))]
    n_improved_vol <- sum(vs_dt[[rmse_col_c]] < vs_dt[[rmse_col_d]], na.rm = TRUE)
    n_vol_total <- nrow(vs_dt)
    cat(sprintf("\n%s RMSE improved: %d / %d variants\n", vol_type, n_improved_vol, n_vol_total))
    cat(sprintf("Overall %s RMSE reduction: %.1f%% | R2: %.3f -> %.3f | d: %.3f -> %.3f\n",
                vol_type, 100 * (1 - rmse_c / rmse_d), r2_d, r2_c, d_d, d_c))
  }
}

# PAI summary (legacy)
if (!is.null(overall$PAI_CFGRS_r2_calib)) {
  cat(sprintf("\nPAI_CFGRS obs=%.1f, calib=%.1f, default=%.1f | R2 c/d=%.3f/%.3f\n",
              overall$PAI_CFGRS_obs_mean_calib, overall$PAI_CFGRS_pred_mean_calib,
              overall$PAI_CFGRS_pred_mean_default,
              overall$PAI_CFGRS_r2_calib, overall$PAI_CFGRS_r2_default))
}

# Gross vs Net PAI summary
if (!is.null(overall$PAI_net_obs_mean)) {
  cat(sprintf("\nPAI decomposition (CFGRS):\n"))
  cat(sprintf("  Net PAI:   obs=%.1f calib=%.1f default=%.1f | R2 c/d=%.3f/%.3f\n",
              overall$PAI_net_obs_mean, overall$PAI_net_pred_mean_calib,
              overall$PAI_net_pred_mean_default,
              overall$PAI_net_r2_calib, overall$PAI_net_r2_default))
  cat(sprintf("  Gross PAI: calib=%.1f default=%.1f (no observed counterpart)\n",
              overall$PAI_gross_pred_mean_calib, overall$PAI_gross_pred_mean_default))
}

# Top height summary
if (!is.null(overall$HT_top_r2_calib)) {
  cat(sprintf("\nTop Height (Lord 40 TPA):\n"))
  cat(sprintf("  RMSE: calib=%.1f default=%.1f ft | R2: %.3f -> %.3f\n",
              overall$HT_top_RMSE_calib, overall$HT_top_RMSE_default,
              overall$HT_top_r2_default, overall$HT_top_r2_calib))
  cat(sprintf("  Bias: calib=%.1f default=%.1f ft | d: %.3f -> %.3f\n",
              overall$HT_top_bias_calib, overall$HT_top_bias_default,
              overall$HT_top_d_default, overall$HT_top_d_calib))
  if ("HT_top_RMSE_calib" %in% names(validation_stats_dt)) {
    vs_ht <- validation_stats_dt[VARIANT != "OVERALL" & !is.na(HT_top_RMSE_calib)]
    n_ht_improved <- sum(vs_ht$HT_top_RMSE_calib < vs_ht$HT_top_RMSE_default, na.rm = TRUE)
    cat(sprintf("  HT_top RMSE improved: %d / %d variants\n", n_ht_improved, nrow(vs_ht)))
  }
}

cat("\nTop 5 variants by calibrated BA R2:\n")
top5 <- validation_stats_dt[VARIANT != "OVERALL"][order(-BA_r2_calib)][1:min(5, .N)]
for (i in 1:nrow(top5)) {
  vol_str <- ""
  if (!is.na(top5$VOL_CFGRS_r2_calib[i])) {
    vol_str <- sprintf(", VOL R2=%.3f (def=%.3f)", top5$VOL_CFGRS_r2_calib[i], top5$VOL_CFGRS_r2_default[i])
  }
  cat(sprintf("  %s: BA R2=%.3f (def=%.3f), RMSE=%.1f (def=%.1f)%s\n",
              top5$VARIANT[i], top5$BA_r2_calib[i], top5$BA_r2_default[i],
              top5$BA_RMSE_calib[i], top5$BA_RMSE_default[i], vol_str))
}

cat("\nBottom 5 variants by calibrated BA R2:\n")
bot5 <- validation_stats_dt[VARIANT != "OVERALL"][order(BA_r2_calib)][1:min(5, .N)]
for (i in 1:nrow(bot5)) {
  cat(sprintf("  %s: BA R2=%.3f (def=%.3f), RMSE=%.1f (def=%.1f)\n",
              bot5$VARIANT[i], bot5$BA_r2_calib[i], bot5$BA_r2_default[i],
              bot5$BA_RMSE_calib[i], bot5$BA_RMSE_default[i]))
}

cat("\n")
cat(strrep("=", 80), "\n")
cat("FIA Benchmark Engine v2: Complete\n")
cat(strrep("=", 80), "\n")
