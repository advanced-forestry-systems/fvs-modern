#!/usr/bin/env Rscript
# =============================================================================
# FVS Calibration: Compare Default vs Calibrated Wykoff Diameter Growth
#
# Produces:
#   1. Prediction accuracy comparison on held out FIA data
#   2. Growth curve comparison across DBH range for key species
#   3. Residual diagnostics for both models
#   4. Summary table of RMSE, bias, R^2 by species
#
# Usage: Rscript calibration/R/07_compare_default_vs_calibrated.R --variant ca
#        Rscript calibration/R/07_compare_default_vs_calibrated.R --all
# =============================================================================

library(tidyverse)
library(jsonlite)
library(data.table)
library(logger)

# =============================================================================
# Parse Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
run_all <- FALSE
variants_to_run <- character(0)

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variants_to_run <- tolower(args[i + 1])
    }
    if (args[i] == "--all") run_all <- TRUE
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
config_dir <- file.path(project_root, "config")
output_base <- file.path(calibration_dir, "output", "variants")
comparison_dir <- file.path(calibration_dir, "output", "comparisons")
dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(calibration_dir, "logs", "07_comparison.log")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)

# Auto detect completed variants if --all
if (run_all || length(variants_to_run) == 0) {
  all_variants <- list.dirs(output_base, recursive = FALSE, full.names = FALSE)
  variants_to_run <- all_variants[sapply(all_variants, function(v) {
    file.exists(file.path(output_base, v, "diameter_growth_posterior.csv")) ||
    file.exists(file.path(output_base, v, "diameter_growth_map.csv"))
  })]
  logger::log_info("Auto detected {length(variants_to_run)} completed variants: {paste(variants_to_run, collapse = ', ')}")
}

if (length(variants_to_run) == 0) {
  stop("No completed variants found. Run calibration first.")
}

# =============================================================================
# Core Functions
# =============================================================================

#' Detect which Wykoff equation form a variant uses
detect_equation_type <- function(config) {
  other <- config$categories$other
  growth <- config$categories$growth

  if ("DGLD" %in% names(other)) return("western")       # CA, BM, CI, CR, EC, EM, IE, SO, UT, WC, WS, PN
  if ("INTERC" %in% names(other)) return("central")      # CS, LS
  if (!is.null(growth$B1) && length(unlist(growth$B1)) > 0) return("ne")  # NE, ACD
  return("unknown")
}


#' Apply the default FVS Wykoff equation on the original (unstandardized) scale
#'
#' Supports three equation forms:
#'   "western" (CA type): uses DGLD, DGDS, DGSITE, DGSLOP, DGEL, etc.
#'   "central" (CS/LS type): uses INTERC, VDBHC, DBHC, SBAC, BALC, etc.
#'   "ne" (NE type): uses B1, B2 potential basal area growth model
predict_default_wykoff <- function(data, config) {

  eq_type <- detect_equation_type(config)
  other <- config$categories$other
  growth <- config$categories$growth
  sp_idx <- data$fvs_sp_idx

  pad_to <- function(x, n) {
    if (is.null(x) || length(x) == 0) return(rep(0, n))
    if (length(x) < n) c(x, rep(0, n - length(x))) else x[1:n]
  }

  if (eq_type == "western") {
    # ------------------------------------------------------------------
    # Western variant (CA, BM, CI, CR, etc.)
    # ln(DDS) = ln(WEIBB1[sp]) + DGLD[sp]*ln(DBH) + DGDS[sp]*DBH^2 +
    #           DGSITE[sp]*ln(SI) + DGSLOP[sp]*SLOPE + DGSLSQ[sp]*SLOPE^2 +
    #           DGSASP[sp]*SLOPE*sin(ASP) + DGCASP[sp]*SLOPE*cos(ASP) +
    #           DGEL[sp]*ELEV_100ft + DGELSQ[sp]*ELEV_100ft^2 +
    #           DGCR[sp]*CR + DGCRSQ[sp]*CR^2 +
    #           DGBAL[sp]*BAL + DGPCCF[sp]*BA
    #
    # Units: ELEV in 100s of feet, SLOPE as decimal (0 to 1), CR as decimal
    # ------------------------------------------------------------------
    weibb1 <- unlist(growth$WEIBB1)
    n_sp <- max(length(weibb1), length(unlist(other$DGLD)), max(sp_idx, na.rm = TRUE))
    if (is.null(weibb1) || length(weibb1) == 0) weibb1 <- rep(1, n_sp)
    if (length(weibb1) < n_sp) weibb1 <- c(weibb1, rep(1, n_sp - length(weibb1)))

    dgld   <- pad_to(unlist(other$DGLD), n_sp)
    dgds   <- pad_to(unlist(other$DGDS), n_sp)
    dgsite <- pad_to(unlist(other$DGSITE), n_sp)
    dgslop <- pad_to(unlist(other$DGSLOP), n_sp)
    dgslsq <- pad_to(unlist(other$DGSLSQ), n_sp)
    dgsasp <- pad_to(unlist(other$DGSASP), n_sp)
    dgcasp <- pad_to(unlist(other$DGCASP), n_sp)
    dgel   <- pad_to(unlist(other$DGEL), n_sp)
    dgelsq <- pad_to(unlist(other$DGELSQ), n_sp)
    dgcr   <- pad_to(unlist(other$DGCR), n_sp)
    dgcrsq <- pad_to(unlist(other$DGCRSQ), n_sp)
    dgbal  <- pad_to(unlist(other$DGBAL), n_sp)
    dgpccf <- pad_to(unlist(other$DGPCCF), n_sp)

    # Convert FIA data units to FVS units
    elev_100ft <- data$ELEV / 100         # feet to 100s of feet
    slope_dec  <- data$SLOPE / 100        # percent to decimal

    # CR: default FVS expects 0 to 1 scale
    # CR_pct in data = FIA_CR * 100, so CR_pct/10000 gives 0 to 1
    # Use CR_01 if available, otherwise compute
    cr_dec <- if ("CR_01" %in% names(data)) data$CR_01 else pmin(pmax(data$CR_pct / 10000, 0), 1)

    # Aspect interaction: data has SLOPE_SASP = SLOPE_pct * sin(ASP)
    # FVS wants slope_decimal * sin(ASP), so divide by 100
    slope_sasp_dec <- data$SLOPE_SASP / 100
    slope_casp_dec <- data$SLOPE_CASP / 100

    ln_dds_pred <- log(pmax(weibb1[sp_idx], 0.001)) +
      dgld[sp_idx]   * data$ln_DBH +
      dgds[sp_idx]   * data$DBH_sq +
      dgsite[sp_idx] * data$ln_SI +
      dgslop[sp_idx] * slope_dec +
      dgslsq[sp_idx] * slope_dec^2 +
      dgsasp[sp_idx] * slope_sasp_dec +
      dgcasp[sp_idx] * slope_casp_dec +
      dgel[sp_idx]   * elev_100ft +
      dgelsq[sp_idx] * elev_100ft^2 +
      dgcr[sp_idx]   * cr_dec +
      dgcrsq[sp_idx] * cr_dec^2 +
      dgbal[sp_idx]  * data$BAL +
      dgpccf[sp_idx] * data$BA

    return(ln_dds_pred)

  } else if (eq_type == "central") {
    # ------------------------------------------------------------------
    # Central/Lake States variant (CS, LS)
    # ln(DDS) = INTERC[sp] + VDBHC[sp]*(1/DBH) + DBHC[sp]*DBH +
    #           DBH2C[sp]*DBH^2 + RDBHC[sp]*RELDBH + RDBHSQC[sp]*RELDBHSQ +
    #           CRWNC[sp]*CR + CRSQC[sp]*CR^2 +
    #           SBAC[sp]*BA + BALC[sp]*BAL + SITEC[sp]*SI
    #
    # Units: CR as decimal, DBH in inches, BA in ft2/ac, BAL in ft2/ac
    # RELDBH = DBH / QMD (relative diameter, approximate with 1.0 if QMD unavailable)
    # ------------------------------------------------------------------
    interc  <- unlist(other$INTERC)
    n_sp <- length(interc)

    vdbhc   <- pad_to(unlist(other$VDBHC), n_sp)
    dbhc    <- pad_to(unlist(other$DBHC), n_sp)
    dbh2c   <- pad_to(unlist(other$DBH2C), n_sp)
    rdbhc   <- pad_to(unlist(other$RDBHC), n_sp)
    rdbhsqc <- pad_to(unlist(other$RDBHSQC), n_sp)
    crwnc   <- pad_to(unlist(other$CRWNC), n_sp)
    crsqc   <- pad_to(unlist(other$CRSQC), n_sp)
    sbac    <- pad_to(unlist(other$SBAC), n_sp)
    balc    <- pad_to(unlist(other$BALC), n_sp)
    sitec   <- pad_to(unlist(other$SITEC), n_sp)

    dbh <- data$DIA_t1
    cr_dec <- if ("CR_01" %in% names(data)) data$CR_01 else pmin(pmax(data$CR_pct / 10000, 0), 1)

    # Approximate QMD from stand level data if available
    qmd <- if ("QMD" %in% names(data)) data$QMD else sqrt(data$BA / (data$TPA * 0.005454))
    qmd <- pmax(qmd, 1, na.rm = TRUE)
    reldbh <- dbh / qmd

    ln_dds_pred <- interc[sp_idx] +
      vdbhc[sp_idx]   * (1.0 / pmax(dbh, 0.1)) +
      dbhc[sp_idx]    * dbh +
      dbh2c[sp_idx]   * dbh^2 +
      rdbhc[sp_idx]   * reldbh +
      rdbhsqc[sp_idx] * reldbh^2 +
      crwnc[sp_idx]   * cr_dec +
      crsqc[sp_idx]   * cr_dec^2 +
      sbac[sp_idx]    * data$BA +
      balc[sp_idx]    * data$BAL +
      sitec[sp_idx]   * data$SI

    return(ln_dds_pred)

  } else if (eq_type == "ne") {
    # ------------------------------------------------------------------
    # Northeastern variant (NE, ACD)
    # Potential basal area growth model:
    #   POTBAG = B1[sp] * SI * (1 - exp(-B2[sp] * DBH))
    #   POTBAG = POTBAG * 0.7  (reduction factor)
    #   DELD = POTBAG * BAGMOD (modifier from BAL competition)
    #   DDS = DELD + 0.005454*DBH^2  => back to squared diameter change
    #
    # Simplified: approximate BAGMOD from CR and BAL
    # ------------------------------------------------------------------
    b1 <- unlist(growth$B1)
    b2 <- unlist(growth$B2)
    n_sp <- length(b1)
    if (length(b2) < n_sp) b2 <- pad_to(b2, n_sp)

    dbh <- data$DIA_t1
    si <- data$SI
    cr_dec <- if ("CR_01" %in% names(data)) data$CR_01 else pmin(pmax(data$CR_pct / 10000, 0), 1)

    # Potential basal area growth
    potbag <- b1[sp_idx] * si * (1.0 - exp(-b2[sp_idx] * dbh))
    potbag <- potbag * 0.7

    # Simplified BAL modifier: scales 0 to 1 based on relative competition
    # In FVS this is from BALMOD subroutine; approximate as CR based modifier
    bagmod <- pmax(cr_dec, 0.1)

    deld <- potbag * bagmod

    # Convert: DDS = new_QBA - old_QBA where QBA = 0.005454 * DBH^2
    # DELD is the basal area increment
    # new_DBH = sqrt((old_QBA + DELD) / 0.005454)
    # DDS = (new_DBH^2 - old_DBH^2) approximately = DELD / 0.005454
    dds <- pmax(deld / 0.005454, 0.01)
    ln_dds_pred <- log(dds)

    return(ln_dds_pred)

  } else {
    warning("Unknown equation type for this variant")
    return(rep(NA_real_, nrow(data)))
  }
}


#' Apply the calibrated Bayesian Wykoff equation
#'
#' Uses MAP or posterior median estimates on the standardized scale.
#' Requires standardization parameters to transform raw predictors.
predict_calibrated_wykoff <- function(data, map_estimates, std_params) {

  # Extract MAP/posterior estimates
  get_param <- function(name) {
    row <- map_estimates %>% filter(variable == name)
    if (nrow(row) == 0) return(0)
    if ("estimate" %in% names(row)) return(row$estimate[1])
    if ("p50" %in% names(row)) return(row$p50[1])
    return(0)
  }

  # Get species intercepts b0[i]
  n_species <- sum(grepl("^b0\\[", map_estimates$variable))
  b0 <- sapply(1:n_species, function(i) get_param(paste0("b0[", i, "]")))

  # If b0 not directly available, reconstruct from z_b0
  if (all(b0 == 0)) {
    mu_b0 <- get_param("mu_b0")
    sigma_b0 <- get_param("sigma_b0")
    z_b0 <- sapply(1:n_species, function(i) get_param(paste0("z_b0[", i, "]")))
    b0 <- mu_b0 + sigma_b0 * z_b0
  }

  # Get shared slope coefficients
  b1  <- get_param("b1")   # ln(DBH)
  b2  <- get_param("b2")   # DBH^2
  b3  <- get_param("b3")   # ln(SI)
  b4  <- get_param("b4")   # SLOPE
  b5  <- get_param("b5")   # SLOPE^2
  b6  <- get_param("b6")   # SLOPE*sin(ASP)
  b7  <- get_param("b7")   # SLOPE*cos(ASP)
  b8  <- get_param("b8")   # ELEV
  b9  <- get_param("b9")   # ELEV^2
  b10 <- get_param("b10")  # CR
  b11 <- get_param("b11")  # CR^2
  b12 <- get_param("b12")  # BAL
  b13 <- get_param("b13")  # BA

  # Standardize predictors using saved parameters
  std <- function(x, param_name) {
    m <- std_params[[paste0(param_name, "_mean")]]
    s <- std_params[[paste0(param_name, "_sd")]]
    if (is.null(m) || is.null(s) || s == 0) return(x)
    (x - m) / s
  }

  ln_DBH_std    <- std(data$ln_DBH, "ln_DBH")
  DBH_sq_std    <- std(data$DBH_sq, "DBH_sq")
  ln_SI_std     <- std(data$ln_SI, "ln_SI")
  SLOPE_std     <- std(data$SLOPE, "SLOPE")
  SLOPE_sq_std  <- SLOPE_std^2
  SLOPE_SASP_std <- std(data$SLOPE_SASP, "SLOPE_SASP")
  SLOPE_CASP_std <- std(data$SLOPE_CASP, "SLOPE_CASP")
  ELEV_std      <- std(data$ELEV, "ELEV")
  ELEV_sq_std   <- ELEV_std^2
  BAL_std       <- std(data$BAL, "BAL")
  BA_std        <- std(data$BA, "BA")
  CR            <- data$CR_pct / 100
  CR_sq         <- CR^2

  # Species index for each observation
  sp_idx <- data$fvs_sp_idx

  # Compute predicted ln(DDS)
  ln_dds_pred <- b0[sp_idx] +
    b1  * ln_DBH_std +
    b2  * DBH_sq_std +
    b3  * ln_SI_std +
    b4  * SLOPE_std +
    b5  * SLOPE_sq_std +
    b6  * SLOPE_SASP_std +
    b7  * SLOPE_CASP_std +
    b8  * ELEV_std +
    b9  * ELEV_sq_std +
    b10 * CR +
    b11 * CR_sq +
    b12 * BAL_std +
    b13 * BA_std

  return(ln_dds_pred)
}


#' Compute prediction accuracy metrics
compute_metrics <- function(observed, predicted, group = NULL) {
  df <- tibble(obs = observed, pred = predicted)
  if (!is.null(group)) df$group <- group

  summarize_fn <- function(d) {
    d <- d %>% filter(is.finite(obs), is.finite(pred))
    n <- nrow(d)
    if (n < 5) return(tibble(n = n, rmse = NA, bias = NA, r2 = NA, mae = NA))

    resid <- d$obs - d$pred
    tibble(
      n = n,
      rmse = sqrt(mean(resid^2)),
      bias = mean(resid),
      mae = mean(abs(resid)),
      r2 = 1 - sum(resid^2) / sum((d$obs - mean(d$obs))^2)
    )
  }

  if (is.null(group)) {
    return(summarize_fn(df))
  } else {
    df %>%
      group_by(group) %>%
      group_modify(~ summarize_fn(.x)) %>%
      ungroup()
  }
}


# =============================================================================
# Process Each Variant
# =============================================================================

all_results <- list()

for (variant in variants_to_run) {

  logger::log_info("Processing variant: {variant}")
  cat("\n=== Variant:", toupper(variant), "===\n")

  # Paths
  variant_output <- file.path(output_base, variant)
  data_file <- file.path(calibration_dir, "data", "processed", variant, "diameter_growth.csv")
  config_file <- file.path(config_dir, paste0(variant, ".json"))
  posterior_file <- file.path(variant_output, "diameter_growth_posterior.csv")
  map_file <- file.path(variant_output, "diameter_growth_map.csv")
  std_file <- file.path(variant_output, "standardization_params.rds")

  # Check required files
  if (!file.exists(data_file)) {
    logger::log_warn("FIA data not found for {variant}: {data_file}")
    cat("  Skipping: FIA data not found\n")
    next
  }
  if (!file.exists(config_file)) {
    logger::log_warn("Config not found for {variant}")
    next
  }

  # Detect equation type early
  config_check <- fromJSON(config_file)
  eq_type <- detect_equation_type(config_check)
  cat("  Equation type:", eq_type, "\n")

  # Prefer MAP estimates (more stable than possibly unconverged HMC)
  if (file.exists(map_file)) {
    cal_estimates <- read_csv(map_file, show_col_types = FALSE)
    cal_source <- "MAP"
  } else if (file.exists(posterior_file)) {
    cal_estimates <- read_csv(posterior_file, show_col_types = FALSE)
    cal_source <- "HMC_posterior"
  } else {
    logger::log_warn("No calibration output for {variant}")
    cat("  Skipping: no calibration output\n")
    next
  }

  if (!file.exists(std_file)) {
    logger::log_warn("Standardization params not found for {variant}")
    cat("  Skipping: standardization params not found\n")
    next
  }

  # Load data
  cat("  Loading data...\n")
  fia_data <- fread(data_file) %>% as_tibble()
  config <- fromJSON(config_file)
  std_params <- readRDS(std_file)

  # Get species lookup
  fia_spcds <- unlist(config$categories$species_definitions$FIAJSP)
  fvs_spcodes <- unlist(config$categories$species_definitions$JSP)
  spcd_to_name <- setNames(trimws(fvs_spcodes), fia_spcds)

  # Add species name and filter to valid species indices
  n_sp_cfg <- switch(eq_type,
    western = max(length(unlist(config$categories$growth$WEIBB1)),
                  length(unlist(config$categories$other$DGLD)),
                  config$maxsp, na.rm = TRUE),
    central = length(unlist(config$categories$other$INTERC)),
    ne = length(unlist(config$categories$growth$B1)),
    config$maxsp
  )

  # Note on CR_pct: FIA stores CR as integer percent (e.g., 35 = 35%).
  # The data prep script does CR_pct = CR * 100, so CR_pct = 3500 means 35%.
  # For the Stan model (02_fit), CR = CR_pct / 100 = 35 (NOT 0 to 1).
  # For the default FVS equation, CR should be 0 to 1, so use CR_pct / 10000.
  fia_data <- fia_data %>%
    mutate(
      sp_name = spcd_to_name[as.character(SPCD)],
      # True crown ratio on 0 to 1 scale for default FVS equation
      CR_01 = pmin(pmax(CR_pct / 10000, 0), 1),
      # CR on 0 to ~99 scale as used by Stan fitting (CR_pct / 100)
      CR_stan = CR_pct / 100
    ) %>%
    filter(
      !is.na(fvs_sp_idx),
      fvs_sp_idx >= 1,
      fvs_sp_idx <= n_sp_cfg,
      is.finite(ln_DDS),
      is.finite(ln_DBH),
      CR_pct >= 0
    )

  # Subsample for speed if very large
  max_obs <- 50000
  if (nrow(fia_data) > max_obs) {
    set.seed(42)
    fia_data <- fia_data %>% slice_sample(n = max_obs)
  }

  cat("  Using", nrow(fia_data), "observations (", cal_source, ")\n")

  if (nrow(fia_data) < 10) {
    cat("  Skipping: insufficient valid observations after filtering\n")
    next
  }

  # --- Apply default FVS equation ---
  cat("  Applying default Wykoff equation...\n")
  fia_data$pred_default <- predict_default_wykoff(fia_data, config)

  # --- Apply calibrated equation ---
  cat("  Applying calibrated Wykoff equation...\n")
  fia_data$pred_calibrated <- predict_calibrated_wykoff(fia_data, cal_estimates, std_params)

  # --- Compute metrics ---
  cat("  Computing accuracy metrics...\n")

  # Overall
  metrics_default <- compute_metrics(fia_data$ln_DDS, fia_data$pred_default) %>%
    mutate(model = "default")
  metrics_calibrated <- compute_metrics(fia_data$ln_DDS, fia_data$pred_calibrated) %>%
    mutate(model = "calibrated")

  overall <- bind_rows(metrics_default, metrics_calibrated) %>%
    mutate(variant = toupper(variant))

  cat("\n  Overall Results:\n")
  cat(sprintf("    Default:    RMSE=%.4f  Bias=%.4f  R2=%.4f  (n=%d)\n",
              metrics_default$rmse, metrics_default$bias, metrics_default$r2, metrics_default$n))
  cat(sprintf("    Calibrated: RMSE=%.4f  Bias=%.4f  R2=%.4f  (n=%d)\n",
              metrics_calibrated$rmse, metrics_calibrated$bias, metrics_calibrated$r2, metrics_calibrated$n))

  # By species
  sp_default <- compute_metrics(fia_data$ln_DDS, fia_data$pred_default, fia_data$sp_name) %>%
    mutate(model = "default")
  sp_calibrated <- compute_metrics(fia_data$ln_DDS, fia_data$pred_calibrated, fia_data$sp_name) %>%
    mutate(model = "calibrated")

  sp_metrics <- bind_rows(sp_default, sp_calibrated) %>%
    rename(species = group) %>%
    mutate(variant = toupper(variant))

  # Save species metrics
  write_csv(sp_metrics, file.path(comparison_dir, paste0(variant, "_species_metrics.csv")))

  # --- Generate comparison plots ---
  cat("  Generating comparison plots...\n")

  # 1. Observed vs predicted scatter
  p_scatter <- fia_data %>%
    select(ln_DDS, pred_default, pred_calibrated) %>%
    pivot_longer(cols = starts_with("pred_"), names_to = "model", values_to = "predicted") %>%
    mutate(model = ifelse(model == "pred_default", "Default FVS", "Calibrated (Bayesian)")) %>%
    ggplot(aes(x = predicted, y = ln_DDS)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(trans = "log10") +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.8) +
    facet_wrap(~model) +
    labs(
      title = paste0(toupper(variant), ": Observed vs Predicted ln(DDS)"),
      subtitle = paste0("n = ", format(nrow(fia_data), big.mark = ","),
                        " | Calibration source: ", cal_source),
      x = "Predicted ln(DDS)", y = "Observed ln(DDS)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave(file.path(comparison_dir, paste0(variant, "_obs_vs_pred.png")),
         p_scatter, width = 12, height = 5.5, dpi = 150)

  # 2. Residual distributions
  fia_data <- fia_data %>%
    mutate(
      resid_default = ln_DDS - pred_default,
      resid_calibrated = ln_DDS - pred_calibrated
    )

  p_resid <- fia_data %>%
    select(resid_default, resid_calibrated) %>%
    pivot_longer(everything(), names_to = "model", values_to = "residual") %>%
    mutate(model = ifelse(model == "resid_default", "Default FVS", "Calibrated (Bayesian)")) %>%
    ggplot(aes(x = residual, fill = model)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Default FVS" = "#E41A1C", "Calibrated (Bayesian)" = "#377EB8")) +
    labs(
      title = paste0(toupper(variant), ": Residual Distributions"),
      x = "Residual [ln(DDS)]", y = "Density", fill = "Model"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(comparison_dir, paste0(variant, "_residuals.png")),
         p_resid, width = 8, height = 5, dpi = 150)

  # 3. Species level comparison bar chart (top 10 species by sample size)
  top_species <- fia_data %>%
    count(sp_name) %>%
    arrange(desc(n)) %>%
    slice_head(n = 10) %>%
    pull(sp_name)

  p_species <- sp_metrics %>%
    filter(species %in% top_species) %>%
    select(species, model, rmse, bias) %>%
    pivot_longer(cols = c(rmse, bias), names_to = "metric", values_to = "value") %>%
    mutate(metric = toupper(metric)) %>%
    ggplot(aes(x = reorder(species, value), y = value, fill = model)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~metric, scales = "free_x") +
    coord_flip() +
    scale_fill_manual(values = c("default" = "#E41A1C", "calibrated" = "#377EB8"),
                      labels = c("Default FVS", "Calibrated")) +
    labs(
      title = paste0(toupper(variant), ": Species Level Accuracy (Top 10 by n)"),
      x = "", y = "Value", fill = "Model"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(comparison_dir, paste0(variant, "_species_comparison.png")),
         p_species, width = 10, height = 6, dpi = 150)

  # 4. Growth curves across DBH range for top 5 species
  top5 <- fia_data %>%
    count(sp_name, fvs_sp_idx) %>%
    arrange(desc(n)) %>%
    slice_head(n = 5)

  # Median stand conditions
  med_si    <- median(fia_data$SI, na.rm = TRUE)
  med_slope <- median(fia_data$SLOPE, na.rm = TRUE)
  med_elev  <- median(fia_data$ELEV, na.rm = TRUE)
  med_bal   <- median(fia_data$BAL, na.rm = TRUE)
  med_ba    <- median(fia_data$BA, na.rm = TRUE)
  med_cr    <- 50  # percent
  med_qmd   <- if ("QMD" %in% names(fia_data)) median(fia_data$QMD, na.rm = TRUE) else 10
  med_tpa   <- if ("TPA" %in% names(fia_data)) median(fia_data$TPA, na.rm = TRUE) else 200

  curve_data <- expand_grid(
    DIA = seq(2, 30, by = 0.5),
    sp_row = 1:nrow(top5)
  ) %>%
    mutate(
      sp_name = top5$sp_name[sp_row],
      fvs_sp_idx = top5$fvs_sp_idx[sp_row],
      DIA_t1 = DIA,
      ln_DBH = log(DIA),
      DBH_sq = DIA^2,
      ln_SI = log(med_si),
      SI = med_si,
      SLOPE = med_slope,
      SLOPE_SASP = med_slope * sin(pi/4),  # assume 45 deg aspect
      SLOPE_CASP = med_slope * cos(pi/4),
      ELEV = med_elev,
      CR_pct = med_cr * 100,  # match data format: FIA_CR(0.5) * 100 = 50, stored as 5000
      CR_01 = med_cr / 100,   # 0 to 1 scale for default equation
      BAL = med_bal,
      BA = med_ba,
      QMD = med_qmd,
      TPA = med_tpa,
      SPCD = 0  # not needed for curve
    )

  curve_data$pred_default <- predict_default_wykoff(curve_data, config)
  curve_data$pred_calibrated <- predict_calibrated_wykoff(curve_data, cal_estimates, std_params)

  # Convert from ln(DDS) to annual diameter increment
  # DDS = (DIA_t2^2 - DIA_t1^2), so DIA_inc ~ sqrt(DDS + DIA^2) - DIA
  curve_long <- curve_data %>%
    select(DIA, sp_name, pred_default, pred_calibrated) %>%
    mutate(
      dds_default = exp(pred_default),
      dds_calibrated = exp(pred_calibrated),
      # Approximate annual DIA increment (assuming 10 year interval in data)
      dinc_default = (sqrt(pmax(dds_default + DIA^2, DIA^2)) - DIA) / 10,
      dinc_calibrated = (sqrt(pmax(dds_calibrated + DIA^2, DIA^2)) - DIA) / 10
    ) %>%
    select(DIA, sp_name, dinc_default, dinc_calibrated) %>%
    pivot_longer(cols = starts_with("dinc_"), names_to = "model", values_to = "dinc_annual") %>%
    mutate(model = ifelse(model == "dinc_default", "Default FVS", "Calibrated (Bayesian)"))

  p_curves <- curve_long %>%
    ggplot(aes(x = DIA, y = dinc_annual, color = model, linetype = model)) +
    geom_line(linewidth = 1) +
    facet_wrap(~sp_name, scales = "free_y") +
    scale_color_manual(values = c("Default FVS" = "#E41A1C", "Calibrated (Bayesian)" = "#377EB8")) +
    labs(
      title = paste0(toupper(variant), ": Annual Diameter Increment Curves"),
      subtitle = paste0("SI=", round(med_si), ", Slope=", round(med_slope, 1),
                        "%, Elev=", round(med_elev), " ft, BA=", round(med_ba),
                        " ft2/ac, CR=", med_cr, "%"),
      x = "DBH (inches)", y = "Annual Diameter Increment (inches/year)",
      color = "Model", linetype = "Model"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(comparison_dir, paste0(variant, "_growth_curves.png")),
         p_curves, width = 12, height = 8, dpi = 150)

  # 5. Residuals vs DBH
  p_resid_dbh <- fia_data %>%
    select(DIA_t1, resid_default, resid_calibrated) %>%
    pivot_longer(cols = starts_with("resid_"), names_to = "model", values_to = "residual") %>%
    mutate(model = ifelse(model == "resid_default", "Default FVS", "Calibrated (Bayesian)")) %>%
    ggplot(aes(x = DIA_t1, y = residual, color = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 1) +
    scale_color_manual(values = c("Default FVS" = "#E41A1C", "Calibrated (Bayesian)" = "#377EB8")) +
    labs(
      title = paste0(toupper(variant), ": Mean Residual by DBH"),
      x = "DBH (inches)", y = "Mean Residual [ln(DDS)]", color = "Model"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(comparison_dir, paste0(variant, "_resid_vs_dbh.png")),
         p_resid_dbh, width = 9, height = 5.5, dpi = 150)

  # Store results
  all_results[[variant]] <- list(
    overall = overall,
    species = sp_metrics,
    n_obs = nrow(fia_data),
    cal_source = cal_source
  )

  cat("  Done.\n")
}

# =============================================================================
# Cross Variant Summary
# =============================================================================

if (length(all_results) > 1) {

  cat("\n\n========================================\n")
  cat("CROSS VARIANT SUMMARY\n")
  cat("========================================\n\n")

  summary_table <- bind_rows(lapply(all_results, function(r) r$overall)) %>%
    select(variant, model, n, rmse, bias, r2, mae) %>%
    arrange(variant, model)

  print(summary_table, n = 100)

  write_csv(summary_table, file.path(comparison_dir, "cross_variant_summary.csv"))

  # Improvement table
  improvement <- summary_table %>%
    select(variant, model, rmse, bias, r2) %>%
    pivot_wider(names_from = model, values_from = c(rmse, bias, r2)) %>%
    mutate(
      rmse_change_pct = 100 * (rmse_calibrated - rmse_default) / rmse_default,
      bias_reduction = abs(bias_default) - abs(bias_calibrated),
      r2_change = r2_calibrated - r2_default
    )

  cat("\nImprovement Summary:\n")
  print(improvement %>% select(variant, rmse_change_pct, bias_reduction, r2_change), n = 100)

  write_csv(improvement, file.path(comparison_dir, "improvement_summary.csv"))

  # Cross variant bar chart
  p_cross <- summary_table %>%
    ggplot(aes(x = variant, y = rmse, fill = model)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("default" = "#E41A1C", "calibrated" = "#377EB8"),
                      labels = c("default" = "Default FVS", "calibrated" = "Calibrated (Bayesian)")) +
    labs(
      title = "Default vs Calibrated Wykoff: RMSE by FVS Variant",
      x = "FVS Variant", y = "RMSE [ln(DDS)]", fill = "Model"
    ) +
    theme_bw(base_size = 13)

  ggsave(file.path(comparison_dir, "cross_variant_rmse.png"),
         p_cross, width = 10, height = 5, dpi = 150)
}

cat("\n\nComparison complete. Outputs saved to:\n  ", comparison_dir, "\n")
logger::log_info("Comparison complete for {length(all_results)} variants")
