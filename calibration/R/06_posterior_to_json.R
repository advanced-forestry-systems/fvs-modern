#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: Convert All Posteriors to JSON Config
#
# Reads posterior samples from ALL fitted component models and writes
# calibrated parameters to JSON matching the original FVS config schema.
#
# Component models integrated:
#   1. Diameter growth (Wykoff) - DG parameters
#   2. Height-diameter (Chapman-Richards) - H-D parameters
#   3. Height increment (where applicable) - HG parameters
#   4. Mortality (logistic) - mortality parameters
#   5. Crown ratio change - BCR parameters
#
# Usage: Rscript calibration/R/06_posterior_to_json.R --variant ca
# =============================================================================

library(tidyverse)
library(jsonlite)
library(data.table)
library(logger)

# =============================================================================
# Parse Command Line Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
variant <- "ca"

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

project_root <- Sys.getenv("FVS_PROJECT_ROOT",
                             normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE))
calibration_dir <- file.path(project_root, "calibration")
output_dir <- file.path(calibration_dir, "output", "variants", variant)
config_dir <- file.path(project_root, "config")
calibrated_dir <- file.path(config_dir, "calibrated")

dir.create(calibrated_dir, showWarnings = FALSE, recursive = TRUE)

# Logging
log_file <- file.path(calibration_dir, "logs",
                       paste0("06_posterior_to_json_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Converting posteriors to JSON for variant {variant}")

# =============================================================================
# Load Original Config as Template
# =============================================================================

logger::log_info("Loading original FVS config...")

original_config_file <- file.path(config_dir, paste0(variant, ".json"))

if (!file.exists(original_config_file)) {
  logger::log_error("Original config not found: {original_config_file}")
  stop("Original config not found: ", original_config_file)
}

original_config <- fromJSON(original_config_file)
n_species <- original_config$maxsp

logger::log_info("Loaded config: {original_config$variant_name}, {n_species} species")

# Track which components were updated
components_updated <- character(0)
components_skipped <- character(0)

# =============================================================================
# Helper Functions
# =============================================================================

# Extract posterior median from a summary dataframe
extract_param <- function(posterior_df, param_name) {
  if (is.null(posterior_df)) return(NA_real_)

  # Try exact match
  match_row <- posterior_df %>% filter(variable == param_name)
  if (nrow(match_row) > 0) {
    # Use p50 if available (CmdStanR format), otherwise Estimate (brms format)
    if ("p50" %in% names(match_row)) return(match_row$p50[1])
    if ("Estimate" %in% names(match_row)) return(match_row$Estimate[1])
    if ("median" %in% names(match_row)) return(match_row$median[1])
  }
  return(NA_real_)
}

# Extract species specific posterior medians as a vector
extract_species_params <- function(posterior_df, base_name, n_sp) {
  if (is.null(posterior_df)) return(rep(NA_real_, n_sp))

  values <- rep(NA_real_, n_sp)
  for (idx in seq_len(n_sp)) {
    # Try indexed form: "b1[3]"
    param_key <- paste0(base_name, "[", idx, "]")
    val <- extract_param(posterior_df, param_key)

    if (is.na(val)) {
      # Try brms naming: "b_<name>_SPCD<code>"
      param_key2 <- paste0("b_", base_name)
      val <- extract_param(posterior_df, param_key2)
    }

    if (!is.na(val)) {
      values[idx] <- val
    }
  }

  # If no species specific values found, try global parameter
  if (all(is.na(values))) {
    global_val <- extract_param(posterior_df, base_name)
    if (!is.na(global_val)) {
      values <- rep(global_val, n_sp)
    }
  }

  return(values)
}

# Safely update a config parameter, preserving original where posterior is NA
safe_update <- function(config, category, param_name, new_values) {
  if (param_name %in% names(config$categories[[category]])) {
    original <- unlist(config$categories[[category]][[param_name]])
    # Replace only where new values are available
    updated <- ifelse(is.na(new_values), original, new_values)
    config$categories[[category]][[param_name]] <- as.list(updated)
    n_updated <- sum(!is.na(new_values))
    logger::log_info("  {param_name}: {n_updated}/{length(new_values)} values updated")
    return(list(config = config, n_updated = n_updated))
  } else {
    logger::log_info("  {param_name}: not found in config category '{category}'")
    return(list(config = config, n_updated = 0))
  }
}

# =============================================================================
# 1. Diameter Growth Parameters
# =============================================================================

logger::log_info("=== Processing diameter growth posteriors ===")

dg_posterior_file <- file.path(output_dir, "diameter_growth_posterior.csv")

if (file.exists(dg_posterior_file)) {
  dg_posterior <- read_csv(dg_posterior_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(dg_posterior)} diameter growth parameters")

  # Map Stan parameter names to FVS config parameter names
  dg_param_mapping <- list(
    DGLD   = "b1",    # ln(DBH) coefficient
    DGDS   = "b2",    # DBH^2 coefficient
    DGSITE = "b3",    # ln(SI) coefficient
    DGSLOP = "b4",    # SLOPE coefficient
    DGSLSQ = "b5",    # SLOPE^2 coefficient
    DGSASP = "b6",    # SLOPE*sin(ASP) coefficient
    DGCASP = "b7",    # SLOPE*cos(ASP) coefficient
    DGEL   = "b8",    # ELEV coefficient
    DGELSQ = "b9",    # ELEV^2 coefficient
    DGCR   = "b10",   # CR coefficient
    DGCRSQ = "b11",   # CR^2 coefficient
    DGBAL  = "b12",   # BAL coefficient
    DGPCCF = "b13"    # BA/CCF coefficient
  )

  total_dg_updated <- 0

  for (fvs_name in names(dg_param_mapping)) {
    stan_name <- dg_param_mapping[[fvs_name]]

    # Species specific parameters (b0, b1, b2, b3) vs global (b4-b13)
    if (stan_name %in% c("b1", "b2", "b3")) {
      new_vals <- extract_species_params(dg_posterior, stan_name, n_species)
    } else {
      # Global parameter: replicate across species
      global_val <- extract_param(dg_posterior, stan_name)
      new_vals <- rep(global_val, n_species)
    }

    result <- safe_update(original_config, "other", fvs_name, new_vals)
    original_config <- result$config
    total_dg_updated <- total_dg_updated + result$n_updated
  }

  # Also update intercept if available
  b0_vals <- extract_species_params(dg_posterior, "b0", n_species)
  if (any(!is.na(b0_vals))) {
    result <- safe_update(original_config, "growth", "WEIBB1", b0_vals)
    original_config <- result$config
  }

  logger::log_info("Diameter growth: {total_dg_updated} parameter values updated")
  components_updated <- c(components_updated, "diameter_growth")

} else {
  logger::log_warn("Diameter growth posterior not found: {dg_posterior_file}")
  components_skipped <- c(components_skipped, "diameter_growth")
}

# =============================================================================
# 2. Height-Diameter Parameters
# =============================================================================

logger::log_info("=== Processing height-diameter posteriors ===")

htdbh_summary_file <- file.path(output_dir, "height_diameter_summary.csv")

if (file.exists(htdbh_summary_file)) {
  htdbh_posterior <- read_csv(htdbh_summary_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(htdbh_posterior)} height-diameter parameters")

  # Chapman-Richards parameters: H = 4.5 + a * (1 - exp(-b * DBH))^c
  # brms stores population level effects with b_ prefix

  a_vals <- extract_species_params(htdbh_posterior, "b_a_Intercept", n_species)
  b_vals <- extract_species_params(htdbh_posterior, "b_b_Intercept", n_species)
  c_val <- extract_param(htdbh_posterior, "b_c_Intercept")

  # Store in height_diameter category
  if (!"height_diameter" %in% names(original_config$categories)) {
    original_config$categories$height_diameter <- list()
  }

  if (any(!is.na(a_vals))) {
    original_config$categories$height_diameter$HD_A <- as.list(a_vals)
    logger::log_info("  HD_A (asymptote): {sum(!is.na(a_vals))} species updated")
  }
  if (any(!is.na(b_vals))) {
    original_config$categories$height_diameter$HD_B <- as.list(b_vals)
    logger::log_info("  HD_B (rate): {sum(!is.na(b_vals))} species updated")
  }
  if (!is.na(c_val)) {
    original_config$categories$height_diameter$HD_C <- c_val
    logger::log_info("  HD_C (shape): updated to {round(c_val, 4)}")
  }

  components_updated <- c(components_updated, "height_diameter")

} else {
  logger::log_warn("Height-diameter posterior not found")
  components_skipped <- c(components_skipped, "height_diameter")
}

# =============================================================================
# 3. Height Increment Parameters (variants with explicit HG params only)
# =============================================================================

logger::log_info("=== Processing height increment posteriors ===")

htinc_posterior_file <- file.path(output_dir, "height_increment_posterior.csv")

if (file.exists(htinc_posterior_file)) {
  htinc_posterior <- read_csv(htinc_posterior_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(htinc_posterior)} height increment parameters")

  # Map Stan parameters to FVS HG config parameters
  hg_param_mapping <- list(
    HGLD = "b1",       # ln(DBH) coefficient (species specific)
    HGSC = "b0"        # Species specific intercept (site/competition)
  )

  for (fvs_name in names(hg_param_mapping)) {
    stan_name <- hg_param_mapping[[fvs_name]]

    if (fvs_name %in% names(original_config$categories$other)) {
      fvs_len <- length(unlist(original_config$categories$other[[fvs_name]]))
      new_vals <- extract_species_params(htinc_posterior, stan_name, fvs_len)
      result <- safe_update(original_config, "other", fvs_name, new_vals)
      original_config <- result$config
    }
  }

  # Update species group parameters (HGHC, HGLDD, HGH2)
  hg_grp_mapping <- list(
    HGHC  = "gamma_shape",
    HGLDD = "b1",
    HGH2  = "b2"
  )

  for (fvs_name in names(hg_grp_mapping)) {
    stan_name <- hg_grp_mapping[[fvs_name]]

    if (fvs_name %in% names(original_config$categories$other)) {
      fvs_len <- length(unlist(original_config$categories$other[[fvs_name]]))
      new_vals <- extract_species_params(htinc_posterior, stan_name, fvs_len)
      result <- safe_update(original_config, "other", fvs_name, new_vals)
      original_config <- result$config
    }
  }

  # Update shared fixed effects (HGLH)
  if ("HGLH" %in% names(original_config$categories$other)) {
    val <- extract_param(htinc_posterior, "b2")
    if (!is.na(val)) {
      fvs_len <- length(unlist(original_config$categories$other$HGLH))
      original_config$categories$other$HGLH <- as.list(rep(val, fvs_len))
      logger::log_info("  HGLH: updated to {round(val, 6)}")
    }
  }

  components_updated <- c(components_updated, "height_increment")

} else {
  logger::log_info("Height increment posterior not found (expected for most variants)")
  components_skipped <- c(components_skipped, "height_increment")
}

# =============================================================================
# 4. Mortality Parameters
# =============================================================================

logger::log_info("=== Processing mortality posteriors ===")

mort_summary_file <- file.path(output_dir, "mortality_summary.csv")

if (file.exists(mort_summary_file)) {
  mort_posterior <- read_csv(mort_summary_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(mort_posterior)} mortality parameters")

  # Store mortality calibration parameters
  if (!"mortality" %in% names(original_config$categories)) {
    original_config$categories$mortality <- list()
  }

  # Extract population level effects from brms output
  mort_params <- list(
    MORT_INT  = "b_Intercept",
    MORT_DBH  = "b_DBH_std",
    MORT_DBH2 = "b_IDBH_stdE2",
    MORT_BAL  = "b_BAL_std",
    MORT_CR   = "b_CR_std",
    MORT_SI   = "b_SI_std",
    MORT_BA   = "b_BA_std"
  )

  for (fvs_name in names(mort_params)) {
    brms_name <- mort_params[[fvs_name]]
    val <- extract_param(mort_posterior, brms_name)
    if (!is.na(val)) {
      original_config$categories$mortality[[fvs_name]] <- val
      logger::log_info("  {fvs_name}: {round(val, 6)}")
    }
  }

  # Extract species random effects
  sp_effects <- mort_posterior %>%
    filter(grepl("r_SPCD", variable)) %>%
    mutate(
      sp_idx = as.integer(gsub(".*\\[(.+),.*", "\\1", variable))
    )

  if (nrow(sp_effects) > 0) {
    est_col <- if ("Estimate" %in% names(sp_effects)) "Estimate" else
               if ("p50" %in% names(sp_effects)) "p50" else
               if ("median" %in% names(sp_effects)) "median" else NULL

    if (!is.null(est_col)) {
      sp_mort_adj <- sp_effects %>%
        select(sp_idx, value = !!sym(est_col)) %>%
        arrange(sp_idx) %>%
        pull(value)

      original_config$categories$mortality$SPECIES_ADJUSTMENT <- as.list(sp_mort_adj)
      logger::log_info("  Species adjustments: {length(sp_mort_adj)} values")
    }
  }

  components_updated <- c(components_updated, "mortality")

} else {
  logger::log_warn("Mortality posterior not found")
  components_skipped <- c(components_skipped, "mortality")
}

# =============================================================================
# 5. Crown Ratio Change Parameters
# =============================================================================

logger::log_info("=== Processing crown ratio posteriors ===")

cr_summary_file <- file.path(output_dir, "crown_ratio_summary.csv")

if (file.exists(cr_summary_file)) {
  cr_posterior <- read_csv(cr_summary_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(cr_posterior)} crown ratio parameters")

  # Map brms fixed effects to FVS BCR parameters
  cr_fixed_mapping <- list(
    b_Intercept = "BCR_INT",
    b_DBH_std   = "BCR1",
    b_BA_std    = "BCR2",
    b_BAL_std   = "BCR_BAL",
    b_CR_std    = "BCR3",
    b_SI_std    = "BCR4"
  )

  for (brms_name in names(cr_fixed_mapping)) {
    fvs_name <- cr_fixed_mapping[[brms_name]]
    val <- extract_param(cr_posterior, brms_name)

    if (!is.na(val)) {
      # BCR1-BCR4 are species indexed arrays in the crown category
      if (fvs_name %in% names(original_config$categories$crown)) {
        fvs_len <- length(unlist(original_config$categories$crown[[fvs_name]]))
        original_vals <- unlist(original_config$categories$crown[[fvs_name]])
        # Additive adjustment to original values
        original_config$categories$crown[[fvs_name]] <- as.list(
          original_vals + val
        )
        logger::log_info("  {fvs_name}: adjusted by {round(val, 6)} ({fvs_len} species)")
      } else {
        # Store as new parameter in crown_change category
        if (!"crown_change" %in% names(original_config$categories)) {
          original_config$categories$crown_change <- list()
        }
        original_config$categories$crown_change[[fvs_name]] <- val
        logger::log_info("  {fvs_name} (new): {round(val, 6)}")
      }
    }
  }

  # Extract species random effects for crown ratio
  sp_cr_effects <- cr_posterior %>%
    filter(grepl("r_SPCD", variable)) %>%
    mutate(
      sp_idx = as.integer(gsub(".*\\[(.+),.*", "\\1", variable))
    )

  if (nrow(sp_cr_effects) > 0) {
    est_col <- if ("Estimate" %in% names(sp_cr_effects)) "Estimate" else
               if ("p50" %in% names(sp_cr_effects)) "p50" else
               if ("median" %in% names(sp_cr_effects)) "median" else NULL

    if (!is.null(est_col)) {
      sp_cr_adj <- sp_cr_effects %>%
        select(sp_idx, value = !!sym(est_col)) %>%
        arrange(sp_idx) %>%
        pull(value)

      if (!"crown_change" %in% names(original_config$categories)) {
        original_config$categories$crown_change <- list()
      }
      original_config$categories$crown_change$SPECIES_ADJUSTMENT <- as.list(sp_cr_adj)
      logger::log_info("  CR species adjustments: {length(sp_cr_adj)} values")
    }
  }

  components_updated <- c(components_updated, "crown_ratio")

} else {
  logger::log_warn("Crown ratio posterior not found")
  components_skipped <- c(components_skipped, "crown_ratio")
}

# =============================================================================
# 6. Stand-Level Density Parameters (SDIMAX, BAMAX, self-thinning slope)
# =============================================================================

logger::log_info("=== Processing stand-level density posteriors ===")

sdimax_file <- file.path(output_dir, "species_sdimax_calibrated.csv")
bamax_file <- file.path(output_dir, "species_bamax_calibrated.csv")
slope_file <- file.path(output_dir, "self_thinning_slopes.csv")
overall_file <- file.path(output_dir, "stand_density_overall.csv")

if (file.exists(sdimax_file)) {
  sdimax_calibrated <- read_csv(sdimax_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(sdimax_calibrated)} species SDIMAX estimates")

  # Determine which SDI parameter this variant uses
  sdi_param <- NULL
  sdi_category <- NULL
  for (cat_name in names(original_config$categories)) {
    cat_data <- original_config$categories[[cat_name]]
    for (k in names(cat_data)) {
      if (k %in% c("SDICON", "R5SDI", "R4SDI", "FMSDI")) {
        sdi_param <- k
        sdi_category <- cat_name
      }
    }
  }

  if (!is.null(sdi_param)) {
    original_sdi <- unlist(original_config$categories[[sdi_category]][[sdi_param]])
    n_sdi <- length(original_sdi)

    # Get FIA species codes from config
    fia_species_codes <- unlist(original_config$categories$species_definitions$FIAJSP)
    spcd_to_idx <- setNames(seq_along(fia_species_codes), fia_species_codes)

    # Update SDIMAX values where calibrated data exists
    updated_sdi <- original_sdi
    n_updated_sdi <- 0

    for (i in seq_len(nrow(sdimax_calibrated))) {
      sp <- sdimax_calibrated$SPCD[i]
      idx <- spcd_to_idx[as.character(sp)]
      if (!is.na(idx) && idx <= n_sdi) {
        updated_sdi[idx] <- round(sdimax_calibrated$sdimax_combined[i])
        n_updated_sdi <- n_updated_sdi + 1
      }
    }

    original_config$categories[[sdi_category]][[sdi_param]] <- as.list(updated_sdi)
    logger::log_info("  {sdi_param}: {n_updated_sdi}/{n_sdi} species updated")
  }

  components_updated <- c(components_updated, "sdimax")
} else {
  logger::log_info("SDIMAX calibration not found (run 09_fit_stand_density.R first)")
  components_skipped <- c(components_skipped, "sdimax")
}

# Update BAMAX if applicable
if (file.exists(bamax_file)) {
  bamax_calibrated <- read_csv(bamax_file, show_col_types = FALSE) %>% as_tibble()
  logger::log_info("Loaded {nrow(bamax_calibrated)} species BAMAX estimates")

  bamax_param <- NULL
  bamax_category <- NULL
  for (cat_name in names(original_config$categories)) {
    cat_data <- original_config$categories[[cat_name]]
    for (k in names(cat_data)) {
      if (k %in% c("BAMAXA", "BAMAX1")) {
        bamax_param <- k
        bamax_category <- cat_name
      }
    }
  }

  if (!is.null(bamax_param)) {
    original_bamax <- unlist(original_config$categories[[bamax_category]][[bamax_param]])
    n_bamax <- length(original_bamax)

    fia_species_codes <- unlist(original_config$categories$species_definitions$FIAJSP)
    spcd_to_idx <- setNames(seq_along(fia_species_codes), fia_species_codes)

    updated_bamax <- original_bamax
    n_updated_bamax <- 0

    for (i in seq_len(nrow(bamax_calibrated))) {
      sp <- bamax_calibrated$SPCD[i]
      idx <- spcd_to_idx[as.character(sp)]
      if (!is.na(idx) && idx <= n_bamax) {
        updated_bamax[idx] <- round(bamax_calibrated$bamax_combined[i])
        n_updated_bamax <- n_updated_bamax + 1
      }
    }

    original_config$categories[[bamax_category]][[bamax_param]] <- as.list(updated_bamax)
    logger::log_info("  {bamax_param}: {n_updated_bamax}/{n_bamax} species updated")
  }

  components_updated <- c(components_updated, "bamax")
}

# Update self-thinning slopes in BNORML if available
if (file.exists(slope_file)) {
  slopes <- read_csv(slope_file, show_col_types = FALSE) %>% as_tibble()

  if ("BNORML" %in% names(original_config$categories$other)) {
    # BNORML is typically a lookup table by diameter class, not by species
    # Log the calibrated slopes for reference but don't overwrite BNORML structure
    # (BNORML format varies by variant and is used for DBH class conversion)
    overall_results <- if (file.exists(overall_file)) {
      read_csv(overall_file, show_col_types = FALSE)
    } else NULL

    if (!is.null(overall_results)) {
      if (!"self_thinning" %in% names(original_config$categories)) {
        original_config$categories$self_thinning <- list()
      }
      original_config$categories$self_thinning$CALIBRATED_SLOPE <- overall_results$pop_slope_bayes[1]
      original_config$categories$self_thinning$REINEKE_STANDARD <- -1.605
      logger::log_info("  Self-thinning slope: {round(overall_results$pop_slope_bayes[1], 4)}")
    }
  }
}

# =============================================================================
# Add Calibration Metadata
# =============================================================================

logger::log_info("Adding calibration metadata...")

original_config$calibration <- list(
  date = as.character(Sys.Date()),
  method = "Bayesian MCMC (CmdStanR + brms) and quantile regression",
  fia_data_source = "rFIA package",
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  prior_type = "Informative (from original FVS parameters)",
  components_updated = as.list(components_updated),
  components_skipped = as.list(components_skipped),
  notes = paste0(
    "Calibration based on FIA data. ",
    length(components_updated), " of 7 component models integrated ",
    "(tree: DG, H-D, HG, mortality, CR; stand: SDIMAX, BAMAX)."
  )
)

# =============================================================================
# Precompute Schema-Free Per-Species Keyword Multipliers
# =============================================================================
#
# Attach a `calibration_multipliers` block that the Python keyword emitter
# (config/config_loader.py :: generate_keywords) consumes directly. This
# decouples the emitter from the per-variant growth/mortality model schema and
# fixes the "at most one keyword block per variant" propagation gap (issue #54).

logger::log_info("Computing per-species keyword multipliers...")

tryCatch({
  source(file.path(calibration_dir, "R", "multipliers.R"), local = TRUE)
  original_config$calibration_multipliers <-
    compute_calibration_multipliers(output_dir, original_config)
  prov <- original_config$calibration_multipliers$provenance
  logger::log_info(
    "Attached calibration_multipliers (DG n={dg}, HG n={hg}, mort RE n={mr})",
    dg = prov$dds_n_species %||% 0,
    hg = prov$htg_n_species %||% 0,
    mr = prov$mort_n_re %||% 0
  )
}, error = function(e) {
  logger::log_warn("calibration_multipliers block skipped: {conditionMessage(e)}")
})

# =============================================================================
# Save Calibrated Config
# =============================================================================

logger::log_info("Saving calibrated config...")

output_file <- file.path(calibrated_dir, paste0(variant, ".json"))
json_output <- toJSON(original_config, pretty = TRUE, auto_unbox = TRUE)
write(json_output, file = output_file)

logger::log_info("Saved calibrated config to {output_file}")

# =============================================================================
# 8. Export Full Posterior Draws for Uncertainty Propagation
# =============================================================================
#
# In addition to the point estimate (median) JSON, we export the full set
# of posterior draws so that the Python runtime can sample parameter vectors
# for Monte Carlo uncertainty propagation through FVS projections.
#
# The draws file stores N posterior samples (default 500) for each component
# model. At runtime, the UncertaintyEngine picks a single draw index and
# assembles a complete parameter vector from correlated draws across all
# components, preserving the joint posterior structure.
# =============================================================================

logger::log_info("=== Exporting posterior draws for uncertainty propagation ===")

n_draws_to_keep <- 500  # subsample from the full MCMC chain
draws_list <- list()

# Helper: subsample posterior draws from a CSV of MCMC output
# Expects columns: .draw (or .iteration), variable, value
# OR wide format with one row per draw and columns per parameter
subsample_draws <- function(draws_df, n_keep = n_draws_to_keep) {
  if (is.null(draws_df) || nrow(draws_df) == 0) return(NULL)

  # Detect format: long (variable, value) vs wide (parameters as columns)
  if ("variable" %in% names(draws_df) && "value" %in% names(draws_df)) {
    # Long format from CmdStanR summary: not raw draws, skip
    return(NULL)
  }

  # Wide format: each row is a posterior draw, columns are parameters
  n_available <- nrow(draws_df)
  if (n_available <= n_keep) {
    return(draws_df)
  }

  # Systematic thinning to preserve chain mixing
  idx <- round(seq(1, n_available, length.out = n_keep))
  return(draws_df[idx, ])
}

# Load raw posterior draws files (these are the actual MCMC samples,
# not the summary statistics loaded earlier)

# Diameter growth draws
dg_draws_file <- file.path(output_dir, "diameter_growth_samples.rds")
if (!file.exists(dg_draws_file)) {
  dg_draws_file <- file.path(output_dir, "diameter_growth_draws.csv")
}

if (file.exists(dg_draws_file)) {
  logger::log_info("Loading diameter growth draws...")
  if (grepl("\\.rds$", dg_draws_file)) {
    dg_draws_raw <- readRDS(dg_draws_file)
  } else {
    dg_draws_raw <- fread(dg_draws_file)
  }
  dg_sub <- subsample_draws(as.data.frame(dg_draws_raw))
  if (!is.null(dg_sub)) {
    # Convert to list of lists for JSON serialization
    # Each element is a named list of parameter values for one draw
    draws_list$diameter_growth <- lapply(seq_len(nrow(dg_sub)), function(i) {
      as.list(dg_sub[i, , drop = FALSE])
    })
    logger::log_info("  Kept {length(draws_list$diameter_growth)} diameter growth draws")
  }
}

# Height diameter draws
hd_draws_file <- file.path(output_dir, "height_diameter_samples.rds")
if (!file.exists(hd_draws_file)) {
  hd_draws_file <- file.path(output_dir, "height_diameter_draws.csv")
}

if (file.exists(hd_draws_file)) {
  logger::log_info("Loading height diameter draws...")
  if (grepl("\\.rds$", hd_draws_file)) {
    hd_draws_raw <- readRDS(hd_draws_file)
  } else {
    hd_draws_raw <- fread(hd_draws_file)
  }
  hd_sub <- subsample_draws(as.data.frame(hd_draws_raw))
  if (!is.null(hd_sub)) {
    draws_list$height_diameter <- lapply(seq_len(nrow(hd_sub)), function(i) {
      as.list(hd_sub[i, , drop = FALSE])
    })
    logger::log_info("  Kept {length(draws_list$height_diameter)} height diameter draws")
  }
}

# Height increment draws (6 variants only)
hi_draws_file <- file.path(output_dir, "height_increment_samples.rds")
if (!file.exists(hi_draws_file)) {
  hi_draws_file <- file.path(output_dir, "height_increment_draws.csv")
}

if (file.exists(hi_draws_file)) {
  logger::log_info("Loading height increment draws...")
  if (grepl("\\.rds$", hi_draws_file)) {
    hi_draws_raw <- readRDS(hi_draws_file)
  } else {
    hi_draws_raw <- fread(hi_draws_file)
  }
  hi_sub <- subsample_draws(as.data.frame(hi_draws_raw))
  if (!is.null(hi_sub)) {
    draws_list$height_increment <- lapply(seq_len(nrow(hi_sub)), function(i) {
      as.list(hi_sub[i, , drop = FALSE])
    })
    logger::log_info("  Kept {length(draws_list$height_increment)} height increment draws")
  }
}

# Mortality draws
mort_draws_file <- file.path(output_dir, "mortality_samples.rds")
if (!file.exists(mort_draws_file)) {
  mort_draws_file <- file.path(output_dir, "mortality_draws.csv")
}

if (file.exists(mort_draws_file)) {
  logger::log_info("Loading mortality draws...")
  if (grepl("\\.rds$", mort_draws_file)) {
    mort_draws_raw <- readRDS(mort_draws_file)
  } else {
    mort_draws_raw <- fread(mort_draws_file)
  }
  mort_sub <- subsample_draws(as.data.frame(mort_draws_raw))
  if (!is.null(mort_sub)) {
    draws_list$mortality <- lapply(seq_len(nrow(mort_sub)), function(i) {
      as.list(mort_sub[i, , drop = FALSE])
    })
    logger::log_info("  Kept {length(draws_list$mortality)} mortality draws")
  }
}

# Crown ratio draws
cr_draws_file <- file.path(output_dir, "crown_ratio_samples.rds")
if (!file.exists(cr_draws_file)) {
  cr_draws_file <- file.path(output_dir, "crown_ratio_draws.csv")
}

if (file.exists(cr_draws_file)) {
  logger::log_info("Loading crown ratio draws...")
  if (grepl("\\.rds$", cr_draws_file)) {
    cr_draws_raw <- readRDS(cr_draws_file)
  } else {
    cr_draws_raw <- fread(cr_draws_file)
  }
  cr_sub <- subsample_draws(as.data.frame(cr_draws_raw))
  if (!is.null(cr_sub)) {
    draws_list$crown_ratio <- lapply(seq_len(nrow(cr_sub)), function(i) {
      as.list(cr_sub[i, , drop = FALSE])
    })
    logger::log_info("  Kept {length(draws_list$crown_ratio)} crown ratio draws")
  }
}

# Stand density draws (SDIMAX posterior from Bayesian quantile regression)
sdi_draws_file <- file.path(output_dir, "stand_density_samples.rds")
if (!file.exists(sdi_draws_file)) {
  sdi_draws_file <- file.path(output_dir, "stand_density_draws.csv")
}

if (file.exists(sdi_draws_file)) {
  logger::log_info("Loading stand density draws...")
  if (grepl("\\.rds$", sdi_draws_file)) {
    sdi_draws_raw <- readRDS(sdi_draws_file)
  } else {
    sdi_draws_raw <- fread(sdi_draws_file)
  }
  sdi_sub <- subsample_draws(as.data.frame(sdi_draws_raw))
  if (!is.null(sdi_sub)) {
    draws_list$stand_density <- lapply(seq_len(nrow(sdi_sub)), function(i) {
      as.list(sdi_sub[i, , drop = FALSE])
    })
    logger::log_info("  Kept {length(draws_list$stand_density)} stand density draws")
  }
}

# Save draws JSON
if (length(draws_list) > 0) {
  draws_output <- list(
    variant = variant,
    n_draws = n_draws_to_keep,
    calibration_date = as.character(Sys.Date()),
    components = names(draws_list),
    draws = draws_list
  )

  draws_file <- file.path(calibrated_dir, paste0(variant, "_draws.json"))
  draws_json <- toJSON(draws_output, pretty = FALSE, auto_unbox = TRUE, digits = 6)
  write(draws_json, file = draws_file)
  logger::log_info("Saved {length(draws_list)} component draw sets to {draws_file}")

  # Also report file size (draws files can be large)
  fsize_mb <- file.size(draws_file) / 1024^2
  logger::log_info("Draws file size: {round(fsize_mb, 1)} MB")
} else {
  logger::log_info("No raw posterior draws found; uncertainty propagation requires")
  logger::log_info("  saving MCMC draws in the fitting scripts (02-05, 09).")
  logger::log_info("  Add posterior::as_draws_df(fit) output to each fitting script.")
}

# =============================================================================
# Summary Report
# =============================================================================

cat("\n")
cat("==========================================\n")
cat("Posterior to JSON Conversion Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("Original config:", original_config_file, "\n")
cat("Calibrated config:", output_file, "\n")
cat("Species:", n_species, "\n")
cat("Calibration date:", as.character(Sys.Date()), "\n")
cat("\nComponents updated:\n")
for (comp in components_updated) {
  cat("  +", comp, "\n")
}
if (length(components_skipped) > 0) {
  cat("\nComponents skipped (posterior not available):\n")
  for (comp in components_skipped) {
    cat("  -", comp, "\n")
  }
}
cat("\n")

logger::log_info("Posterior to JSON conversion complete: {length(components_updated)} components integrated")
