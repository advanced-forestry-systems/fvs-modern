#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: Fetch and Prepare Stand-Level FIA Data
#
# Computes stand level metrics from FIA tree lists for calibrating
# density dependent parameters: SDIMAX, BAMAX, self-thinning slope,
# and the Reineke exponent.
#
# For each FIA plot/condition, computes:
#   - Trees per acre (TPA) and basal area (BA, ft2/ac)
#   - Quadratic mean diameter (QMD, in)
#   - Stand Density Index (SDI) using Reineke's equation
#   - Relative density (RD = BA / sqrt(QMD))
#   - Top height (average of 40 tallest TPA)
#   - Species composition (leading species, species group SDI)
#   - Mortality indicators (annual mortality rate by BA)
#
# Output: calibration/data/processed/<variant>/stand_density.csv
#
# Usage: Rscript calibration/R/08_fetch_stand_data.R --variant ne
# =============================================================================

library(tidyverse)
library(data.table)
library(rFIA)
library(jsonlite)
library(logger)

# =============================================================================
# Parse Command Line Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
variant <- "ne"
fia_dir <- NULL

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variant <- args[i + 1]
    }
    if (args[i] == "--fia-dir" & i < length(args)) {
      fia_dir <- args[i + 1]
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

# If fia_dir not passed as arg, check environment variable
if (is.null(fia_dir)) {
  env_fia <- Sys.getenv("FVS_FIA_DATA_DIR", "")
  if (nchar(env_fia) > 0) fia_dir <- env_fia
}
processed_dir <- file.path(calibration_dir, "data", "processed", variant)
config_dir <- file.path(project_root, "config")

dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(calibration_dir, "logs",
                       paste0("08_stand_data_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Preparing stand-level data for variant {variant}")

# =============================================================================
# Variant to State Mapping
# =============================================================================

# Same mapping used by 01_fetch_fia_data.R
variant_states <- list(
  ne = c("ME", "NH", "VT", "MA", "CT", "RI", "NY", "PA", "NJ", "DE", "MD", "WV"),
  ls = c("MI", "MN", "WI"),
  cs = c("IN", "IL", "IA", "MO", "OH", "KY"),
  sn = c("VA", "NC", "SC", "GA", "FL", "AL", "MS", "TN", "AR", "LA", "TX", "OK"),
  ca = c("CA"),
  oc = c("CA"),   # Oregon/California overlap
  pn = c("WA", "OR"),
  wc = c("WA", "OR"),
  nc = c("WA", "OR"),
  ie = c("ID", "MT"),
  ci = c("ID", "MT"),
  kt = c("MT"),
  em = c("ID", "MT", "WY", "CO", "UT", "NV", "AZ", "NM"),
  cr = c("CO", "AZ", "NM"),
  so = c("OR", "CA"),
  bm = c("OR"),
  ec = c("WA", "OR"),
  tt = c("UT", "NV"),
  ut = c("UT", "NV"),
  ws = c("CA", "OR", "WA"),
  op = c("OR"),
  ak = c("AK"),
  acd = c("ME", "NH", "VT"),
  bc = c("WA", "OR", "ID", "MT"),  # Placeholder, BC uses Canadian data
  on = c("MN", "WI", "MI")   # Placeholder, ON uses Canadian data
)

variant_upper <- toupper(variant)

states <- variant_states[[variant]]
if (is.null(states)) {
  logger::log_error("Unknown variant: {variant}")
  stop("Unknown variant: ", variant)
}

logger::log_info("Variant {variant} maps to states: {paste(states, collapse = ', ')}")
logger::log_info("Will filter on FIA variant column = '{variant_upper}' when available")

# =============================================================================
# Load FVS Config for Species Mapping
# =============================================================================

config_file <- file.path(config_dir, paste0(variant, ".json"))
config <- fromJSON(config_file)

# Get FIA species codes from the config
fia_species <- unlist(config$categories$species_definitions$FIAJSP)
n_species <- config$maxsp

logger::log_info("Config has {n_species} species, {sum(fia_species > 0)} with FIA codes")

# =============================================================================
# Download and Process FIA Data
# =============================================================================

logger::log_info("Loading FIA data for {length(states)} states...")

# Determine FIA source
if (is.null(fia_dir)) {
  fia_source_dir <- file.path(calibration_dir, "data", "fia_raw")
  dir.create(fia_source_dir, showWarnings = FALSE, recursive = TRUE)
  logger::log_info("Download mode: FIA data cached in {fia_source_dir}")
} else {
  fia_source_dir <- fia_dir
  logger::log_info("Local FIA mode: reading from {fia_source_dir}")
}

all_stand_data <- list()

# Detect ENTIRE_FIA flat file format
entire_tree_file <- file.path(fia_source_dir, "ENTIRE_TREE.csv")
entire_plot_file <- file.path(fia_source_dir, "ENTIRE_PLOT.csv")
entire_cond_file <- file.path(fia_source_dir, "ENTIRE_COND.csv")
use_entire_mode <- !is.null(fia_dir) &&
                   file.exists(entire_tree_file) &&
                   file.exists(entire_plot_file) &&
                   file.exists(entire_cond_file)

# State FIPS lookup
state_fips <- c(
  AL=1, AK=2, AZ=4, AR=5, CA=6, CO=8, CT=9, DE=10, FL=12, GA=13,
  HI=15, ID=16, IL=17, IN=18, IA=19, KS=20, KY=21, LA=22, ME=23,
  MD=24, MA=25, MI=26, MN=27, MS=28, MO=29, MT=30, NE=31, NV=32,
  NH=33, NJ=34, NM=35, NY=36, NC=37, ND=38, OH=39, OK=40, OR=41,
  PA=42, RI=44, SC=45, SD=46, TN=47, TX=48, UT=49, VT=50, VA=51,
  WA=53, WV=54, WI=55, WY=56
)
target_statecds <- unname(state_fips[states])
target_statecds <- target_statecds[!is.na(target_statecds)]

# Pre load ENTIRE files once if in ENTIRE mode (avoid re reading 13 GB per state)
entire_tree_dt <- NULL
entire_cond_dt <- NULL
entire_plot_dt <- NULL

if (use_entire_mode) {
  logger::log_info("Detected ENTIRE_FIA flat file format. Reading files once...")
  entire_tree_dt <- fread(entire_tree_file)[STATECD %in% target_statecds]
  entire_cond_dt <- fread(entire_cond_file)[STATECD %in% target_statecds]
  entire_plot_dt <- fread(entire_plot_file)[STATECD %in% target_statecds]
  logger::log_info("  Loaded {nrow(entire_tree_dt)} trees, {nrow(entire_cond_dt)} conditions, {nrow(entire_plot_dt)} plots")
}

for (state in states) {
  logger::log_info("Processing state: {state}")

  tryCatch({
    if (use_entire_mode) {
      # ENTIRE_FIA flat file mode (pre loaded above)
      stcd <- state_fips[state]
      if (is.na(stcd)) stop(paste("Unknown state:", state))
      db <- list(
        TREE = as_tibble(entire_tree_dt[STATECD == stcd]),
        COND = as_tibble(entire_cond_dt[STATECD == stcd]),
        PLOT = as_tibble(entire_plot_dt[STATECD == stcd])
      )
    } else if (!is.null(fia_dir)) {
      # State subdirectory mode
      state_subdir <- file.path(fia_source_dir, state)
      if (dir.exists(state_subdir) && length(list.files(state_subdir)) > 0) {
        db <- readFIA(dir = state_subdir, common = TRUE)
      } else {
        db <- readFIA(dir = fia_source_dir, states = state, common = TRUE)
      }
    } else {
      # Download mode
      state_dir <- file.path(fia_source_dir, state)
      if (!dir.exists(state_dir) || length(list.files(state_dir)) == 0) {
        logger::log_info("Downloading FIA data for {state}...")
        getFIA(states = state, dir = fia_source_dir, load = FALSE)
      }
      db <- readFIA(state_dir, common = TRUE)
    }

    # Extract tree-level data with plot conditions
    tree_dat <- db$TREE %>%
      as_tibble() %>%
      filter(
        STATUSCD %in% c(1, 2),     # Live and recently dead trees
        DIA >= 1.0,                 # Minimum 1 inch DBH
        !is.na(TPA_UNADJ),
        TPA_UNADJ > 0
      ) %>%
      select(
        CN, PLT_CN, CONDID, INVYR, CYCLE, SUBCYCLE,
        STATUSCD, SPCD, DIA, HT, CR, TPA_UNADJ,
        ACTUALHT, BHAGE, TOTAGE
      )

    # Get plot/condition info
    cond_dat <- db$COND %>%
      as_tibble() %>%
      select(any_of(c(
        "PLT_CN", "CONDID", "INVYR",
        "FORTYPCD", "FLDTYPCD", "STDAGE", "STDSZCD", "STDORGCD",
        "SITECLCD", "SICOND", "SISP", "SLOPE", "ASPECT", "ELEV",
        "BALIVE", "LIVE_CANOPY_CVR_PCT"
      )))

    # Extract PLOT table (grab variant column if it exists)
    plot_raw <- db$PLOT %>% as_tibble()
    plot_cols <- c("CN", "INVYR", "LAT", "LON", "ELEV")

    var_col <- intersect(
      names(plot_raw),
      c("VARIANT", "FVS_VARIANT", "FVSVARIANT", "FVSVARIANTCD",
        "variant", "fvs_variant")
    )
    if (length(var_col) > 0) {
      plot_cols <- c(plot_cols, var_col[1])
      logger::log_info("Found FIA variant column: {var_col[1]}")
    }

    plot_dat <- plot_raw %>%
      select(any_of(plot_cols)) %>%
      rename(PLT_CN = CN)

    if (length(var_col) > 0) {
      plot_dat <- plot_dat %>% rename(FIA_VARIANT = !!sym(var_col[1]))
    }

    # Join tree data to plot/condition
    combined <- tree_dat %>%
      left_join(cond_dat, by = c("PLT_CN", "CONDID", "INVYR")) %>%
      left_join(plot_dat, by = c("PLT_CN", "INVYR"))

    # Filter to plots assigned to this FVS variant (precise)
    if ("FIA_VARIANT" %in% names(combined)) {
      combined <- combined %>%
        mutate(FIA_VARIANT = trimws(toupper(FIA_VARIANT))) %>%
        filter(FIA_VARIANT == variant_upper)
      logger::log_info("Filtered by FIA variant: {nrow(combined)} records for {variant_upper}")
    }

    # =========================================================================
    # Compute Stand-Level Summaries per Plot/Condition/Inventory
    # =========================================================================

    # Normalize ELEV column name (ELEV.x if from COND join, ELEV if only from PLOT)
    if ("ELEV.x" %in% names(combined)) {
      combined$.elev_col <- combined$ELEV.x
    } else if ("ELEV" %in% names(combined)) {
      combined$.elev_col <- combined$ELEV
    } else {
      combined$.elev_col <- NA_real_
    }

    stand_metrics <- combined %>%
      filter(STATUSCD == 1) %>%  # Live trees only for density
      group_by(PLT_CN, CONDID, INVYR, CYCLE) %>%
      summarise(
        # Basic stand metrics
        tpa = sum(TPA_UNADJ, na.rm = TRUE),
        ba_ft2ac = sum(TPA_UNADJ * 0.005454 * DIA^2, na.rm = TRUE),
        qmd = sqrt(ba_ft2ac / tpa / 0.005454),
        n_trees = n(),

        # SDI using Reineke's equation: SDI = TPA * (QMD / 10)^1.605
        sdi = tpa * (qmd / 10)^1.605,

        # Relative density (Curtis): RD = BA / sqrt(QMD)
        rd_curtis = ba_ft2ac / sqrt(qmd),

        # Top height: average height of 40 largest TPA equivalent trees
        top_ht = {
          ht_sorted <- sort(HT[!is.na(HT)], decreasing = TRUE)
          if (length(ht_sorted) >= 5) mean(ht_sorted[1:min(5, length(ht_sorted))])
          else if (length(ht_sorted) > 0) mean(ht_sorted)
          else NA_real_
        },

        # Species composition
        leading_spcd = {
          sp_ba <- tapply(TPA_UNADJ * 0.005454 * DIA^2, SPCD, sum, na.rm = TRUE)
          as.integer(names(which.max(sp_ba)))
        },
        n_species = n_distinct(SPCD),
        leading_sp_pct = {
          sp_ba <- tapply(TPA_UNADJ * 0.005454 * DIA^2, SPCD, sum, na.rm = TRUE)
          100 * max(sp_ba, na.rm = TRUE) / sum(sp_ba, na.rm = TRUE)
        },

        # Site information
        si_cond = first(c(na.omit(SICOND), NA_real_)),
        site_class = first(c(na.omit(SITECLCD), NA_real_)),
        fortypcd = first(c(na.omit(FORTYPCD), NA_real_)),
        slope = first(c(na.omit(SLOPE), NA_real_)),
        aspect = first(c(na.omit(ASPECT), NA_real_)),
        elev = first(c(na.omit(.elev_col), NA_real_)),
        lat = first(c(na.omit(LAT), NA_real_)),
        lon = first(c(na.omit(LON), NA_real_)),
        std_age = first(c(na.omit(STDAGE), NA_real_)),

        .groups = "drop"
      ) %>%
      filter(
        tpa > 10,           # Minimum tree density
        ba_ft2ac > 5,       # Minimum basal area
        qmd > 1,            # Minimum QMD
        n_trees >= 5        # At least 5 measured trees per condition
      )

    # =========================================================================
    # Compute Species-Level SDI within Each Plot
    # =========================================================================

    species_sdi <- combined %>%
      filter(STATUSCD == 1) %>%
      group_by(PLT_CN, CONDID, INVYR, SPCD) %>%
      summarise(
        sp_tpa = sum(TPA_UNADJ, na.rm = TRUE),
        sp_ba = sum(TPA_UNADJ * 0.005454 * DIA^2, na.rm = TRUE),
        sp_qmd = sqrt(sp_ba / sp_tpa / 0.005454),
        sp_sdi = sp_tpa * (sp_qmd / 10)^1.605,
        sp_n_trees = n(),
        .groups = "drop"
      ) %>%
      filter(sp_n_trees >= 3, sp_tpa > 5)

    # =========================================================================
    # Compute Mortality Rates from Remeasurement
    # =========================================================================

    # Find remeasurement pairs for plots
    plot_cycles <- combined %>%
      distinct(PLT_CN, CONDID, INVYR, CYCLE) %>%
      group_by(PLT_CN, CONDID) %>%
      filter(n() >= 2) %>%
      arrange(INVYR) %>%
      mutate(
        meas_order = row_number(),
        next_invyr = lead(INVYR),
        interval = lead(INVYR) - INVYR
      ) %>%
      filter(!is.na(next_invyr)) %>%
      ungroup()

    if (nrow(plot_cycles) > 0) {
      # Compute mortality BA for remeasured plots
      mort_ba <- combined %>%
        filter(STATUSCD == 2) %>%  # Dead trees
        inner_join(plot_cycles %>% select(PLT_CN, CONDID, next_invyr),
                   by = c("PLT_CN", "CONDID", "INVYR" = "next_invyr")) %>%
        group_by(PLT_CN, CONDID, INVYR) %>%
        summarise(
          mort_ba_ft2ac = sum(TPA_UNADJ * 0.005454 * DIA^2, na.rm = TRUE),
          mort_tpa = sum(TPA_UNADJ, na.rm = TRUE),
          .groups = "drop"
        )

      # Merge mortality with stand metrics
      stand_with_mort <- stand_metrics %>%
        left_join(
          plot_cycles %>% select(PLT_CN, CONDID, INVYR, interval),
          by = c("PLT_CN", "CONDID", "INVYR")
        ) %>%
        left_join(mort_ba, by = c("PLT_CN", "CONDID", "INVYR")) %>%
        mutate(
          mort_ba_ft2ac = replace_na(mort_ba_ft2ac, 0),
          mort_tpa = replace_na(mort_tpa, 0),
          annual_mort_rate_ba = ifelse(
            !is.na(interval) & interval > 0 & ba_ft2ac > 0,
            1 - (1 - mort_ba_ft2ac / (ba_ft2ac + mort_ba_ft2ac))^(1 / interval),
            NA_real_
          )
        )
    } else {
      stand_with_mort <- stand_metrics %>%
        mutate(
          interval = NA_real_,
          mort_ba_ft2ac = NA_real_,
          mort_tpa = NA_real_,
          annual_mort_rate_ba = NA_real_
        )
    }

    stand_with_mort$state <- state
    all_stand_data[[state]] <- stand_with_mort

    # Also save species SDI data
    species_sdi$state <- state
    all_stand_data[[paste0(state, "_species")]] <- species_sdi

    logger::log_info("State {state}: {nrow(stand_with_mort)} plot-conditions, {nrow(species_sdi)} species-level records")

  }, error = function(e) {
    logger::log_error("Failed to process state {state}: {conditionMessage(e)}")
    logger::log_error("Traceback: {paste(capture.output(traceback()), collapse = ' | ')}")
  })
}

# =============================================================================
# Combine and Save
# =============================================================================

# Stand-level data
stand_data <- bind_rows(
  all_stand_data[!grepl("_species$", names(all_stand_data))]
) %>%
  mutate(variant = variant)

# Species-level SDI data
species_data <- bind_rows(
  all_stand_data[grepl("_species$", names(all_stand_data))]
) %>%
  mutate(variant = variant)

# Save
stand_file <- file.path(processed_dir, "stand_density.csv")
write_csv(stand_data, stand_file)
logger::log_info("Saved {nrow(stand_data)} stand observations to {stand_file}")

species_file <- file.path(processed_dir, "species_sdi.csv")
write_csv(species_data, species_file)
logger::log_info("Saved {nrow(species_data)} species-level SDI records to {species_file}")

# =============================================================================
# Quick Summary Statistics
# =============================================================================

cat("\n")
cat("==========================================\n")
cat("Stand-Level Data Extraction Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("States:", paste(states, collapse = ", "), "\n")
cat("Plot conditions:", nrow(stand_data), "\n")
cat("Species SDI records:", nrow(species_data), "\n")
cat("\nStand density summary:\n")
cat("  SDI range:", round(min(stand_data$sdi)), "to", round(max(stand_data$sdi)), "\n")
cat("  BA range:", round(min(stand_data$ba_ft2ac)), "to", round(max(stand_data$ba_ft2ac)), "ft2/ac\n")
cat("  QMD range:", round(min(stand_data$qmd), 1), "to", round(max(stand_data$qmd), 1), "in\n")
n_mort <- sum(!is.na(stand_data$annual_mort_rate_ba))
cat("  Plots with mortality data:", n_mort, "\n")
if (n_mort > 0) {
  cat("  Mean annual BA mortality rate:",
      round(mean(stand_data$annual_mort_rate_ba, na.rm = TRUE) * 100, 2), "%\n")
}
cat("\nOutput saved to:", processed_dir, "\n\n")

logger::log_info("Stand-level data preparation complete for variant {variant}")
