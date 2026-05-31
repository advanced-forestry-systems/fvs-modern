#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: FIA Data Fetch and Preparation
#
# Downloads or reads locally stored FIA data and prepares remeasurement pairs
# for calibrating tree level models (diameter growth, height, mortality, crown).
#
# Supports two modes:
#   1. Download via rFIA (default, for local workstation use)
#   2. Local FIA directory (for HPC with pre downloaded ENTIRE_FIA.zip)
#
# Usage:
#   Rscript calibration/R/01_fetch_fia_data.R --variant ne
#   Rscript calibration/R/01_fetch_fia_data.R --variant ne --fia-dir /path/to/fia
# =============================================================================

library(rFIA)
library(tidyverse)
library(data.table)
library(jsonlite)
library(logger)
library(sf)

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
processed_dir <- file.path(calibration_dir, "data", "processed", variant)
config_dir <- file.path(project_root, "config")

# If fia_dir not passed as arg, check environment variable
if (is.null(fia_dir)) {
  env_fia <- Sys.getenv("FVS_FIA_DATA_DIR", "")
  if (nchar(env_fia) > 0) fia_dir <- env_fia
}

dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(calibration_dir, "logs",
                      paste0("01_fetch_fia_data_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting FIA data fetch for variant {variant}")
if (!is.null(fia_dir)) {
  logger::log_info("Using local FIA directory: {fia_dir}")
}

# =============================================================================
# FVS Variant to State Mapping
#
# Used to determine which state level FIA files to load. This is a rough
# geographic filter for efficiency; the precise variant assignment comes
# from the FIA PLOT table's variant column (preferred) or falls back to
# this mapping when that column is absent in older FIA vintages.
# =============================================================================

variant_states <- list(
  ne  = c("ME", "NH", "VT", "MA", "CT", "RI", "NY", "PA", "NJ", "DE", "MD", "WV"),
  ls  = c("MI", "MN", "WI"),
  cs  = c("IN", "IL", "IA", "MO", "OH", "KY"),
  sn  = c("VA", "NC", "SC", "GA", "FL", "AL", "MS", "TN", "AR", "LA", "TX", "OK"),
  ca  = c("CA"),
  oc  = c("CA"),
  pn  = c("WA", "OR"),
  wc  = c("WA", "OR"),
  nc  = c("WA", "OR"),
  ie  = c("ID", "MT"),
  ci  = c("ID", "MT"),
  kt  = c("MT"),
  em  = c("ID", "MT", "WY", "CO", "UT", "NV", "AZ", "NM"),
  cr  = c("CO", "AZ", "NM"),
  so  = c("OR", "CA"),
  bm  = c("OR"),
  ec  = c("WA", "OR"),
  tt  = c("UT", "NV"),
  ut  = c("UT", "NV"),
  ws  = c("CA", "OR", "WA"),
  op  = c("OR"),
  ak  = c("AK"),
  acd = c("ME", "NH", "VT"),
  bc  = c("WA", "OR", "ID", "MT"),   # Proxy: Pacific Northwest US states for British Columbia
  on  = c("MN", "WI", "MI")    # Proxy: Great Lakes US states for Ontario
)

# FVS variant codes as they appear in FIA (uppercase two letter)
variant_upper <- toupper(variant)

states <- variant_states[[variant]]
if (is.null(states)) {
  logger::log_error("Unknown variant: {variant}")
  stop("Unknown variant: ", variant)
}

logger::log_info("Variant {variant} maps to states: {paste(states, collapse = ', ')}")
logger::log_info("Will filter on FIA variant column = '{variant_upper}' when available")

# =============================================================================
# Load FVS Variant Spatial Boundaries and Ecodivision Mapping
#
# The GeoJSON file contains dissolved FVS variant polygons in WGS84 for
# spatial assignment of FIA plots to variants. The ecodivision map links
# each variant to its constituent ecodivisions (Cleland et al. 2007) with
# area weights, enabling selection of ecodivision specific NSVB bark ratio
# coefficients rather than cross regional defaults.
#
# Variants not in the shapefile (ACD, BC, ON, OP, KT, OC) fall back to
# state based filtering.
# =============================================================================

variant_geojson <- file.path(calibration_dir, "data", "fvs_variant_boundaries.geojson")
ecodiv_map_file <- file.path(calibration_dir, "data", "fvs_variant_ecodivision_map.csv")

has_spatial <- file.exists(variant_geojson)
has_ecodiv_map <- file.exists(ecodiv_map_file)

if (has_spatial) {
  variant_polys <- st_read(variant_geojson, quiet = TRUE)
  logger::log_info("Loaded FVS variant boundaries: {nrow(variant_polys)} variants")
  logger::log_info("Available variants in shapefile: {paste(sort(variant_polys$FVSVariant), collapse = ', ')}")

  # Check if this variant has a spatial boundary
  variant_in_shapefile <- variant_upper %in% variant_polys$FVSVariant
  if (!variant_in_shapefile) {
    logger::log_warn("Variant {variant_upper} not in shapefile; will use state based filtering")
  }
} else {
  variant_in_shapefile <- FALSE
  logger::log_warn("Variant boundaries GeoJSON not found ({variant_geojson}); using state based filtering")
}

# Load ecodivision mapping for this variant (used in NSVB coefficient selection)
variant_ecodivisions <- NULL
if (has_ecodiv_map) {
  ecodiv_map <- read.csv(ecodiv_map_file, stringsAsFactors = FALSE)
  variant_ecodivisions <- ecodiv_map %>%
    filter(FVSVariant == variant_upper, !grepl("^Wate", DIVISION)) %>%
    arrange(desc(area_pct))

  if (nrow(variant_ecodivisions) > 0) {
    dominant_division <- variant_ecodivisions$DIVISION[1]
    logger::log_info("Variant {variant_upper} ecodivisions: {paste(variant_ecodivisions$DIVISION, collapse = ', ')}")
    logger::log_info("Dominant ecodivision: {dominant_division} ({round(variant_ecodivisions$area_pct[1], 1)}% of variant area)")
  } else {
    dominant_division <- NULL
    logger::log_warn("No ecodivision mapping found for variant {variant_upper}")
  }
} else {
  dominant_division <- NULL
  logger::log_warn("Ecodivision map not found ({ecodiv_map_file}); will use cross regional NSVB defaults")
}

# =============================================================================
# Load FVS Config for Species Mapping
# =============================================================================

config_file <- file.path(config_dir, paste0(variant, ".json"))
if (!file.exists(config_file)) {
  logger::log_error("Config file not found: {config_file}")
  stop("Config file not found: ", config_file)
}

config <- fromJSON(config_file)

# Extract species definitions
fia_species_codes <- unlist(config$categories$species_definitions$FIAJSP)
n_species <- config$maxsp
logger::log_info("Config: {n_species} species loaded for variant {variant}")

# =============================================================================
# Load NSVB Bark Ratio Tables (Westfall et al. 2024, GTR WO 104)
#
# Uses the NSVB taper equation to compute inside bark diameter at breast
# height for every FIA tree, providing a universal bark ratio (DIB/DOB)
# independent of FVS variant config. Coefficients come from the supplemental
# CSV files of the national scale volume and biomass framework.
#
# Taper equation (derivative of volume ratio model):
#   d_ib^2 = V_ib / (0.005454154 * H) * alpha * beta *
#            (1 - h/H)^(alpha - 1) * [1 - (1 - h/H)^alpha]^(beta - 1)
#
# where V_ib is gross total inside bark volume from the appropriate model
# (Schumacher Hall, segmented, continuously variable, or modified Wiley),
# and alpha/beta are volume ratio coefficients. Evaluated at h = 4.5 ft
# (breast height) this yields inside bark diameter, so bark_ratio = d_ib / D.
# =============================================================================

nsvb_dir <- file.path(calibration_dir, "data", "nsvb")

if (dir.exists(nsvb_dir)) {
  logger::log_info("Loading NSVB bark ratio tables from {nsvb_dir}")

  # Inside bark volume coefficients by SPCD (no ecodivision = cross regional default)
  volib_spcd <- tryCatch(
    read.csv(file.path(nsvb_dir, "Table S1a_volib_coefs_spcd.csv"),
             stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  # Inside bark volume ratio coefficients by SPCD
  rcumib_spcd <- tryCatch(
    read.csv(file.path(nsvb_dir, "Table S5a_rcumib_coefs_spcd.csv"),
             stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  # Jenkins group fallback tables
  volib_jenkins <- tryCatch(
    read.csv(file.path(nsvb_dir, "Table S1b_volib_coefs_jenkins.csv"),
             stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  rcumib_jenkins <- tryCatch(
    read.csv(file.path(nsvb_dir, "Table S5b_rcumib_coefs_jenkins.csv"),
             stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  # REF_SPECIES for JENKINS_SPGRPCD lookup
  ref_species <- tryCatch(
    read.csv(file.path(nsvb_dir, "REF_SPECIES.csv"),
             stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  has_nsvb <- !is.null(volib_spcd) && !is.null(rcumib_spcd)

  if (has_nsvb) {
    # ---- Ecodivision specific coefficient selection ----
    # If we have a full coefficient table with DIVISION column AND an ecodivision
    # mapping for this variant, prioritize division specific coefficients.
    # Fall back to cross regional defaults for species without division entries.
    nsvb_full_file <- file.path(nsvb_dir, "nsvb_bark_coefs_full.csv")
    use_ecodiv_coefs <- !is.null(dominant_division) && file.exists(nsvb_full_file)

    if (use_ecodiv_coefs) {
      logger::log_info("NSVB: using ecodivision specific coefficients (primary: {dominant_division})")
      nsvb_full <- read.csv(nsvb_full_file, stringsAsFactors = FALSE)

      # For each species, try the variant's ecodivisions in area weighted order,
      # then fall back to cross regional default
      div_priority <- c(variant_ecodivisions$DIVISION, "default")

      nsvb_coefs <- nsvb_full %>%
        filter(DIVISION %in% div_priority) %>%
        mutate(div_rank = match(DIVISION, div_priority)) %>%
        arrange(SPCD, div_rank) %>%
        distinct(SPCD, .keep_all = TRUE) %>%
        select(SPCD, model = model_vol, a, a1, b, b1, c, c1, alpha, beta) %>%
        mutate(across(c(a1, b1, c1), ~ replace_na(.x, NA_real_)))

      n_div_specific <- sum(nsvb_coefs$SPCD %in%
        (nsvb_full %>% filter(DIVISION == dominant_division))$SPCD)
      logger::log_info("NSVB: {n_div_specific} species with {dominant_division} specific coefs, ",
                       "{nrow(nsvb_coefs)} total species covered")
    } else {
      # Cross regional (no division) defaults
      logger::log_info("NSVB: using cross regional default coefficients")
      volib_default <- volib_spcd %>%
        filter(is.na(DIVISION) | DIVISION == "") %>%
        select(SPCD, model, a, a1, b, b1, c, c1) %>%
        distinct(SPCD, .keep_all = TRUE)

      rcumib_default <- rcumib_spcd %>%
        filter(is.na(DIVISION) | DIVISION == "") %>%
        select(SPCD, alpha, beta) %>%
        distinct(SPCD, .keep_all = TRUE)

      nsvb_coefs <- volib_default %>%
        inner_join(rcumib_default, by = "SPCD")
    }

    logger::log_info("NSVB: loaded bark ratio coefficients for {nrow(nsvb_coefs)} species")

    # Jenkins group fallback
    if (!is.null(volib_jenkins) && !is.null(rcumib_jenkins) && !is.null(ref_species)) {
      jenkins_coefs <- volib_jenkins %>%
        rename(model = model, a = a, b = b, c = c) %>%
        mutate(a1 = NA_real_, b1 = NA_real_, c1 = NA_real_) %>%
        inner_join(rcumib_jenkins %>% select(JENKINS_SPGRPCD, alpha, beta),
                   by = "JENKINS_SPGRPCD")

      jenkins_lookup <- ref_species %>%
        select(SPCD, JENKINS_SPGRPCD) %>%
        filter(SPCD <= 999)

      logger::log_info("NSVB: Jenkins fallback covers {nrow(jenkins_coefs)} groups")
    } else {
      jenkins_coefs <- NULL
      jenkins_lookup <- NULL
    }
  } else {
    logger::log_warn("NSVB tables not found; will use config bark ratios or default 1.0")
  }
} else {
  has_nsvb <- FALSE
  logger::log_warn("NSVB data directory not found ({nsvb_dir}); using fallback bark ratios")
}

# Function to compute inside bark volume for NSVB model types
# D = DBH (inches, outside bark), H = total height (feet)
nsvb_volume_ib <- function(model, SPCD, D, H, a, a1, b, b1, c_coef, c1) {
  vol <- rep(NA_real_, length(D))

  # Model 1: Schumacher Hall
  m1 <- model == 1
  vol[m1] <- a[m1] * D[m1]^b[m1] * H[m1]^c_coef[m1]

  # Model 2: Segmented
  m2 <- model == 2
  if (any(m2)) {
    k <- ifelse(SPCD[m2] < 300, 9, 11)
    small <- D[m2] < k
    vol[m2] <- ifelse(small,
      a[m2] * D[m2]^b[m2] * H[m2]^c_coef[m2],
      a[m2] * k^(b[m2] - b1[m2]) * D[m2]^b1[m2] * H[m2]^c_coef[m2]
    )
  }

  # Model 3: Continuously variable
  m3 <- model == 3
  if (any(m3)) {
    vol[m3] <- a[m3] * D[m3]^(a1[m3] * (1 - exp(-b[m3] * D[m3]))^c1[m3]) * H[m3]^c_coef[m3]
  }

  # Model 4: Modified Wiley
  m4 <- model == 4
  if (any(m4)) {
    vol[m4] <- a[m4] * D[m4]^b[m4] * H[m4]^c_coef[m4] * exp(-b1[m4] * D[m4])
  }

  return(vol)
}

# Function to compute bark ratio at breast height using NSVB taper equation
# Returns DIB / DOB for each tree
compute_nsvb_bark_ratio <- function(D, H, model, SPCD,
                                     a, a1, b, b1, c_coef, c1,
                                     alpha, beta) {
  bh <- 4.5  # breast height in feet

  # Total inside bark volume
  V_ib <- nsvb_volume_ib(model, SPCD, D, H, a, a1, b, b1, c_coef, c1)

  # Relative height
  rh <- bh / H  # h/H

  # Volume ratio derivative components at breast height
  term1 <- (1 - rh)^(alpha - 1)
  term2 <- (1 - (1 - rh)^alpha)^(beta - 1)

  # Inside bark diameter squared at breast height
  d_ib_sq <- V_ib / (0.005454154 * H) * alpha * beta * term1 * term2
  d_ib_sq <- pmax(d_ib_sq, 0)  # safety

  # Bark ratio = DIB / DOB
  d_ib <- sqrt(d_ib_sq)
  bark_ratio <- d_ib / D

  # Clamp to reasonable range (0.7 to 1.0; bark ratios outside this are suspect)
  bark_ratio <- pmin(pmax(bark_ratio, 0.70), 1.0)

  return(bark_ratio)
}

# Legacy config bark ratio as final fallback
bkrat_params <- config$categories$bark_ratio$BKRAT
has_config_bark_ratio <- !is.null(bkrat_params) && length(bkrat_params) > 0
if (has_config_bark_ratio) {
  bkrat_params <- unlist(bkrat_params)
  avg_bark_ratio <- mean(bkrat_params, na.rm = TRUE)
  logger::log_info("Config fallback bark ratio available: avg = {round(avg_bark_ratio, 3)}")
} else {
  avg_bark_ratio <- 1.0
}

# =============================================================================
# Load FIA Data (Download or Read Local)
# =============================================================================

logger::log_info("Loading FIA data for {length(states)} states...")

# Determine FIA source directory
if (!is.null(fia_dir)) {
  # Local mode: ENTIRE_FIA.zip was unzipped to fia_dir
  # Expect state subdirectories like fia_dir/ME/, fia_dir/NH/, etc.
  # OR flat structure with state prefix files like ME_TREE.csv, ME_PLOT.csv
  fia_source_dir <- fia_dir
  logger::log_info("Reading from local FIA directory: {fia_source_dir}")
} else {
  # Download mode: use rFIA to fetch data

  fia_source_dir <- file.path(calibration_dir, "data", "fia_raw")
  dir.create(fia_source_dir, showWarnings = FALSE, recursive = TRUE)
  logger::log_info("Download mode: FIA data cached in {fia_source_dir}")
}

# Process each state
all_tree_data <- list()
all_cond_data <- list()
all_plot_data <- list()

# Detect ENTIRE_FIA flat file format (ENTIRE_TREE.csv, ENTIRE_PLOT.csv, etc.)
# This format has all states in single files; rFIA::readFIA() cannot filter
# states from it reliably, so we use fread() directly with STATECD filtering.
entire_tree_file <- file.path(fia_source_dir, "ENTIRE_TREE.csv")
entire_plot_file <- file.path(fia_source_dir, "ENTIRE_PLOT.csv")
entire_cond_file <- file.path(fia_source_dir, "ENTIRE_COND.csv")
use_entire_mode <- !is.null(fia_dir) &&
                   file.exists(entire_tree_file) &&
                   file.exists(entire_plot_file) &&
                   file.exists(entire_cond_file)

if (use_entire_mode) {
  logger::log_info("Detected ENTIRE_FIA flat file format. Using fread() directly.")
}

# Build STATECD lookup from state abbreviations
state_fips <- c(
  AL=1, AK=2, AZ=4, AR=5, CA=6, CO=8, CT=9, DE=10, FL=12, GA=13,
  HI=15, ID=16, IL=17, IN=18, IA=19, KS=20, KY=21, LA=22, ME=23,
  MD=24, MA=25, MI=26, MN=27, MS=28, MO=29, MT=30, NE=31, NV=32,
  NH=33, NJ=34, NM=35, NY=36, NC=37, ND=38, OH=39, OK=40, OR=41,
  PA=42, RI=44, SC=45, SD=46, TN=47, TX=48, UT=49, VT=50, VA=51,
  WA=53, WV=54, WI=55, WY=56, AS=60, GU=66, MP=69, PR=72, VI=78
)
target_statecds <- unname(state_fips[states])
target_statecds <- target_statecds[!is.na(target_statecds)]

if (use_entire_mode) {
  # ==========================================================================
  # ENTIRE_FIA mode: read flat CSV files with fread(), filter by STATECD
  # ==========================================================================
  logger::log_info("Reading ENTIRE_TREE.csv (filtering {length(target_statecds)} states)...")

  tree_raw <- fread(entire_tree_file, select = c(
    "CN", "PLT_CN", "CONDID", "SUBP", "TREE", "INVYR", "CYCLE", "SUBCYCLE",
    "STATUSCD", "SPCD", "DIA", "HT", "CR", "ACTUALHT",
    "TPA_UNADJ", "BHAGE", "TOTAGE", "PREV_TRE_CN", "STATECD"
  ))
  tree_raw <- tree_raw[STATECD %in% target_statecds]
  logger::log_info("  Loaded {nrow(tree_raw)} tree records for target states")

  logger::log_info("Reading ENTIRE_COND.csv...")
  cond_raw <- fread(entire_cond_file, select = c(
    "PLT_CN", "CONDID", "INVYR", "STATECD",
    "FORTYPCD", "FLDTYPCD", "STDAGE", "STDSZCD", "STDORGCD",
    "SITECLCD", "SICOND", "SISP", "SLOPE", "ASPECT",
    "BALIVE", "LIVE_CANOPY_CVR_PCT"
  ))
  cond_raw <- cond_raw[STATECD %in% target_statecds]
  cond_raw[, STATECD := NULL]

  logger::log_info("Reading ENTIRE_PLOT.csv...")
  plot_raw_dt <- fread(entire_plot_file)
  plot_raw_dt <- plot_raw_dt[STATECD %in% target_statecds]

  # State code to abbreviation lookup
  fips_to_abbr <- setNames(names(state_fips), state_fips)

  # Process into the same format as the per state loop below
  for (state in states) {
    stcd <- state_fips[state]
    if (is.na(stcd)) next

    tree_df <- tree_raw[STATECD == stcd] %>%
      as_tibble() %>%
      filter(
        STATUSCD %in% c(1, 2),
        DIA >= 1.0,
        !is.na(TPA_UNADJ),
        TPA_UNADJ > 0
      ) %>%
      select(CN, PLT_CN, CONDID, SUBP, TREE, INVYR, CYCLE, SUBCYCLE,
             STATUSCD, SPCD, DIA, HT, CR, ACTUALHT,
             TPA_UNADJ, BHAGE, TOTAGE, PREV_TRE_CN) %>%
      mutate(state = state)

    cond_df <- cond_raw[cond_raw$PLT_CN %in% tree_df$PLT_CN] %>%
      as_tibble() %>%
      select(PLT_CN, CONDID, INVYR,
             FORTYPCD, FLDTYPCD, STDAGE, STDSZCD, STDORGCD,
             SITECLCD, SICOND, SISP, SLOPE, ASPECT,
             BALIVE, LIVE_CANOPY_CVR_PCT)

    # Plot data with variant column detection
    plot_state <- plot_raw_dt[STATECD == stcd] %>% as_tibble()
    plot_cols <- c("CN", "INVYR", "LAT", "LON", "ELEV", "MEASYEAR")

    var_col <- intersect(
      names(plot_state),
      c("VARIANT", "FVS_VARIANT", "FVSVARIANT", "FVSVARIANTCD",
        "variant", "fvs_variant")
    )
    if (length(var_col) > 0) {
      plot_cols <- c(plot_cols, var_col[1])
      if (state == states[1]) logger::log_info("Found FIA variant column: {var_col[1]}")
    }

    plot_df <- plot_state %>%
      select(any_of(plot_cols)) %>%
      rename(PLT_CN = CN)

    if (length(var_col) > 0) {
      plot_df <- plot_df %>% rename(FIA_VARIANT = !!sym(var_col[1]))
    }

    all_tree_data[[state]] <- tree_df
    all_cond_data[[state]] <- cond_df
    all_plot_data[[state]] <- plot_df

    logger::log_info("State {state}: {nrow(tree_df)} tree records")
  }

  # Free memory from the big data.tables
  rm(tree_raw, cond_raw, plot_raw_dt)
  gc()

} else {

# ==========================================================================
# Per state mode: use rFIA or state subdirectories
# ==========================================================================

for (state in states) {
  logger::log_info("Processing state: {state}")

  tryCatch({
    if (!is.null(fia_dir)) {
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
      db <- readFIA(dir = file.path(fia_source_dir, state), common = TRUE)
    }

    # Extract TREE table
    tree_df <- db$TREE %>%
      as_tibble() %>%
      filter(
        STATUSCD %in% c(1, 2),     # Live and recently dead
        DIA >= 1.0,                 # Minimum 1 inch DBH
        !is.na(TPA_UNADJ),
        TPA_UNADJ > 0
      ) %>%
      select(
        CN, PLT_CN, CONDID, SUBP, TREE, INVYR, CYCLE, SUBCYCLE,
        STATUSCD, SPCD, DIA, HT, CR, ACTUALHT,
        TPA_UNADJ, BHAGE, TOTAGE, PREV_TRE_CN
      ) %>%
      mutate(state = state)

    # Extract COND table
    cond_df <- db$COND %>%
      as_tibble() %>%
      select(
        PLT_CN, CONDID, INVYR,
        FORTYPCD, FLDTYPCD, STDAGE, STDSZCD, STDORGCD,
        SITECLCD, SICOND, SISP, SLOPE, ASPECT,
        BALIVE, LIVE_CANOPY_CVR_PCT
      )

    # Extract PLOT table (grab variant column if it exists)
    plot_raw <- db$PLOT %>% as_tibble()
    plot_cols <- c("CN", "INVYR", "LAT", "LON", "ELEV", "MEASYEAR")

    # FIA databases now include an FVS variant assignment per plot
    # Column names vary: VARIANT, FVS_VARIANT, FVSVARIANT, FVSVARIANTCD
    var_col <- intersect(
      names(plot_raw),
      c("VARIANT", "FVS_VARIANT", "FVSVARIANT", "FVSVARIANTCD",
        "variant", "fvs_variant")
    )
    if (length(var_col) > 0) {
      plot_cols <- c(plot_cols, var_col[1])
      logger::log_info("Found FIA variant column: {var_col[1]}")
    }

    plot_df <- plot_raw %>%
      select(any_of(plot_cols)) %>%
      rename(PLT_CN = CN)

    # Standardize variant column name to FIA_VARIANT
    if (length(var_col) > 0) {
      plot_df <- plot_df %>% rename(FIA_VARIANT = !!sym(var_col[1]))
    }

    all_tree_data[[state]] <- tree_df
    all_cond_data[[state]] <- cond_df
    all_plot_data[[state]] <- plot_df

    logger::log_info("State {state}: {nrow(tree_df)} tree records loaded")

  }, error = function(e) {
    logger::log_error("Failed to process state {state}: {e$message}")
  })
}

}  # end else (per state mode)

# Combine across states
tree_all <- bind_rows(all_tree_data)
cond_all <- bind_rows(all_cond_data)
plot_all <- bind_rows(all_plot_data)

logger::log_info("Combined: {nrow(tree_all)} tree records across {length(all_tree_data)} states")

# =============================================================================
# Join Tree + Condition + Plot
# =============================================================================

# Coerce join keys to character to avoid int64 type mismatch from fread
tree_all <- tree_all %>% mutate(PLT_CN = as.character(PLT_CN), CN = as.character(CN))
cond_all <- cond_all %>% mutate(PLT_CN = as.character(PLT_CN))
plot_all <- plot_all %>% mutate(PLT_CN = as.character(PLT_CN))

combined <- tree_all %>%
  left_join(cond_all, by = c("PLT_CN", "CONDID", "INVYR")) %>%
  left_join(plot_all, by = c("PLT_CN", "INVYR"))

logger::log_info("After join: {ncol(combined)} columns. Cond columns present: {paste(intersect(names(cond_all), names(combined)), collapse = ', ')}")

# =============================================================================
# Assign FIA Plots to FVS Variant via Spatial Join
#
# Three strategies in priority order:
#   1. Spatial join using FVS variant boundary polygons (most precise)
#   2. FIA database variant column (FIA_VARIANT) as fallback
#   3. State based mapping only (least precise, for missing variants)
#
# Variants not in the shapefile (ACD, BC, ON, OP, KT, OC) fall back to
# method 2 or 3. ACD is restricted to ME, NH, VT by state mapping.
# =============================================================================

n_before_filter <- nrow(combined)

if (variant_in_shapefile && has_spatial) {
  logger::log_info("Assigning plots to variant {variant_upper} via spatial join...")

  # Get unique plot locations
  plot_locs <- combined %>%
    distinct(PLT_CN, LAT, LON) %>%
    filter(!is.na(LAT), !is.na(LON))

  # Convert to sf points (FIA coordinates are NAD83, effectively WGS84)
  plot_sf <- st_as_sf(plot_locs, coords = c("LON", "LAT"), crs = 4326)

  # Spatial join: which variant polygon does each plot fall in?
  plot_variant <- st_join(plot_sf, variant_polys %>% select(FVSVariant), left = TRUE)
  st_geometry(plot_variant) <- NULL  # drop geometry to get a plain data frame

  # Keep plots that fall within this variant's polygon
  plots_in_variant <- plot_variant %>%
    filter(FVSVariant == variant_upper) %>%
    pull(PLT_CN) %>%
    unique()

  combined <- combined %>%
    filter(PLT_CN %in% plots_in_variant)

  logger::log_info(
    "Spatial assignment: {n_before_filter} -> {nrow(combined)} tree records for {variant_upper} ",
    "({length(plots_in_variant)} plots matched)"
  )

  if (nrow(combined) == 0) {
    logger::log_error("No plots spatially matched variant {variant_upper}")
    stop("No FIA plots spatially matched variant ", variant_upper)
  }

} else {
  # Fallback: use FIA variant column if available
  has_fia_variant <- "FIA_VARIANT" %in% names(combined)

  if (has_fia_variant) {
    combined <- combined %>%
      mutate(FIA_VARIANT = trimws(toupper(FIA_VARIANT))) %>%
      filter(FIA_VARIANT == variant_upper)

    logger::log_info(
      "Filtered by FIA variant column: {n_before_filter} -> {nrow(combined)} tree records for {variant_upper}"
    )

    if (nrow(combined) == 0) {
      logger::log_error("No plots found with FIA_VARIANT == '{variant_upper}'")
      stop("No FIA plots matched variant ", variant_upper)
    }
  } else {
    logger::log_warn(
      "No spatial boundaries or FIA variant column for {variant_upper}; using state mapping only. ",
      "This is less precise for variants that share states (e.g., IE/CI in Idaho)."
    )
  }
}

# =============================================================================
# Build Remeasurement Pairs (Same Tree, Two Inventories)
# =============================================================================

logger::log_info("Building remeasurement pairs from PREV_TRE_CN linkages...")

# Method 1: Use PREV_TRE_CN for exact tree linkage (preferred)
trees_with_prev <- combined %>%
  mutate(PREV_TRE_CN = as.character(PREV_TRE_CN)) %>%
  filter(!is.na(PREV_TRE_CN), PREV_TRE_CN != "0", PREV_TRE_CN != "NA")

# Match current tree to its previous measurement
prev_trees <- combined %>%
  select(CN, PLT_CN, INVYR, STATUSCD, DIA, HT, CR,
         TPA_UNADJ, MEASYEAR) %>%
  rename_with(~ paste0(.x, "_t1"), .cols = c(INVYR, STATUSCD, DIA, HT, CR,
                                              TPA_UNADJ, MEASYEAR))

remeas <- trees_with_prev %>%
  rename(
    INVYR_t2 = INVYR,
    STATUSCD_t2 = STATUSCD,
    DIA_t2 = DIA,
    HT_t2 = HT,
    CR_t2 = CR,
    TPA_UNADJ_t2 = TPA_UNADJ,
    MEASYEAR_t2 = MEASYEAR
  ) %>%
  left_join(
    prev_trees,
    by = c("PREV_TRE_CN" = "CN")
  ) %>%
  filter(
    !is.na(DIA_t1), !is.na(DIA_t2),
    DIA_t1 > 0, DIA_t2 > 0
  )

# Compute measurement interval
remeas <- remeas %>%
  mutate(
    years_interval = INVYR_t2 - INVYR_t1
  ) %>%
  filter(
    years_interval > 0,
    years_interval <= 15
  )

logger::log_info("Found {nrow(remeas)} remeasurement pairs via PREV_TRE_CN")

# If PREV_TRE_CN yields too few pairs, fall back to tree_id matching
if (nrow(remeas) < 100) {
  logger::log_warn("Few PREV_TRE_CN pairs; falling back to tree_id matching")

  combined <- combined %>%
    mutate(tree_id = paste(PLT_CN, SUBP, TREE, sep = "_"))

  remeas <- combined %>%
    filter(STATUSCD == 1) %>%
    group_by(tree_id) %>%
    filter(n() >= 2) %>%
    arrange(INVYR) %>%
    mutate(
      DIA_t1 = DIA,
      DIA_t2 = lead(DIA),
      HT_t1 = HT,
      HT_t2 = lead(HT),
      CR_t1 = CR,
      CR_t2 = lead(CR),
      INVYR_t1 = INVYR,
      INVYR_t2 = lead(INVYR),
      STATUSCD_t1 = STATUSCD,
      STATUSCD_t2 = lead(STATUSCD),
      years_interval = lead(INVYR) - INVYR
    ) %>%
    filter(
      !is.na(DIA_t2),
      DIA_t1 > 0, DIA_t2 > 0,
      years_interval > 0, years_interval <= 15
    ) %>%
    ungroup()

  logger::log_info("Fallback yielded {nrow(remeas)} remeasurement pairs")
}

# =============================================================================
# Compute Inside Bark Diameters and Diameter Growth
# =============================================================================

# Map FIA species codes to FVS species indices
fia_to_fvs <- setNames(seq_along(fia_species_codes), fia_species_codes)

remeas <- remeas %>%
  mutate(
    # Map species to FVS index (use 1 as fallback for unknown species)
    fvs_sp_idx = as.integer(fia_to_fvs[as.character(SPCD)]),
    fvs_sp_idx = ifelse(is.na(fvs_sp_idx), 1L, fvs_sp_idx)
  )

# ---- Compute bark ratio via NSVB taper equation (Westfall et al. 2024) ----
if (has_nsvb) {
  logger::log_info("Computing bark ratios using NSVB taper equation...")

  # Join NSVB coefficients to remeasurement data by SPCD
  remeas <- remeas %>%
    left_join(nsvb_coefs, by = "SPCD", suffix = c("", "_nsvb"))

  # For unmatched species, try Jenkins group fallback
  n_missing <- sum(is.na(remeas$model))
  if (n_missing > 0 && !is.null(jenkins_coefs) && !is.null(jenkins_lookup)) {
    logger::log_info("NSVB: {n_missing} records lack species level coefficients; trying Jenkins groups")

    remeas <- remeas %>%
      left_join(jenkins_lookup, by = "SPCD", suffix = c("", "_jk"))

    # Fill missing coefficients from Jenkins group
    missing_idx <- is.na(remeas$model) & !is.na(remeas$JENKINS_SPGRPCD)
    if (any(missing_idx)) {
      jk_fill <- jenkins_coefs[match(remeas$JENKINS_SPGRPCD[missing_idx],
                                      jenkins_coefs$JENKINS_SPGRPCD), ]
      remeas$model[missing_idx]  <- jk_fill$model
      remeas$a[missing_idx]      <- jk_fill$a
      remeas$a1[missing_idx]     <- jk_fill$a1
      remeas$b[missing_idx]      <- jk_fill$b
      remeas$b1[missing_idx]     <- jk_fill$b1
      remeas$c[missing_idx]      <- jk_fill$c
      remeas$c1[missing_idx]     <- jk_fill$c1
      remeas$alpha[missing_idx]  <- jk_fill$alpha
      remeas$beta[missing_idx]   <- jk_fill$beta

      logger::log_info("NSVB: filled {sum(!is.na(jk_fill$model))} records via Jenkins groups")
    }

    remeas <- remeas %>% select(-any_of("JENKINS_SPGRPCD"))
  }

  # Compute bark ratio for trees with valid coefficients and height
  has_coefs <- !is.na(remeas$model) & !is.na(remeas$alpha) & !is.na(remeas$beta)
  has_ht <- !is.na(remeas$HT_t1) & remeas$HT_t1 > 4.5

  can_compute <- has_coefs & has_ht
  logger::log_info("NSVB: computing bark ratio for {sum(can_compute)} of {nrow(remeas)} records")

  remeas$bark_ratio <- avg_bark_ratio  # default fallback

  if (any(can_compute)) {
    remeas$bark_ratio[can_compute] <- compute_nsvb_bark_ratio(
      D      = remeas$DIA_t1[can_compute],
      H      = remeas$HT_t1[can_compute],
      model  = remeas$model[can_compute],
      SPCD   = remeas$SPCD[can_compute],
      a      = remeas$a[can_compute],
      a1     = remeas$a1[can_compute],
      b      = remeas$b[can_compute],
      b1     = remeas$b1[can_compute],
      c_coef = remeas$c[can_compute],
      c1     = remeas$c1[can_compute],
      alpha  = remeas$alpha[can_compute],
      beta   = remeas$beta[can_compute]
    )
  }

  # Log summary
  br_vals <- remeas$bark_ratio[can_compute]
  logger::log_info("NSVB bark ratio: mean = {round(mean(br_vals), 4)}, ",
                   "range = [{round(min(br_vals), 4)}, {round(max(br_vals), 4)}]")

  # Clean up NSVB columns
  remeas <- remeas %>%
    select(-any_of(c("model", "a", "a1", "b", "b1", "c", "c1", "alpha", "beta")))

} else {
  # Fallback: use config bark ratios or 1.0
  logger::log_info("Using config/default bark ratios (NSVB tables not available)")

  remeas <- remeas %>%
    mutate(
      bark_ratio = if (has_config_bark_ratio) {
        ifelse(
          fvs_sp_idx <= length(bkrat_params) &
            !is.na(bkrat_params[pmin(fvs_sp_idx, length(bkrat_params))]) &
            bkrat_params[pmin(fvs_sp_idx, length(bkrat_params))] > 0,
          bkrat_params[pmin(fvs_sp_idx, length(bkrat_params))],
          avg_bark_ratio
        )
      } else {
        1.0
      }
    )
}

# Compute derived diameter variables
remeas <- remeas %>%
  mutate(
    # Inside bark diameters
    DIB_t1 = DIA_t1 * bark_ratio,
    DIB_t2 = DIA_t2 * bark_ratio,

    # Squared diameter increment (FVS Wykoff model predicts ln(DDS))
    DDS = DIB_t2^2 - DIB_t1^2,

    # Annual height increment
    HT_annual = ifelse(
      !is.na(HT_t1) & !is.na(HT_t2) & HT_t2 > HT_t1,
      (HT_t2 - HT_t1) / years_interval,
      NA_real_
    ),

    # Annual crown ratio change
    CR_annual = ifelse(
      !is.na(CR_t1) & !is.na(CR_t2),
      (CR_t2 - CR_t1) / years_interval,
      NA_real_
    )
  )

# =============================================================================
# Compute Stand Level Covariates (BA, BAL, SDI, etc.)
# =============================================================================

logger::log_info("Computing stand level covariates...")

# Stand metrics from live trees at time 1
stand_metrics <- combined %>%
  filter(STATUSCD == 1, !is.na(DIA), DIA > 0) %>%
  group_by(PLT_CN, CONDID, INVYR) %>%
  summarise(
    BA_ft2ac = sum(TPA_UNADJ * 0.005454 * DIA^2, na.rm = TRUE),
    TPA = sum(TPA_UNADJ, na.rm = TRUE),
    QMD = sqrt(BA_ft2ac / TPA / 0.005454),
    SDI = TPA * (QMD / 10)^1.605,
    top_ht = {
      ht_sorted <- sort(HT[!is.na(HT)], decreasing = TRUE)
      if (length(ht_sorted) >= 5) mean(ht_sorted[1:min(5, length(ht_sorted))])
      else if (length(ht_sorted) > 0) mean(ht_sorted)
      else NA_real_
    },
    .groups = "drop"
  )

# Basal Area in Larger trees (BAL) for each tree
bal_data <- combined %>%
  filter(STATUSCD == 1, !is.na(DIA), DIA > 0) %>%
  group_by(PLT_CN, CONDID, INVYR) %>%
  arrange(desc(DIA)) %>%
  mutate(
    tree_ba = TPA_UNADJ * 0.005454 * DIA^2,
    BAL_ft2ac = cumsum(tree_ba) - tree_ba  # BA of trees larger than this one
  ) %>%
  select(PLT_CN, CONDID, INVYR, SUBP, TREE, SPCD, BAL_ft2ac) %>%
  ungroup()

# Site index: use SICOND from FIA condition table where available
site_index <- cond_all %>%
  filter(!is.na(SICOND), SICOND > 0) %>%
  distinct(PLT_CN, CONDID, INVYR, SICOND)

# =============================================================================
# Join Covariates to Remeasurement Pairs
# =============================================================================

# Use time 1 inventory year for covariate matching
remeas <- remeas %>%
  left_join(
    stand_metrics %>% rename(BA = BA_ft2ac),
    by = c("PLT_CN.y" = "PLT_CN", "CONDID" = "CONDID", "INVYR_t1" = "INVYR")
  ) %>%
  left_join(
    bal_data,
    by = c("PLT_CN.x" = "PLT_CN", "CONDID" = "CONDID",
           "INVYR_t2" = "INVYR", "SUBP" = "SUBP", "TREE" = "TREE", "SPCD" = "SPCD")
  ) %>%
  mutate(
    # Use SICOND as site index (from condition table); fall back to top_ht
    SI = coalesce(SICOND, top_ht),
    # Use computed BAL; fall back to BALIVE from condition table
    BAL = coalesce(BAL_ft2ac, BALIVE)
  )

# =============================================================================
# Compute Derived Model Variables
# =============================================================================

remeas <- remeas %>%
  mutate(
    # Slope / aspect transformations
    ASPECT_rad = ASPECT * pi / 180,
    SLOPE_SASP = SLOPE * sin(ASPECT_rad),
    SLOPE_CASP = SLOPE * cos(ASPECT_rad),
    ELEV_100ft = ELEV / 100,

    # Logarithmic transforms for Wykoff model
    ln_DBH = log(DIA_t1),
    DBH_sq = DIA_t1^2,
    ln_SI = log(pmax(SI, 1)),
    CR_pct = CR_t1 * 100,  # FIA stores 0 to 1, some models want 0 to 100

    # Annual DDS
    DDS_annual = DDS / years_interval,
    ln_DDS = log(pmax(DDS, 0.001))
  )

# =============================================================================
# Quality Filters
# =============================================================================

logger::log_info("Applying quality filters...")

n_before <- nrow(remeas)

remeas_clean <- remeas %>%
  filter(
    !is.na(DDS), DDS > 0,
    !is.na(DIA_t1), DIA_t1 >= 1.0,
    !is.na(SI), SI > 0,
    !is.na(BA), BA > 0,
    !is.na(SLOPE),
    !is.na(ELEV)
  )

logger::log_info("Quality filter: {n_before} -> {nrow(remeas_clean)} observations ({round(100 * nrow(remeas_clean)/n_before, 1)}% retained)")

# =============================================================================
# Save Diameter Growth Dataset
# =============================================================================

dg_out <- remeas_clean %>%
  transmute(
    PLT_CN = coalesce(PLT_CN.x, PLT_CN.y),
    tree_id = paste(PLT_CN, SUBP, TREE, sep = "_"),
    SPCD, fvs_sp_idx, variant = !!variant,
    # Measurements
    DIA_t1, DIA_t2, DIB_t1, DIB_t2, DDS, ln_DDS,
    HT_t1, HT_t2, HT_annual,
    CR_t1, CR_t2, CR_annual,
    years_interval,
    STATUSCD_t1, STATUSCD_t2,
    # Covariates
    ln_DBH, DBH_sq, ln_SI, SI,
    SLOPE, SLOPE_SASP, SLOPE_CASP,
    ELEV, ELEV_100ft,
    CR_pct, BAL, BA, SDI, QMD, TPA,
    # Location
    LAT, LON
  )

dg_file <- file.path(processed_dir, "diameter_growth.csv")
write_csv(dg_out, dg_file)
logger::log_info("Saved diameter growth data: {nrow(dg_out)} obs to {dg_file}")

# =============================================================================
# Save Height Growth Dataset (trees with valid height remeasurements)
# =============================================================================

ht_out <- remeas_clean %>%
  filter(!is.na(HT_t1), !is.na(HT_t2), HT_t1 > 4.5, HT_t2 > HT_t1) %>%
  transmute(
    PLT_CN = coalesce(PLT_CN.x, PLT_CN.y),
    tree_id = paste(PLT_CN, SUBP, TREE, sep = "_"),
    SPCD, fvs_sp_idx, variant = !!variant,
    DIA_t1, HT_t1, HT_t2, HT_annual,
    years_interval,
    ln_DBH, SI, BA, BAL, CR_pct,
    SLOPE, ELEV, LAT, LON
  )

ht_file <- file.path(processed_dir, "height_growth.csv")
write_csv(ht_out, ht_file)
logger::log_info("Saved height growth data: {nrow(ht_out)} obs to {ht_file}")

# =============================================================================
# Save Height Diameter Dataset (all live trees with HT and DIA)
# =============================================================================

hd_out <- combined %>%
  filter(STATUSCD == 1, !is.na(DIA), DIA >= 1.0, !is.na(HT), HT > 4.5) %>%
  transmute(
    PLT_CN, CONDID, INVYR,
    SPCD,
    fvs_sp_idx = as.integer(fia_to_fvs[as.character(SPCD)]),
    fvs_sp_idx = ifelse(is.na(fvs_sp_idx), 1L, fvs_sp_idx),
    variant = !!variant,
    DIA, HT,
    SI = coalesce(SICOND, 50),  # Default SI = 50 if missing (SICOND from cond join)
    LAT, LON, ELEV
  )

hd_file <- file.path(processed_dir, "height_diameter.csv")
write_csv(hd_out, hd_file)
logger::log_info("Saved height diameter data: {nrow(hd_out)} obs to {hd_file}")

# =============================================================================
# Save Mortality Dataset (live trees + trees that died between measurements)
# =============================================================================

mort_out <- remeas %>%
  filter(
    STATUSCD_t1 == 1,          # Alive at time 1
    !is.na(DIA_t1), DIA_t1 >= 1.0,
    !is.na(SI), SI > 0,
    !is.na(BA), BA > 0,
    years_interval > 0
  ) %>%
  transmute(
    PLT_CN = coalesce(PLT_CN.x, PLT_CN.y),
    tree_id = paste(PLT_CN, SUBP, TREE, sep = "_"),
    SPCD, fvs_sp_idx, variant = !!variant,
    DIA = DIA_t1, HT = HT_t1, CR = CR_t1,
    died = as.integer(STATUSCD_t2 == 2),  # 1 = died, 0 = survived
    years_interval,
    ln_DBH = log(DIA_t1), SI, BA, BAL,
    CR_pct = CR_t1 * 100,
    SLOPE, ELEV, LAT, LON
  )

mort_file <- file.path(processed_dir, "mortality.csv")
write_csv(mort_out, mort_file)
logger::log_info("Saved mortality data: {nrow(mort_out)} obs ({sum(mort_out$died)} deaths)")

# =============================================================================
# Save Crown Ratio Change Dataset
# =============================================================================

cr_out <- remeas_clean %>%
  filter(!is.na(CR_t1), !is.na(CR_t2), CR_t1 > 0) %>%
  transmute(
    PLT_CN = coalesce(PLT_CN.x, PLT_CN.y),
    tree_id = paste(PLT_CN, SUBP, TREE, sep = "_"),
    SPCD, fvs_sp_idx, variant = !!variant,
    DIA = DIA_t1, DIA_sq = DIA_t1^2,
    CR_init = CR_t1, CR_final = CR_t2,
    delta_CR = CR_annual,
    years_interval,
    SI, BA, BAL,
    SLOPE, ELEV, LAT, LON
  )

cr_file <- file.path(processed_dir, "crown_ratio_change.csv")
write_csv(cr_out, cr_file)
logger::log_info("Saved crown ratio change data: {nrow(cr_out)} obs to {cr_file}")

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat("==========================================\n")
cat("FIA Data Preparation Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("States:", paste(states, collapse = ", "), "\n")
cat("FIA source:", ifelse(!is.null(fia_dir), fia_dir, "rFIA download"), "\n")
cat("\nDatasets created:\n")
cat("  Diameter growth:", nrow(dg_out), "observations\n")
cat("  Height growth:  ", nrow(ht_out), "observations\n")
cat("  Height diameter:", nrow(hd_out), "observations\n")
cat("  Mortality:      ", nrow(mort_out), "observations",
    "(", sum(mort_out$died), "deaths )\n")
cat("  Crown ratio:    ", nrow(cr_out), "observations\n")
cat("\nOutput directory:", processed_dir, "\n\n")

logger::log_info("FIA data preparation complete for variant {variant}")
