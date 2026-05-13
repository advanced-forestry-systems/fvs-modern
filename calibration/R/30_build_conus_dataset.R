## =============================================================================
## 30_build_conus_dataset.R
##
## FVS-CONUS: Build unified remeasurement dataset from all state zip files.
##
## Reads CHANGEdata CSVs from David Marshall's compiled state zips, applies
## filtering rules (April 14, 2026 instructions), assigns EPA Level 3
## ecoregion (replacing Bailey's ecodivision), joins site productivity
## indices (asymptotic biomass, BGI, ClimateSI, composite CSPI), and
## produces the unified input RDS consumed by scripts 31-34.
##
## Filtering logic (per D. Marshall, 2026-04-14):
##   DG:        startSTATUSCD=1, endSTATUSCD=1, DIA not NA at both ends
##   HG:        same + startHTCD=1, endHTCD=1 (field measured heights)
##   Mortality: startSTATUSCD=1, endSTATUSCD=2, endAGENTCD<80 (no harvest)
##   Plot:      COND_STATUS_CD=1, CONDPROP_UNADJ>=0.75
##
## Output:
##   calibration/data/conus_remeasurement_pairs.rds
##   calibration/data/conus_hcb_cross_section.rds     (for script 33)
##   calibration/data/conus_data_summary.csv
##
## Usage on Cardinal:
##   Rscript calibration/R/30_build_conus_dataset.R
##   Rscript calibration/R/30_build_conus_dataset.R --zip_dir calibration/data/fia_zips
##
## Author: A. Weiskittel, G. Johnson
## Date: 2026-04-14
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(sf)
  library(glue)
})

## ---- Configuration ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

## Defaults (Cardinal paths)
ZIP_DIR       <- "calibration/data/fia_zips"
OUT_DIR       <- "calibration/data"
BGI_RDS       <- "~/FIA/asym_agb_analysis/results/stage6_bgi_results.rds"
BLUP_CSV      <- "~/FIA/asym_agb_analysis/results/blup_plot_asymptotes.csv"
SDIMAX_CSV    <- "~/SDImax/brms.SDImax.1-27-24.csv"
RASTER_CSV    <- "calibration/data/plot_raster_lookup.csv"
EPA_L3_SHP    <- "~/GIS/us_eco_l3.shp"
CSPI_RDS      <- "~/FIA/asym_agb_analysis/results/stage9_composite_productivity.rds"
BLUP_RDS      <- "~/FIA/asym_agb_analysis/results/blup_plot_asymptotes.rds"
CLIM_PCA_TIF  <- "~/SiteIndex/CONUS_Climate_Embedding_Top4.tif"
FIA_COND_CSV  <- "~/FIA/ENTIRE_COND.csv"   # for ECOSUBCD (NSBE anchor)

## Legacy Bailey's ecodivision (retained as fallback, but EPA L3 preferred)
VARIANT_GEOJSON <- "calibration/data/fvs_variant_boundaries.geojson"
ECODIV_MAP    <- "calibration/data/fvs_variant_ecodivision_map.csv"

## Parse CLI overrides
for (i in seq_along(args)) {
  if (args[i] == "--zip_dir"    && i < length(args)) ZIP_DIR       <- args[i + 1]
  if (args[i] == "--out_dir"    && i < length(args)) OUT_DIR       <- args[i + 1]
  if (args[i] == "--bgi_rds"    && i < length(args)) BGI_RDS       <- args[i + 1]
  if (args[i] == "--blup_csv"   && i < length(args)) BLUP_CSV      <- args[i + 1]
  if (args[i] == "--sdimax_csv" && i < length(args)) SDIMAX_CSV    <- args[i + 1]
  if (args[i] == "--raster_csv" && i < length(args)) RASTER_CSV    <- args[i + 1]
  if (args[i] == "--epa_l3"     && i < length(args)) EPA_L3_SHP    <- args[i + 1]
  if (args[i] == "--cspi_rds"   && i < length(args)) CSPI_RDS      <- args[i + 1]
  if (args[i] == "--clim_pca"   && i < length(args)) CLIM_PCA_TIF  <- args[i + 1]
  if (args[i] == "--cond_csv"   && i < length(args)) FIA_COND_CSV  <- args[i + 1]
}

CONDPROP_MIN  <- 0.75    # minimum condition proportion (Feb 2026 decision)
MIN_REMPER    <- 1       # minimum remeasurement period (years)
MAX_REMPER    <- 20      # maximum remeasurement period (years)
EXCLUDE_AK    <- TRUE    # exclude AK from CONUS fitting

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("FVS-CONUS: Build Unified Remeasurement Dataset\n")
cat("==============================================================\n")
cat("Zip directory:   ", ZIP_DIR, "\n")
cat("Output directory: ", OUT_DIR, "\n")
cat("CONDPROP minimum: ", CONDPROP_MIN, "\n")
cat("Start time:       ", format(Sys.time()), "\n")
cat("==============================================================\n\n")

## ---- Step 1: Read all CHANGEdata files from state zips --------------------

cat("Step 1: Reading CHANGEdata from state zip files...\n")

zip_files <- list.files(ZIP_DIR, pattern = "_FIAdata_\\d{8}\\.zip$",
                        full.names = TRUE)
cat("  Found", length(zip_files), "zip files\n")

if (length(zip_files) == 0) stop("No zip files found in ", ZIP_DIR)

## Function to read CHANGEdata from a single zip
read_change_from_zip <- function(zip_path) {
  st <- substr(basename(zip_path), 1, 2)

  ## Skip AK if configured
  if (EXCLUDE_AK && st == "AK") {
    cat("  ", st, ": skipping (non-CONUS)\n")
    return(NULL)
  }

  ## List files in the zip and find CHANGEdata
  contents <- tryCatch(unzip(zip_path, list = TRUE)$Name,
                       error = function(e) character())
  change_file <- contents[grepl("CHANGEdata", contents, ignore.case = TRUE)]

  if (length(change_file) == 0) {
    cat("  ", st, ": WARNING - no CHANGEdata found, skipping\n")
    return(NULL)
  }

  ## Extract to temp directory and read
  tmp_dir <- tempdir()
  unzip(zip_path, files = change_file[1], exdir = tmp_dir, overwrite = TRUE)
  extracted <- file.path(tmp_dir, change_file[1])

  dt <- tryCatch(
    fread(extracted, showProgress = FALSE, na.strings = c("NA", "")),
    error = function(e) {
      cat("  ", st, ": ERROR reading - ", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  ## Add state abbreviation
  dt[, STATE_ABBR := st]
  nrow_raw <- nrow(dt)

  cat("  ", st, ":", format(nrow_raw, big.mark = ","), "records\n")
  unlink(extracted)

  return(dt)
}

## Read all states
change_list <- map(zip_files, read_change_from_zip)
change_list <- compact(change_list)

cat("\n  Combining", length(change_list), "states...\n")
change_all <- rbindlist(change_list, use.names = TRUE, fill = TRUE)
rm(change_list); gc(verbose = FALSE)

cat("  Total raw records:", format(nrow(change_all), big.mark = ","), "\n")
cat("  States represented:", uniqueN(change_all$STATE_ABBR), "\n\n")

## ---- Step 2: Plot-level filtering (COND_STATUS, CONDPROP) ----------------

cat("Step 2: Applying plot-level filters...\n")

n_before <- nrow(change_all)

## Filter: accessible forest (COND_STATUS_CD = 1) and sufficient
## condition proportion (>= 0.75 per Feb 2026 decision)
change_all <- change_all[
  !is.na(start.COND_STATUS_CD) & start.COND_STATUS_CD == 1 &
  !is.na(startCONDPROP_UNADJ) & startCONDPROP_UNADJ >= CONDPROP_MIN
]

cat("  COND_STATUS_CD = 1 & CONDPROP >= ", CONDPROP_MIN, ": ",
    format(nrow(change_all), big.mark = ","), " records",
    " (removed ", format(n_before - nrow(change_all), big.mark = ","), ")\n")

## Filter: valid remeasurement period
n_before <- nrow(change_all)
change_all <- change_all[
  !is.na(REMPER) & REMPER >= MIN_REMPER & REMPER <= MAX_REMPER
]
cat("  REMPER ", MIN_REMPER, "-", MAX_REMPER, " years: ",
    format(nrow(change_all), big.mark = ","), " records\n")

## Remove ingrowth records (startTreeStatus = 4) from the modeling dataset
## These have NA start-period attributes and exist only for end-of-period
## stand summary computation (already done by David)
n_before <- nrow(change_all)
change_all <- change_all[is.na(startTreeStatus) | startTreeStatus != 4]
cat("  Remove ingrowth (startTreeStatus=4): ",
    format(nrow(change_all), big.mark = ","), " records",
    " (removed ", format(n_before - nrow(change_all), big.mark = ","),
    " ingrowth)\n\n")


## ---- Step 3: EPA Level 3 Ecoregion Assignment ----------------------------

cat("Step 3: Assigning EPA Level 3 ecoregions...\n")

## Create plot-level coordinates (unique plots)
plot_coords <- change_all[
  !is.na(LAT) & !is.na(LON),
  .(LAT = LAT[1], LON = LON[1]),
  by = .(STATECD, UNITCD, COUNTYCD, PLOT)
]

## Uses NA_L3CODE (hierarchical format like "8.3.5") consistently in both
## the BLUP pre-assigned codes and the spatial-join fallback. This gives 84
## unique codes for CONUS matching the EPA Level III Ecoregions shapefile.

## 3a. Primary method: BLUP data already has EPA_L3_CODE (in NA_L3CODE format)
##     for ~135K plots (from Aaron's asymptotic biomass analysis)
blup_eco <- NULL
if (file.exists(BLUP_CSV)) {
  cat("  Attempting EPA L3 assignment from BLUP data (pre-assigned)...\n")
  blup_raw <- fread(BLUP_CSV, showProgress = FALSE,
                    select = c("PlotID", "EPA_L3_CODE"))
  blup_raw[, plot_key := gsub("_", "-", PlotID)]
  ## Treat "UNKNOWN" as NA so the spatial-join fallback can fix them
  blup_raw[EPA_L3_CODE == "UNKNOWN", EPA_L3_CODE := NA_character_]
  plot_coords[, plot_key := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "-")]
  blup_eco <- merge(plot_coords[, .(plot_key)],
                    blup_raw[, .(plot_key, EPA_L3_CODE)],
                    by = "plot_key", all.x = TRUE)
  n_from_blup <- sum(!is.na(blup_eco$EPA_L3_CODE))
  cat("    Matched", format(n_from_blup, big.mark = ","), "of",
      format(nrow(plot_coords), big.mark = ","), "plots from BLUP\n")
  plot_coords <- merge(plot_coords, blup_eco, by = "plot_key", all.x = TRUE)
}

## 3b. For remaining plots without EPA L3: spatial join from shapefile
##     using NA_L3CODE (not US_L3CODE) for consistency with BLUP format
n_missing_eco <- sum(is.na(plot_coords$EPA_L3_CODE))

if (n_missing_eco > 0 && file.exists(EPA_L3_SHP)) {
  cat("  Spatial join for", format(n_missing_eco, big.mark = ","),
      "plots without EPA L3 assignment...\n")

  ## Load EPA L3 shapefile (Lambert Azimuthal Equal Area)
  epa_sf <- st_read(EPA_L3_SHP, quiet = TRUE)

  ## Get plots needing assignment
  needs_eco <- plot_coords[is.na(EPA_L3_CODE)]
  pts_sf <- st_as_sf(needs_eco, coords = c("LON", "LAT"), crs = 4326)

  ## Transform plot points to match shapefile CRS
  pts_sf <- st_transform(pts_sf, st_crs(epa_sf))

  ## Point-in-polygon spatial join; use NA_L3CODE for CONUS-wide hierarchy
  joined <- st_join(pts_sf, epa_sf[, c("NA_L3CODE", "NA_L3NAME")],
                    join = st_within)

  ## For any points still unmatched, use nearest feature
  still_missing <- which(is.na(joined$NA_L3CODE))
  if (length(still_missing) > 0) {
    cat("    ", length(still_missing),
        " plots outside all polygons, assigning nearest\n")
    nearest_idx <- st_nearest_feature(pts_sf[still_missing, ], epa_sf)
    joined$NA_L3CODE[still_missing] <- epa_sf$NA_L3CODE[nearest_idx]
    joined$NA_L3NAME[still_missing] <- epa_sf$NA_L3NAME[nearest_idx]
  }

  ## Update plot_coords
  needs_eco$EPA_L3_CODE <- as.character(joined$NA_L3CODE)
  needs_eco$EPA_L3_NAME <- as.character(joined$NA_L3NAME)
  plot_coords[is.na(EPA_L3_CODE),
              c("EPA_L3_CODE", "EPA_L3_NAME") := needs_eco[, .(EPA_L3_CODE, EPA_L3_NAME)]]

  cat("    EPA L3 spatial join complete\n")

} else if (n_missing_eco > 0) {
  cat("  WARNING: EPA L3 shapefile not found at ", EPA_L3_SHP, "\n")
  cat("  ", n_missing_eco, " plots remain without ecoregion assignment\n")
}

## 3c. Also assign FVS variant via variant boundaries (for reference/comparison)
if (file.exists(VARIANT_GEOJSON)) {
  cat("  Assigning FVS variant labels (reference only)...\n")
  variant_sf <- st_read(VARIANT_GEOJSON, quiet = TRUE)
  plot_pts <- st_as_sf(plot_coords, coords = c("LON", "LAT"), crs = 4326)
  plot_variant <- st_join(plot_pts, variant_sf, join = st_within)
  plot_coords$fvs_variant <- plot_variant$FVSVariant

  ## Nearest for unmatched
  missing_idx <- which(is.na(plot_coords$fvs_variant))
  if (length(missing_idx) > 0) {
    nearest_idx <- st_nearest_feature(plot_pts[missing_idx, ], variant_sf)
    plot_coords$fvs_variant[missing_idx] <- variant_sf$FVSVariant[nearest_idx]
  }
} else {
  plot_coords[, fvs_variant := NA_character_]
}

## Rename EPA L3 to the standard column name used by scripts 31-34
## Scripts reference "ecodiv_code" so we map EPA_L3_CODE to that field.
## This is the key change: ecodiv_code is now EPA Level 3, not Bailey's division.
plot_coords[, ecodiv_code := EPA_L3_CODE]

## Join back to main data
join_cols <- intersect(c("STATECD", "UNITCD", "COUNTYCD", "PLOT",
                         "fvs_variant", "ecodiv_code", "EPA_L3_CODE"),
                       names(plot_coords))
change_all <- merge(change_all,
                    plot_coords[, ..join_cols],
                    by = c("STATECD", "UNITCD", "COUNTYCD", "PLOT"),
                    all.x = TRUE)

cat("  EPA L3 ecoregions assigned: ",
    uniqueN(change_all$ecodiv_code[!is.na(change_all$ecodiv_code)]),
    " unique\n")
cat("  FVS variants assigned: ",
    uniqueN(change_all$fvs_variant[!is.na(change_all$fvs_variant)]),
    " unique\n")

## Fallback for any remaining NAs (e.g. plots with no LAT/LON)
if (any(is.na(change_all$ecodiv_code))) {
  n_na <- sum(is.na(change_all$ecodiv_code))
  cat("  WARNING:", n_na, "records without ecoregion (missing LAT/LON?)\n")
  state_mode <- change_all[!is.na(ecodiv_code),
                           .(mode_eco = names(sort(table(ecodiv_code),
                                                   decreasing = TRUE))[1]),
                           by = STATECD]
  change_all <- merge(change_all, state_mode, by = "STATECD", all.x = TRUE)
  change_all[is.na(ecodiv_code), ecodiv_code := mode_eco]
  change_all[, mode_eco := NULL]
}

cat("\n")


## ---- Step 3d: Join FIA ECOSUBCD from ENTIRE_COND.csv ---------------------
## ECOSUBCD is FIA's Bailey-based ecoregion subsection. Retained for
## cross-comparison with the EPA L3 grouping (the NSBE equations are
## anchored to ECOSUBCD, so we want both in the dataset to evaluate whether
## EPA L3 substantively changes the species x ecoregion hierarchy.)

cat("Step 3d: Joining FIA ECOSUBCD from COND table...\n")

if (file.exists(FIA_COND_CSV)) {
  cond_sel <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT", "CONDID",
                "INVYR", "ECOSUBCD", "COND_STATUS_CD", "CONDPROP_UNADJ")
  cond_dt <- fread(FIA_COND_CSV, showProgress = FALSE, select = cond_sel)

  ## Take the dominant (largest-CONDPROP) forested condition per plot
  cond_dt <- cond_dt[COND_STATUS_CD == 1 & !is.na(ECOSUBCD)]
  setorder(cond_dt, STATECD, UNITCD, COUNTYCD, PLOT, -CONDPROP_UNADJ)
  cond_u <- cond_dt[, .SD[1], by = .(STATECD, UNITCD, COUNTYCD, PLOT)]

  change_all <- merge(change_all,
                      cond_u[, .(STATECD, UNITCD, COUNTYCD, PLOT, ECOSUBCD)],
                      by = c("STATECD", "UNITCD", "COUNTYCD", "PLOT"),
                      all.x = TRUE)
  n_ecosub <- sum(!is.na(change_all$ECOSUBCD))
  cat("    ECOSUBCD joined:", format(n_ecosub, big.mark = ","),
      "of", format(nrow(change_all), big.mark = ","), "records\n")
  cat("    Unique ECOSUBCDs:", uniqueN(change_all$ECOSUBCD), "\n")
} else {
  cat("  WARNING: ENTIRE_COND.csv not found at", FIA_COND_CSV, "\n")
  change_all[, ECOSUBCD := NA_character_]
}
cat("\n")


## ---- Step 4: Join Site Productivity Indices ------------------------------

cat("Step 4: Joining site productivity indices...\n")

## Create plot key for joins
change_all[, plot_key := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "-")]

## 4a. ClimateSI from raster lookup
if (file.exists(RASTER_CSV)) {
  cat("  Reading ClimateSI raster lookup...\n")
  raster_lkp <- fread(RASTER_CSV, showProgress = FALSE)

  ## The raster lookup uses PLT_CN which we don't have directly.
  ## Join via LAT/LON (rounded) as a spatial proxy.
  ## This is approximate but adequate for the fuzzed FIA coordinates.
  raster_lkp[, LAT_r := round(LAT, 2)]
  raster_lkp[, LON_r := round(LON, 2)]
  raster_agg <- raster_lkp[, .(
    ClimateSI_ft = mean(ClimateSI_ft, na.rm = TRUE),
    SDIMAX_raster = mean(as.numeric(SDIMAX_imperial), na.rm = TRUE)
  ), by = .(LAT_r, LON_r)]

  change_all[, LAT_r := round(LAT, 2)]
  change_all[, LON_r := round(LON, 2)]
  change_all <- merge(change_all, raster_agg,
                      by = c("LAT_r", "LON_r"), all.x = TRUE)
  change_all[, c("LAT_r", "LON_r") := NULL]

  ## Convert ClimateSI from feet to meters for modeling
  change_all[, climate_si := ClimateSI_ft * 0.3048]
  cat("    ClimateSI joined:", sum(!is.na(change_all$climate_si)),
      "of", nrow(change_all), "records\n")
} else {
  cat("  WARNING: Raster lookup not found, skipping ClimateSI\n")
  change_all[, climate_si := NA_real_]
}

## 4b. BGI from stage6 results
if (file.exists(BGI_RDS)) {
  cat("  Reading BGI data...\n")
  bgi_raw <- readRDS(BGI_RDS)

  ## BGI is in $plot_bgi; join via rounded LAT/LON
  if ("plot_bgi" %in% names(bgi_raw)) {
    bgi_dt <- as.data.table(bgi_raw$plot_bgi)
  } else {
    bgi_dt <- as.data.table(bgi_raw)
  }

  if ("BGI_A_mean" %in% names(bgi_dt)) {
    bgi_dt[, LAT_r := round(LAT, 2)]
    bgi_dt[, LON_r := round(LON, 2)]
    bgi_agg <- bgi_dt[, .(bgi = mean(BGI_A_mean, na.rm = TRUE)),
                       by = .(LAT_r, LON_r)]

    change_all[, LAT_r := round(LAT, 2)]
    change_all[, LON_r := round(LON, 2)]
    change_all <- merge(change_all, bgi_agg,
                        by = c("LAT_r", "LON_r"), all.x = TRUE)
    change_all[, c("LAT_r", "LON_r") := NULL]

    cat("    BGI joined:", sum(!is.na(change_all$bgi)),
        "of", nrow(change_all), "records\n")
  }
} else {
  cat("  WARNING: BGI RDS not found, skipping\n")
  change_all[, bgi := NA_real_]
}

## 4c. BLUP asymptotes (max_biomass = Asym_mean, byi = k)
if (file.exists(BLUP_CSV)) {
  cat("  Reading BLUP asymptotes...\n")
  blup_dt <- fread(BLUP_CSV, showProgress = FALSE)

  ## PlotID format: "STATECD_UNITCD_COUNTYCD_PLOT"
  if ("PlotID" %in% names(blup_dt)) {
    blup_dt[, plot_key := gsub("_", "-", PlotID)]
    blup_join <- blup_dt[, .(plot_key,
                              max_biomass = Asym_mean,
                              byi = k)]
    change_all <- merge(change_all, blup_join,
                        by = "plot_key", all.x = TRUE)
    cat("    BLUP joined:", sum(!is.na(change_all$max_biomass)),
        "of", nrow(change_all), "records\n")
  }
} else {
  cat("  WARNING: BLUP CSV not found, skipping\n")
  change_all[, c("max_biomass", "byi") := NA_real_]
}

## 4d. SDImax (brms) for RD_sdimax computation
if (file.exists(SDIMAX_CSV)) {
  cat("  Reading SDImax (brms)...\n")
  sdi_dt <- fread(SDIMAX_CSV, showProgress = FALSE)
  sdi_dt[, plot_key := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "-")]
  sdi_join <- sdi_dt[, .(plot_key, SDImax_brms = SDImax.mean)]

  change_all <- merge(change_all, sdi_join,
                      by = "plot_key", all.x = TRUE)
  cat("    SDImax joined:", sum(!is.na(change_all$SDImax_brms)),
      "of", nrow(change_all), "records\n")
} else {
  cat("  WARNING: SDImax CSV not found, skipping\n")
  change_all[, SDImax_brms := NA_real_]
}

## 4e. Composite Site Productivity Index (CSPI) from stage9 analysis
##     Includes CSPI_PCA, CSPI_GEOM, and component metrics.
##     The composite_data in stage9 has 89K rows but no plot keys,
##     so we compute our own composite from the per-plot BLUP + BGI + ClimateSI
##     that are already joined above.
cat("  Computing composite site productivity (CSPI)...\n")

## Method: z-score normalization of available metrics, then combine via
## PCA loading weights from stage9 analysis: BGI 0.40, Asym 0.35, SI 0.25
## (weights from Aaron's stage9_composite_productivity.rds)
CSPI_W_BGI  <- 0.40
CSPI_W_ASYM <- 0.35
CSPI_W_SI   <- 0.25

## Z-score each component where available
if ("bgi" %in% names(change_all) && any(!is.na(change_all$bgi))) {
  change_all[, bgi_z := (bgi - mean(bgi, na.rm = TRUE)) /
                         sd(bgi, na.rm = TRUE)]
} else {
  change_all[, bgi_z := NA_real_]
}
if ("max_biomass" %in% names(change_all) && any(!is.na(change_all$max_biomass))) {
  change_all[, asym_z := (max_biomass - mean(max_biomass, na.rm = TRUE)) /
                          sd(max_biomass, na.rm = TRUE)]
} else {
  change_all[, asym_z := NA_real_]
}
if ("climate_si" %in% names(change_all) && any(!is.na(change_all$climate_si))) {
  change_all[, si_z := (climate_si - mean(climate_si, na.rm = TRUE)) /
                        sd(climate_si, na.rm = TRUE)]
} else {
  change_all[, si_z := NA_real_]
}

## Weighted composite: available components only, renormalize weights
change_all[, cspi := {
  w_bgi  <- fifelse(!is.na(bgi_z),  CSPI_W_BGI,  0)
  w_asym <- fifelse(!is.na(asym_z), CSPI_W_ASYM, 0)
  w_si   <- fifelse(!is.na(si_z),   CSPI_W_SI,   0)
  w_tot  <- w_bgi + w_asym + w_si
  fifelse(w_tot > 0,
    (fifelse(!is.na(bgi_z),  bgi_z * CSPI_W_BGI,  0) +
     fifelse(!is.na(asym_z), asym_z * CSPI_W_ASYM, 0) +
     fifelse(!is.na(si_z),   si_z * CSPI_W_SI,     0)) / w_tot,
    NA_real_)
}]

n_cspi <- sum(!is.na(change_all$cspi))
cat("    CSPI computed:", format(n_cspi, big.mark = ","),
    "of", format(nrow(change_all), big.mark = ","), "records\n")

## Clean up z-score intermediates
change_all[, c("bgi_z", "asym_z", "si_z") := NULL]

cat("\n")


## ---- Step 4f: Climate PCA extraction from top-4 embedding raster ----------

cat("Step 4f: Extracting climate PCA components from embedding raster...\n")

if (file.exists(CLIM_PCA_TIF)) {
  suppressPackageStartupMessages(library(terra))
  clim_r <- rast(CLIM_PCA_TIF)
  cat("  Raster layers:", nlyr(clim_r), " resolution:", res(clim_r)[1], "m\n")

  ## Unique plot coordinates
  pc <- change_all[!is.na(LAT) & !is.na(LON),
                   .(LAT = LAT[1], LON = LON[1]),
                   by = .(STATECD, UNITCD, COUNTYCD, PLOT)]

  ## Project plots to raster CRS (Lambert Azimuthal Equal Area)
  pts <- vect(pc[, .(LON, LAT)], geom = c("LON", "LAT"), crs = "EPSG:4326")
  pts <- project(pts, crs(clim_r))

  ## Extract all 4 layers
  ex <- extract(clim_r, pts)
  pc[, clim_pca1 := ex[[2]]]
  pc[, clim_pca2 := ex[[3]]]
  pc[, clim_pca3 := ex[[4]]]
  pc[, clim_pca4 := ex[[5]]]

  n_clim <- sum(!is.na(pc$clim_pca1))
  cat("    Climate PCA extracted:", format(n_clim, big.mark = ","),
      "of", format(nrow(pc), big.mark = ","), "plots\n")

  ## Join back to main data
  change_all <- merge(change_all,
                      pc[, .(STATECD, UNITCD, COUNTYCD, PLOT,
                             clim_pca1, clim_pca2, clim_pca3, clim_pca4)],
                      by = c("STATECD", "UNITCD", "COUNTYCD", "PLOT"),
                      all.x = TRUE)
  cat("    Joined to records:",
      format(sum(!is.na(change_all$clim_pca1)), big.mark = ","),
      "of", format(nrow(change_all), big.mark = ","), "\n")
} else {
  cat("  WARNING: Climate PCA raster not found at", CLIM_PCA_TIF, "\n")
  change_all[, c("clim_pca1", "clim_pca2", "clim_pca3", "clim_pca4") := 0]
}

cat("\n")


## ---- Step 5: Standardize column names for scripts 31-34 ------------------

cat("Step 5: Standardizing columns for downstream scripts...\n")

## Map David's CHANGEdata columns to the names expected by scripts 31-34.
## Scripts expect: SPCD, ecodiv_code, DBH1, DBH2, HT1, HT2, CR1, CR2,
##   BAL1, BAL2, BA1, BA2, QMD1, QMD2, CCFL1, CCH1, YEARS, STATUS1,
##   STATUS2, HTCD1, HTCD2, climate_si, bgi, SLOPE, ASPECT, ELEV, etc.

std <- change_all[, .(
  ## Identifiers
  STATECD     = STATECD,
  UNITCD      = UNITCD,
  COUNTYCD    = COUNTYCD,
  PLOT        = PLOT,
  SUBP        = SUBP,
  TREE        = TREE,
  SPCD        = SPCD,
  SPGRPCD     = SPGRPCD,
  plot_key    = plot_key,
  fvs_variant = fvs_variant,
  ecodiv_code = ecodiv_code,

  ## Start of period tree attributes
  DBH1        = startDIA,
  HT1         = startHT,
  CR1         = startCR / 100,       # convert 0-100 to 0-1 proportion
  HTCD1       = startHTCD,
  STATUS1     = startSTATUSCD,
  AGENTCD1    = startAGENTCD,
  TREESTATUS1 = startTreeStatus,

  ## End of period tree attributes
  DBH2        = endDIA,
  HT2         = endHT,
  CR2         = endCR / 100,         # convert 0-100 to 0-1 proportion
  HTCD2       = endHTCD,
  STATUS2     = endSTATUSCD,
  AGENTCD2    = endAGENTCD,
  TREESTATUS2 = endTreeStatus,

  ## Start of period stand summary variables
  BAL1        = startBAL,
  BA1         = startBAPA,
  QMD1        = startQMD,
  TPA1        = startTPA,
  SDI1        = startSDI,
  RD1         = startRD,
  CCF1        = startCCF,
  CCFL1       = startCCFL,
  CCH1        = startCCH,
  HT40_1      = startHT40,
  HTlorey1    = startHTlorey,
  NLIVE1      = startNLive,

  ## End of period stand summary variables
  BAL2        = endBAL,
  BA2         = endBAPA,
  QMD2        = endQMD,
  TPA2        = endTPA,
  SDI2        = endSDI,
  RD2         = endRD,
  CCF2        = endCCF,
  CCFL2       = endCCFL,
  CCH2        = endCCH,

  ## Site and location
  LAT         = LAT,
  LON         = LON,
  ELEV        = ELEV,
  SLOPE       = SLOPE,
  ASPECT      = ASPECT,
  FORTYPCD    = FORTYPCD,
  SICOND      = SICOND,
  SIBASE      = SIBASE,
  SISP        = SISP,
  SICOND_FVS  = SICOND_FVS,
  STDAGE      = STDAGE,

  ## Measurement period
  YEARS       = REMPER,
  YEARS_ADJ   = REMPERadj,
  INVYR1      = startINVYR,
  INVYR2      = endINVYR,
  CONDPROP    = startCONDPROP_UNADJ,

  ## Expansion factors
  TPA_UNADJ1  = startTPA_UNADJ,
  EXPAN1       = startEXPAN,

  ## Data quality flags
  DIA_MISS1    = startDIAmissing,
  HT_MISS1     = startHTmissing,
  CR_MISS1     = startCRmissing,
  CONDID_CHANGE = CONDIDchange,
  DESIGN_CHANGE = DESIGNCDchange,

  ## Site productivity indices (from joins)
  climate_si  = climate_si,
  bgi         = bgi,
  max_biomass = max_biomass,
  byi         = byi,
  SDImax_brms = SDImax_brms,
  cspi        = cspi,

  ## Climate PCA (4 components from CONUS_Climate_Embedding_Top4.tif)
  clim_pca1   = clim_pca1,
  clim_pca2   = clim_pca2,
  clim_pca3   = clim_pca3,
  clim_pca4   = clim_pca4,

  ## EPA Level 3 ecoregion (original code, distinct from ecodiv_code mapping)
  EPA_L3_CODE = EPA_L3_CODE,

  ## FIA Bailey ecoregion subsection (NSBE equation anchor; for comparison)
  ECOSUBCD    = ECOSUBCD
)]

## Compute RD_sdimax where SDImax is available (legacy, from David's startSDI)
std[, rd_sdimax := fifelse(!is.na(SDImax_brms) & SDImax_brms > 0,
                           SDI1 / SDImax_brms, NA_real_)]

## (Climate PCA components clim_pca1-4 populated in Step 4f from the
##  CONUS_Climate_Embedding_Top4.tif raster; no placeholder needed.)

cat("  Standardized dataset:", format(nrow(std), big.mark = ","), "records,",
    ncol(std), "columns\n\n")


## ---- Step 5b: Metric conversion, additive SDI, interactions, EPA L1/L2 ---
##
##   FIA data (David's CHANGEdata) is natively Imperial. All downstream Stan
##   models fit in metric to be comparable with Weiskittel 2016 and Kuehne
##   2020/2022. We convert in place on `std` so scripts 31..34 consume metric
##   directly without any inline conversion logic.
##
##   RD uses the additive Stage SDI formulation requested by Aaron:
##     sdi_tree = (DBH_cm / 25.4)^1.605 * TPH_UNADJ
##   summed to plot then divided by brms.SDImax.mean (already plot-keyed via
##   the join in Step 4d).
##
##   EPA ecoregion codes are hierarchical by construction. NA_L3CODE of the
##   form "X.Y.Z" decomposes into L1 = "X" and L2 = "X.Y" so nested random
##   effects come for free once the L3 column is populated.
## --------------------------------------------------------------------------

cat("Step 5b: Metric conversion + additive SDI + EPA L1/L2 ...\n")

## 5b.1 Conversion factors
IN_TO_CM     <- 2.54
FT_TO_M      <- 0.3048
FTAC_TO_M2HA <- 0.229568     # ft^2/ac -> m^2/ha
TRAC_TO_HA   <- 2.47105       # tr/ac   -> tr/ha

## 5b.2 Tree-level length conversions
for (v in c("DBH1", "DBH2")) {
  if (v %in% names(std)) std[, (v) := get(v) * IN_TO_CM]
}
for (v in c("HT1", "HT2", "CCH1", "CCH2", "HT40_1", "HTlorey1",
            "SICOND", "SICOND_FVS")) {
  if (v %in% names(std)) std[, (v) := get(v) * FT_TO_M]
}

## 5b.3 Stand-level density conversions
for (v in c("BA1", "BA2", "BAL1", "BAL2", "CCFL1", "CCFL2")) {
  if (v %in% names(std)) std[, (v) := get(v) * FTAC_TO_M2HA]
}
if ("QMD1" %in% names(std)) std[, QMD1 := QMD1 * IN_TO_CM]
if ("QMD2" %in% names(std)) std[, QMD2 := QMD2 * IN_TO_CM]

## 5b.4 Trees per hectare (rename to TPH so Imperial / metric never confused)
if ("TPA1" %in% names(std))        std[, TPH1       := TPA1       * TRAC_TO_HA]
if ("TPA2" %in% names(std))        std[, TPH2       := TPA2       * TRAC_TO_HA]
if ("TPA_UNADJ1" %in% names(std))  std[, TPH_UNADJ1 := TPA_UNADJ1 * TRAC_TO_HA]

## 5b.5 Additive SDI per tree (Aaron's formula, cm and tr/ha)
##        sdi_contrib = (DBH_cm / 25.4)^1.605 * TPH_UNADJ
## The 25.4 cm baseline preserves the Reineke 10 inch reference tree even
## though DBH is now in cm. SDImax_brms is on the same scale.
std[, sdi_contrib1 := fifelse(
  !is.na(DBH1) & DBH1 > 0 & !is.na(TPH_UNADJ1) & TPH_UNADJ1 > 0,
  (DBH1 / 25.4)^1.605 * TPH_UNADJ1,
  NA_real_
)]

## Aggregate to plot-level additive SDI (live trees at start only)
plot_sdi <- std[STATUS1 == 1 & !is.na(sdi_contrib1),
                .(sdi_additive1 = sum(sdi_contrib1, na.rm = TRUE)),
                by = plot_key]
std <- merge(std, plot_sdi, by = "plot_key", all.x = TRUE)

## 5b.6 Relative density = additive SDI / brms SDImax
std[, rd_add := fifelse(!is.na(SDImax_brms) & SDImax_brms > 0 &
                        !is.na(sdi_additive1),
                        sdi_additive1 / SDImax_brms, NA_real_)]

cat("  rd_add computed for ",
    format(sum(!is.na(std$rd_add)), big.mark = ","), " of ",
    format(nrow(std), big.mark = ","), " records\n", sep = "")
cat("  rd_add distribution: median=",
    round(median(std$rd_add, na.rm = TRUE), 3),
    "  q95=", round(quantile(std$rd_add, 0.95, na.rm = TRUE), 3), "\n",
    sep = "")

## 5b.7 BAPA*RD and BAL*RD interactions (both now metric)
std[, ba_x_rd  := BA1  * rd_add]
std[, bal_x_rd := BAL1 * rd_add]

## 5b.8 Softwood / hardwood BAL split (for Kuehne 2022 form)
##   If David's CHANGEdata already carries BAL_SW/BAL_HW use them, otherwise
##   use SPGRPCD >= 100 as the softwood flag (FIA convention is SPGRPCD 1..24
##   softwoods; exact numeric ranges vary by RPA release but we use SPCD
##   ranges documented in REF_SPECIES: SPCD 10..299 softwood, 300+ hardwood).
if (!"BAL_SW1" %in% names(std)) {
  std[, is_sw := SPCD < 300]
  ## Plot-level BA by softwood / hardwood
  bawt <- std[STATUS1 == 1,
              .(BA_sw = sum(fifelse(is_sw, BA1 / TPH1 * TPH_UNADJ1, 0),
                             na.rm = TRUE),
                BA_hw = sum(fifelse(!is_sw, BA1 / TPH1 * TPH_UNADJ1, 0),
                             na.rm = TRUE)),
              by = plot_key]
  ## BAL: use David's BAL column directly, then split by species group share.
  ## This is an approximation (true BAL_SW requires ranking all larger trees
  ## by species group); revisit if the Kuehne 2022 form ends up chosen.
  sw_share <- bawt[, .(plot_key,
                       sw_share = BA_sw / pmax(BA_sw + BA_hw, 1e-6))]
  std <- merge(std, sw_share, by = "plot_key", all.x = TRUE)
  std[, BAL_SW1 := BAL1 *      sw_share]
  std[, BAL_HW1 := BAL1 * (1 - sw_share)]
  std[, is_sw := NULL]
  std[, sw_share := NULL]
}

## 5b.9 Parse EPA L1 and L2 from hierarchical NA_L3CODE (e.g. "8.4.1")
std[, EPA_L1_CODE := sub("^([^.]+).*$",           "\\1", EPA_L3_CODE)]
std[, EPA_L2_CODE := sub("^([^.]+\\.[^.]+).*$",   "\\1", EPA_L3_CODE)]
cat("  EPA L1 unique: ", uniqueN(std$EPA_L1_CODE[!is.na(std$EPA_L1_CODE)]),
    "\n")
cat("  EPA L2 unique: ", uniqueN(std$EPA_L2_CODE[!is.na(std$EPA_L2_CODE)]),
    "\n")
cat("  EPA L3 unique: ", uniqueN(std$EPA_L3_CODE[!is.na(std$EPA_L3_CODE)]),
    "\n\n")


## ---- Step 6: Apply model-specific filters and save -----------------------

cat("Step 6: Applying model-specific filters...\n")

## 6a. Base filter: live trees at start with valid core covariates
base <- std[
  STATUS1 == 1 &
  !is.na(DBH1) & DBH1 > 0 &
  !is.na(BA1) & BA1 > 0 &
  !is.na(QMD1) & QMD1 > 0 &
  !is.na(BAL1) &
  !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
  !is.na(ecodiv_code) &
  YEARS >= MIN_REMPER & YEARS <= MAX_REMPER
]

cat("  Base filter (live start, valid covariates):",
    format(nrow(base), big.mark = ","), "\n")

## 6b. DG subset: live at both ends with valid DIA
dg <- base[
  STATUS2 == 1 &
  !is.na(DBH2) & DBH2 > 0
]
cat("  DG records (live-live, valid DIA):", format(nrow(dg), big.mark = ","), "\n")

## 6c. HG subset: live at both ends with field-measured HT
hg <- base[
  STATUS2 == 1 &
  !is.na(HT1) & HT1 > 0 &
  !is.na(HT2) & HT2 > 0 &
  HTCD1 == 1 & HTCD2 == 1
]
cat("  HG records (live-live, HTCD=1):", format(nrow(hg), big.mark = ","), "\n")

## 6d. Mortality subset: live at start, died or survived, no harvest
mort <- base[
  (STATUS2 == 1 | (STATUS2 == 2 & (is.na(AGENTCD2) | AGENTCD2 < 80)))
]
cat("  Mortality records (exclude harvest):",
    format(nrow(mort), big.mark = ","),
    " (", round(100 * mean(mort$STATUS2 == 2, na.rm = TRUE), 2), "% dead)\n")

## Save the unified dataset (base contains everything needed;
## scripts 31-34 apply their own component-specific filters from this)
out_rds <- file.path(OUT_DIR, "conus_remeasurement_pairs.rds")
cat("\n  Saving unified dataset to:", out_rds, "\n")
saveRDS(base, out_rds)

cat("  File size:", round(file.size(out_rds) / 1024^2, 1), "MB\n")


## ---- Step 7: Build HCB cross-sectional dataset from TREEdata ------------

cat("\nStep 7: Building HCB cross-sectional dataset from TREEdata...\n")

read_tree_from_zip <- function(zip_path) {
  st <- substr(basename(zip_path), 1, 2)
  if (EXCLUDE_AK && st == "AK") return(NULL)

  contents <- tryCatch(unzip(zip_path, list = TRUE)$Name,
                       error = function(e) character())
  tree_file <- contents[grepl("TREEdata", contents, ignore.case = TRUE)]
  if (length(tree_file) == 0) return(NULL)

  tmp_dir <- tempdir()
  unzip(zip_path, files = tree_file[1], exdir = tmp_dir, overwrite = TRUE)
  extracted <- file.path(tmp_dir, tree_file[1])

  dt <- tryCatch(
    fread(extracted, showProgress = FALSE, na.strings = c("NA", ""),
          select = c("STATECD", "UNITCD", "COUNTYCD", "PLOT", "SUBP",
                     "TREE", "INVYR", "SPCD", "STATUSCD", "DIA", "HT",
                     "HTCD", "CR", "BAL", "BAPA", "QMD", "SDI", "RD",
                     "CCF", "CCFL", "CCH", "TPA", "LAT", "LON", "ELEV",
                     "SLOPE", "ASPECT", "FORTYPCD", "SICOND", "COND_STATUS_CD",
                     "CONDPROP_UNADJ")),
    error = function(e) NULL
  )

  unlink(extracted)
  if (!is.null(dt)) dt[, STATE_ABBR := st]
  return(dt)
}

cat("  Reading TREEdata files...\n")
tree_list <- map(zip_files, read_tree_from_zip)
tree_list <- compact(tree_list)
tree_all <- rbindlist(tree_list, use.names = TRUE, fill = TRUE)
rm(tree_list); gc(verbose = FALSE)

## Filter for HCB: live trees, valid HT and CR, accessible forest
hcb <- tree_all[
  STATUSCD == 1 &
  COND_STATUS_CD == 1 &
  !is.na(CONDPROP_UNADJ) & CONDPROP_UNADJ >= CONDPROP_MIN &
  !is.na(DIA) & DIA > 0 &
  !is.na(HT) & HT > 0 &
  !is.na(CR) & CR > 0 & CR <= 100 &
  !is.na(BAPA) & BAPA > 0
]

## Compute HCB
hcb[, HCB := HT * (1 - CR / 100)]
hcb[, CR_prop := CR / 100]

## Join ecodivision to HCB data
if (exists("plot_coords") && "ecodiv_code" %in% names(plot_coords)) {
  hcb <- merge(hcb,
               plot_coords[, .(STATECD, UNITCD, COUNTYCD, PLOT,
                               fvs_variant, ecodiv_code)],
               by = c("STATECD", "UNITCD", "COUNTYCD", "PLOT"),
               all.x = TRUE)
}

hcb_rds <- file.path(OUT_DIR, "conus_hcb_cross_section.rds")
cat("  HCB records:", format(nrow(hcb), big.mark = ","), "\n")
cat("  Saving to:", hcb_rds, "\n")
saveRDS(hcb, hcb_rds)

rm(tree_all, hcb); gc(verbose = FALSE)


## ---- Step 8: Summary statistics ------------------------------------------

cat("\nStep 8: Generating summary...\n")

## Reload the base dataset for summary
base <- readRDS(out_rds)

summary_by_state <- base[, .(
  n_records     = .N,
  n_plots       = uniqueN(paste(STATECD, UNITCD, COUNTYCD, PLOT)),
  n_species     = uniqueN(SPCD),
  mean_dbh      = round(mean(DBH1, na.rm = TRUE), 2),
  mean_ba       = round(mean(BA1, na.rm = TRUE), 1),
  mean_remper   = round(mean(YEARS, na.rm = TRUE), 1),
  pct_has_ht    = round(100 * mean(!is.na(HT1)), 1),
  pct_htcd1     = round(100 * mean(HTCD1 == 1, na.rm = TRUE), 1),
  pct_mortality = round(100 * mean(STATUS2 == 2, na.rm = TRUE), 2),
  pct_has_bgi   = round(100 * mean(!is.na(bgi)), 1),
  pct_has_clim  = round(100 * mean(!is.na(climate_si)), 1),
  pct_has_sdimax = round(100 * mean(!is.na(SDImax_brms)), 1)
), by = STATECD]

summary_csv <- file.path(OUT_DIR, "conus_data_summary.csv")
fwrite(summary_by_state, summary_csv)

## Overall summary
cat("\n==============================================================\n")
cat("SUMMARY\n")
cat("==============================================================\n")
cat("Total records (base filter):  ", format(nrow(base), big.mark = ","), "\n")
cat("Total plots:                  ",
    format(uniqueN(base$plot_key), big.mark = ","), "\n")
cat("Total species:                ", uniqueN(base$SPCD), "\n")
cat("EPA L3 ecoregions:            ", uniqueN(base$ecodiv_code), "\n")
cat("States:                       ", uniqueN(base$STATECD), "\n")
cat("Mean remeasurement period:    ", round(mean(base$YEARS), 1), " years\n")
cat("DG-eligible:                  ", format(nrow(base[STATUS2 == 1 & !is.na(DBH2)]),
                                             big.mark = ","), "\n")
cat("HG-eligible (HTCD=1):         ",
    format(nrow(base[STATUS2 == 1 & !is.na(HT1) & !is.na(HT2) &
                     HTCD1 == 1 & HTCD2 == 1]), big.mark = ","), "\n")
cat("Mortality-eligible:           ",
    format(nrow(base[(STATUS2 == 1) |
                     (STATUS2 == 2 & (is.na(AGENTCD2) | AGENTCD2 < 80))]),
           big.mark = ","), "\n")

cat("\nSite productivity coverage:\n")
cat("  ClimateSI:   ", round(100 * mean(!is.na(base$climate_si)), 1), "%\n")
cat("  BGI:         ", round(100 * mean(!is.na(base$bgi)), 1), "%\n")
cat("  max_biomass: ", round(100 * mean(!is.na(base$max_biomass)), 1), "%\n")
cat("  BYI:         ", round(100 * mean(!is.na(base$byi)), 1), "%\n")
cat("  SDImax_brms: ", round(100 * mean(!is.na(base$SDImax_brms)), 1), "%\n")
cat("  CSPI:        ", round(100 * mean(!is.na(base$cspi)), 1), "%\n")
cat("  Climate PCA: ", round(100 * mean(!is.na(base$clim_pca1)), 1), "%\n")

cat("\nDerived density / competition (metric, post Step 5b):\n")
cat("  rd_add coverage:   ",
    round(100 * mean(!is.na(base$rd_add)), 1), "%   median=",
    round(median(base$rd_add, na.rm = TRUE), 3), "\n", sep = "")
cat("  ba_x_rd  mean:     ",
    round(mean(base$ba_x_rd,  na.rm = TRUE), 2), " m2/ha * RD\n")
cat("  bal_x_rd mean:     ",
    round(mean(base$bal_x_rd, na.rm = TRUE), 2), " m2/ha * RD\n")
cat("  BAL_SW1  mean:     ",
    round(mean(base$BAL_SW1,  na.rm = TRUE), 2), " m2/ha\n")
cat("  BAL_HW1  mean:     ",
    round(mean(base$BAL_HW1,  na.rm = TRUE), 2), " m2/ha\n")

cat("\nEPA ecoregion nesting:\n")
cat("  EPA L1 unique: ", uniqueN(base$EPA_L1_CODE), "\n")
cat("  EPA L2 unique: ", uniqueN(base$EPA_L2_CODE), "\n")
cat("  EPA L3 unique: ", uniqueN(base$EPA_L3_CODE), "\n")

cat("\nOutputs:\n")
cat("  ", out_rds, "\n")
cat("  ", hcb_rds, "\n")
cat("  ", summary_csv, "\n")
cat("\n==============================================================\n")
cat("Done. Time:", format(Sys.time()), "\n")
cat("==============================================================\n")
