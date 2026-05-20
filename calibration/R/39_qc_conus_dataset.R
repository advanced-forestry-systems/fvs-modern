## =============================================================================
## 39_qc_conus_dataset.R
##
## FVS-CONUS: Data quality checks on conus_remeasurement_pairs.rds
##
## Reports:
##   1. Coverage by state, EPA L3 ecoregion, species
##   2. Outlier detection on DBH, HT, CR, growth increments
##   3. Distribution of remeasurement periods
##   4. Site productivity joint coverage and correlations
##   5. Records per species x ecoregion cell (identifies sparse combinations)
##   6. Spatial coverage sanity check
##   7. Data integrity flags (impossible values, negative growth, etc.)
##
## Usage:
##   Rscript calibration/R/39_qc_conus_dataset.R
##   Rscript calibration/R/39_qc_conus_dataset.R --data path/to.rds --out_dir qc/
##
## Outputs (to OUT_DIR):
##   conus_qc_summary.csv       numerical summary by state
##   conus_qc_species_cells.csv records per species x ecoregion
##   conus_qc_outliers.csv      flagged records with context
##   conus_qc_report.txt        human-readable log
##
## Author: A. Weiskittel, 2026-04-15
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

## ---- Configuration ---------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
DATA    <- "calibration/data/conus_remeasurement_pairs.rds"
OUT_DIR <- "calibration/output/conus/qc"

for (i in seq_along(args)) {
  if (args[i] == "--data"    && i < length(args)) DATA    <- args[i + 1]
  if (args[i] == "--out_dir" && i < length(args)) OUT_DIR <- args[i + 1]
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
report_file <- file.path(OUT_DIR, "conus_qc_report.txt")
sink(report_file, split = TRUE)

cat("==============================================================\n")
cat("FVS-CONUS Data QC Report\n")
cat("Generated:", format(Sys.time()), "\n")
cat("Input:    ", DATA, "\n")
cat("==============================================================\n\n")

## ---- Load data -------------------------------------------------------------
stopifnot(file.exists(DATA))
dat <- readRDS(DATA)
setDT(dat)

cat("Records:  ", format(nrow(dat), big.mark = ","), "\n")
cat("Columns:  ", ncol(dat), "\n")
cat("Plots:    ", format(uniqueN(dat$plot_key), big.mark = ","), "\n\n")

## ---- 1. Summary by state ---------------------------------------------------
cat("---- 1. Summary by state ----\n")
by_state <- dat[, .(
  n_records     = .N,
  n_plots       = uniqueN(plot_key),
  n_species     = uniqueN(SPCD),
  n_ecoregions  = uniqueN(ecodiv_code),
  mean_dbh      = round(mean(DBH1, na.rm = TRUE), 2),
  mean_ht       = round(mean(HT1, na.rm = TRUE), 1),
  mean_remper   = round(mean(YEARS, na.rm = TRUE), 2),
  pct_mortality = round(100 * mean(STATUS2 == 2, na.rm = TRUE), 2),
  pct_has_bgi   = round(100 * mean(!is.na(bgi)), 1),
  pct_has_cspi  = round(100 * mean(!is.na(cspi)), 1),
  pct_has_clim  = round(100 * mean(!is.na(clim_pca1)), 1)
), by = STATECD][order(STATECD)]

print(by_state)
fwrite(by_state, file.path(OUT_DIR, "conus_qc_summary.csv"))
cat("\n")

## ---- 2. Outlier detection --------------------------------------------------
cat("---- 2. Outlier detection ----\n")

flags <- dat[, .(
  plot_key, SPCD, STATECD, SUBP, TREE,
  DBH1, DBH2, HT1, HT2, CR1, YEARS, STATUS1, STATUS2
)]

## Flag impossible values and implausible growth
flags[, flag_impossible_dbh := (DBH1 <= 0 | (!is.na(DBH2) & DBH2 <= 0))]
flags[, flag_impossible_ht  := (!is.na(HT1) & HT1 <= 0) |
                               (!is.na(HT2) & HT2 <= 0)]
flags[, flag_cr_oob         := !is.na(CR1) & (CR1 < 0 | CR1 > 1)]

## Negative DBH growth in surviving trees (acceptable small amount from
## recalibration, but large shrinkage is suspect)
flags[, dg_annual := (DBH2 - DBH1) / YEARS]
flags[, flag_dbh_shrink := STATUS2 == 1 & !is.na(DBH2) &
                          (DBH2 - DBH1) < -1.0]
## Extreme positive growth (> 3 cm/yr sustained)
flags[, flag_dbh_extreme := STATUS2 == 1 & !is.na(DBH2) &
                           dg_annual > 3.0]

## HT shrinkage > 2 m (surviving trees, field measured both periods only)
flags[, flag_ht_shrink := STATUS2 == 1 & !is.na(HT2) &
                         (HT2 - HT1) < -2.0]

flag_counts <- data.table(
  flag = c("impossible DBH", "impossible HT", "CR out of bounds",
           "DBH shrinkage > 1 cm", "DBH growth > 3 cm/yr",
           "HT shrinkage > 2 m"),
  n    = c(sum(flags$flag_impossible_dbh, na.rm = TRUE),
           sum(flags$flag_impossible_ht,  na.rm = TRUE),
           sum(flags$flag_cr_oob,         na.rm = TRUE),
           sum(flags$flag_dbh_shrink,     na.rm = TRUE),
           sum(flags$flag_dbh_extreme,    na.rm = TRUE),
           sum(flags$flag_ht_shrink,      na.rm = TRUE))
)
flag_counts[, pct := round(100 * n / nrow(flags), 4)]
print(flag_counts)

## Save first 500 flagged outliers for inspection
flagged <- flags[flag_impossible_dbh | flag_impossible_ht | flag_cr_oob |
                 flag_dbh_shrink | flag_dbh_extreme | flag_ht_shrink]
fwrite(head(flagged, 500), file.path(OUT_DIR, "conus_qc_outliers.csv"))
cat("\nSaved first 500 flagged records to conus_qc_outliers.csv\n\n")

## ---- 3. Remeasurement period distribution ----------------------------------
cat("---- 3. Remeasurement period (YEARS) ----\n")
cat("  min:    ", min(dat$YEARS, na.rm = TRUE), "\n")
cat("  q25:    ", quantile(dat$YEARS, 0.25, na.rm = TRUE), "\n")
cat("  median: ", median(dat$YEARS, na.rm = TRUE), "\n")
cat("  mean:   ", round(mean(dat$YEARS, na.rm = TRUE), 2), "\n")
cat("  q75:    ", quantile(dat$YEARS, 0.75, na.rm = TRUE), "\n")
cat("  max:    ", max(dat$YEARS, na.rm = TRUE), "\n\n")

## ---- 4. Site productivity coverage and correlations ------------------------
cat("---- 4. Site productivity coverage ----\n")
for (v in c("climate_si", "bgi", "max_biomass", "byi", "SDImax_brms", "cspi")) {
  n_na <- sum(is.na(dat[[v]]))
  cat(sprintf("  %-14s: %5.1f%% missing (n=%s)\n",
              v, 100 * n_na / nrow(dat), format(n_na, big.mark = ",")))
}
cat("\nPairwise Pearson correlations (site productivity):\n")
sp_cols <- c("climate_si", "bgi", "max_biomass", "byi", "cspi")
sp_dat <- dat[, ..sp_cols]
cor_mat <- round(cor(sp_dat, use = "pairwise.complete.obs"), 3)
print(cor_mat)
cat("\n")

## ---- 5. Species x ecoregion cell occupancy ---------------------------------
cat("---- 5. Species x ecoregion cell occupancy ----\n")
cells <- dat[!is.na(SPCD) & !is.na(ecodiv_code),
             .(n = .N), by = .(SPCD, ecodiv_code)][order(-n)]
cat("  Total cells:           ",
    format(nrow(cells), big.mark = ","), "\n")
cat("  Cells with <  50 obs:  ",
    format(sum(cells$n < 50), big.mark = ","),
    "(", round(100 * mean(cells$n < 50), 1), "%)\n")
cat("  Cells with < 500 obs:  ",
    format(sum(cells$n < 500), big.mark = ","),
    "(", round(100 * mean(cells$n < 500), 1), "%)\n")
cat("  Median obs per cell:   ", median(cells$n), "\n")
cat("  Max obs per cell:      ", max(cells$n), "\n\n")
fwrite(cells, file.path(OUT_DIR, "conus_qc_species_cells.csv"))

## By species total and species present in >= 3 ecoregions
by_sp <- dat[, .(n_records = .N,
                 n_eco = uniqueN(ecodiv_code),
                 n_plots = uniqueN(plot_key)),
             by = SPCD][order(-n_records)]
cat("Top 20 species by record count:\n")
print(head(by_sp, 20))
cat("\n")
cat("Species with < 5,000 records (rare species pooling candidates):",
    sum(by_sp$n_records < 5000), "\n\n")

## ---- 6. Spatial extent sanity ----------------------------------------------
cat("---- 6. Spatial extent ----\n")
cat("  LAT range:", round(range(dat$LAT, na.rm = TRUE), 3), "\n")
cat("  LON range:", round(range(dat$LON, na.rm = TRUE), 3), "\n")
cat("  Records with LAT/LON NA:", sum(is.na(dat$LAT) | is.na(dat$LON)), "\n\n")

## ---- 7. ECOSUBCD (FIA Bailey) vs EPA L3 comparison -------------------------
cat("---- 7. ECOSUBCD (FIA Bailey) vs EPA L3 comparison ----\n")
if ("ECOSUBCD" %in% names(dat) && any(!is.na(dat$ECOSUBCD))) {
  cat("  Unique ECOSUBCD values: ", uniqueN(dat$ECOSUBCD), "\n")
  cat("  Unique EPA L3 values:   ", uniqueN(dat$ecodiv_code), "\n")

  ## Records per cell in each scheme
  eco_bailey <- dat[, .N, by = ECOSUBCD]
  eco_epa    <- dat[, .N, by = ecodiv_code]
  cat(sprintf("  Median records per cell -- Bailey: %d, EPA L3: %d\n",
              median(eco_bailey$N), median(eco_epa$N)))

  ## Cross-tab: how many EPA L3 codes map to each Bailey subsection?
  xtab <- dat[!is.na(ECOSUBCD) & !is.na(ecodiv_code),
              .(n = .N), by = .(ECOSUBCD, ecodiv_code)]
  map_depth <- xtab[, .(n_epa_codes = uniqueN(ecodiv_code)),
                    by = ECOSUBCD][order(-n_epa_codes)]
  cat("  Each ECOSUBCD splits across N EPA L3 codes (top 10):\n")
  print(head(map_depth, 10))

  inv_depth <- xtab[, .(n_bailey = uniqueN(ECOSUBCD)),
                    by = ecodiv_code][order(-n_bailey)]
  cat("  Each EPA L3 splits across N Bailey subsections (top 10):\n")
  print(head(inv_depth, 10))

  fwrite(xtab, file.path(OUT_DIR, "conus_qc_eco_crosswalk.csv"))
} else {
  cat("  ECOSUBCD not present or all NA. Skipping comparison.\n")
}
cat("\n")


## ---- 8. Summary / overall status -------------------------------------------
cat("==============================================================\n")
cat("QC Summary\n")
cat("==============================================================\n")
cat("Unique EPA L3 ecoregions:", uniqueN(dat$ecodiv_code),
    " (expect 84 after fix)\n")
if ("ECOSUBCD" %in% names(dat)) {
  cat("Unique FIA ECOSUBCDs:    ", uniqueN(dat$ECOSUBCD), "\n")
}
cat("Unique species (SPCD):   ", uniqueN(dat$SPCD), "\n")
cat("Total flagged outliers:  ",
    format(nrow(flagged), big.mark = ","),
    " (", round(100 * nrow(flagged) / nrow(dat), 3), "%)\n")
cat("==============================================================\n")

sink()
cat("\nQC report saved to:", report_file, "\n")
cat("Outputs in:         ", OUT_DIR, "\n")
