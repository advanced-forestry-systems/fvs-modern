## =============================================================================
## 30b_build_metric_patch.R
##
## FVS-CONUS: in-place metric + additive-SDI patch for an already-built
## conus_remeasurement_pairs.rds. This runs only Step 5b of
## 30_build_conus_dataset.R and writes a new RDS alongside the original, so
## the multi-hour re-read of 48 state zips is avoided.
##
## Handles the David-Marshall CHANGEdata column convention where TPA lives in
## `EXPAN1` / `EXPAN2` rather than `TPA_UNADJ1` / `TPA_UNADJ2`.
##
## Idempotent: if rd_add already exists, exits without changes.
##
## Usage:
##   Rscript scripts/30b_build_metric_patch.R \
##       --in  calibration/data/conus_remeasurement_pairs.rds \
##       --out calibration/data/conus_remeasurement_pairs_metric.rds
##
## Author: A. Weiskittel, 2026-04-15
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--in",  type = "character",
              default = "calibration/data/conus_remeasurement_pairs.rds",
              dest    = "in_path"),
  make_option("--out", type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds",
              dest    = "out_path")
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS metric + additive SDI patch\n")
cat("Input :", opts$in_path, "\n")
cat("Output:", opts$out_path, "\n")
cat("==============================================================\n")

stopifnot(file.exists(opts$in_path))

std <- readRDS(opts$in_path)
setDT(std)
cat(sprintf("Loaded %s rows, %d cols\n",
            format(nrow(std), big.mark = ","), ncol(std)))

## Idempotency: already patched?
if ("rd_add" %in% names(std) && "ba_x_rd" %in% names(std) &&
    "EPA_L1_CODE" %in% names(std)) {
  cat("rd_add / ba_x_rd / EPA_L1_CODE already present. No patch needed.\n")
  if (!file.exists(opts$out_path)) saveRDS(std, opts$out_path)
  quit(status = 0)
}

## ---- 0. Rename EXPAN* -> TPA_UNADJ* where needed ---------------------------
if (!"TPA_UNADJ1" %in% names(std) && "EXPAN1" %in% names(std)) {
  cat("Renaming EXPAN1 -> TPA_UNADJ1\n")
  setnames(std, "EXPAN1", "TPA_UNADJ1")
}
if (!"TPA_UNADJ2" %in% names(std) && "EXPAN2" %in% names(std)) {
  setnames(std, "EXPAN2", "TPA_UNADJ2")
}

## ---- 1. Conversion factors -------------------------------------------------
IN_TO_CM     <- 2.54
FT_TO_M      <- 0.3048
FTAC_TO_M2HA <- 0.229568        # ft^2/ac -> m^2/ha
TRAC_TO_HA   <- 2.47105         # tr/ac   -> tr/ha

## Guard against double-conversion. If DBH1 max > 80, probably cm already.
dbh_max <- suppressWarnings(max(std$DBH1, na.rm = TRUE))
if (is.finite(dbh_max) && dbh_max > 100) {
  cat(sprintf("DBH1 max = %.1f (looks metric already). Skipping unit conversions.\n",
              dbh_max))
  metric_already <- TRUE
} else {
  cat(sprintf("DBH1 max = %.1f (looks Imperial). Applying conversions.\n",
              dbh_max))
  metric_already <- FALSE
}

if (!metric_already) {
  ## 2. Length conversions
  for (v in c("DBH1", "DBH2")) {
    if (v %in% names(std)) std[, (v) := get(v) * IN_TO_CM]
  }
  for (v in c("HT1", "HT2", "CCH1", "CCH2", "HT40_1", "HTlorey1",
              "SICOND", "SICOND_FVS")) {
    if (v %in% names(std)) std[, (v) := get(v) * FT_TO_M]
  }
  ## 3. Stand-level density conversions
  for (v in c("BA1", "BA2", "BAL1", "BAL2", "CCFL1", "CCFL2")) {
    if (v %in% names(std)) std[, (v) := get(v) * FTAC_TO_M2HA]
  }
  for (v in c("QMD1", "QMD2")) {
    if (v %in% names(std)) std[, (v) := get(v) * IN_TO_CM]
  }
}

## 4. TPH columns (create if absent)
if ("TPA1" %in% names(std) && !"TPH1" %in% names(std))
  std[, TPH1 := TPA1 * TRAC_TO_HA]
if ("TPA2" %in% names(std) && !"TPH2" %in% names(std))
  std[, TPH2 := TPA2 * TRAC_TO_HA]
if ("TPA_UNADJ1" %in% names(std) && !"TPH_UNADJ1" %in% names(std))
  std[, TPH_UNADJ1 := TPA_UNADJ1 * TRAC_TO_HA]

## ---- 5. Additive SDI -------------------------------------------------------
stopifnot("TPH_UNADJ1" %in% names(std),
          "DBH1"       %in% names(std),
          "plot_key"   %in% names(std),
          "STATUS1"    %in% names(std))

std[, sdi_contrib1 := fifelse(
  !is.na(DBH1) & DBH1 > 0 & !is.na(TPH_UNADJ1) & TPH_UNADJ1 > 0,
  (DBH1 / 25.4)^1.605 * TPH_UNADJ1,
  NA_real_
)]

plot_sdi <- std[STATUS1 == 1 & !is.na(sdi_contrib1),
                .(sdi_additive1 = sum(sdi_contrib1, na.rm = TRUE)),
                by = plot_key]
std <- merge(std, plot_sdi, by = "plot_key", all.x = TRUE)

if (!"SDImax_brms" %in% names(std)) {
  stop("SDImax_brms column missing. Step 4d join in the full build was skipped.")
}
std[, rd_add := fifelse(!is.na(SDImax_brms) & SDImax_brms > 0 &
                        !is.na(sdi_additive1),
                        sdi_additive1 / SDImax_brms, NA_real_)]

cat(sprintf("  rd_add: non-NA = %s of %s   median=%.3f  q95=%.3f\n",
            format(sum(!is.na(std$rd_add)), big.mark = ","),
            format(nrow(std),               big.mark = ","),
            median(std$rd_add, na.rm = TRUE),
            quantile(std$rd_add, 0.95, na.rm = TRUE)))

## ---- 6. Interactions -------------------------------------------------------
std[, ba_x_rd  := BA1  * rd_add]
std[, bal_x_rd := BAL1 * rd_add]

## ---- 7. BAL_SW / BAL_HW (SPCD<300 approx softwood) ------------------------
if (!"BAL_SW1" %in% names(std)) {
  std[, is_sw := SPCD < 300]
  bawt <- std[STATUS1 == 1 & !is.na(TPH1) & TPH1 > 0,
              .(BA_sw = sum(fifelse(is_sw, BA1 / TPH1 * TPH_UNADJ1, 0),
                             na.rm = TRUE),
                BA_hw = sum(fifelse(!is_sw, BA1 / TPH1 * TPH_UNADJ1, 0),
                             na.rm = TRUE)),
              by = plot_key]
  sw_share <- bawt[, .(plot_key,
                       sw_share = BA_sw / pmax(BA_sw + BA_hw, 1e-6))]
  std <- merge(std, sw_share, by = "plot_key", all.x = TRUE)
  std[, BAL_SW1 := BAL1 *      sw_share]
  std[, BAL_HW1 := BAL1 * (1 - sw_share)]
  std[, is_sw := NULL]
  std[, sw_share := NULL]
  cat(sprintf("  BAL_SW1 / BAL_HW1: mean=%.2f / %.2f m^2/ha\n",
              mean(std$BAL_SW1, na.rm = TRUE),
              mean(std$BAL_HW1, na.rm = TRUE)))
}

## ---- 8. EPA L1, L2 ---------------------------------------------------------
if (!"EPA_L3_CODE" %in% names(std) && "ecodiv_code" %in% names(std)) {
  setnames(std, "ecodiv_code", "EPA_L3_CODE")
}
stopifnot("EPA_L3_CODE" %in% names(std))

## EPA_L3_CODE contains a mix of proper hierarchical "X.Y.Z" codes and
## dotless numeric L3-only codes or "UNKNOWN". Parse only the hierarchical
## rows; leave the rest NA so prepare_dg_data() drops them downstream.
has_hier <- grepl("^[^.]+\\.[^.]+\\.[^.]+$", std$EPA_L3_CODE)
std[, EPA_L1_CODE := NA_character_]
std[, EPA_L2_CODE := NA_character_]
std[has_hier, EPA_L1_CODE := sub("^([^.]+)\\..*$",         "\\1", EPA_L3_CODE)]
std[has_hier, EPA_L2_CODE := sub("^([^.]+\\.[^.]+)\\..*$", "\\1", EPA_L3_CODE)]
## Also null out the fake L3 codes so all three indices are consistent
std[!has_hier, EPA_L3_CODE := NA_character_]

cat(sprintf("  Hierarchical EPA codes: %s (%.1f%% of rows)\n",
            format(sum(has_hier), big.mark = ","),
            100 * mean(has_hier)))

cat(sprintf("  EPA L1 unique: %d,  L2 unique: %d,  L3 unique: %d\n",
            uniqueN(std$EPA_L1_CODE[!is.na(std$EPA_L1_CODE)]),
            uniqueN(std$EPA_L2_CODE[!is.na(std$EPA_L2_CODE)]),
            uniqueN(std$EPA_L3_CODE[!is.na(std$EPA_L3_CODE)])))

## ---- 9. Save --------------------------------------------------------------
cat(sprintf("Writing patched RDS to %s ...\n", opts$out_path))
saveRDS(std, opts$out_path)
cat(sprintf("Done. %s rows, %d cols\n",
            format(nrow(std), big.mark = ","), ncol(std)))
