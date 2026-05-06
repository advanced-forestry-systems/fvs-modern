## =============================================================================
## 30d_fix_units.R
##
## FVS-CONUS: surgical unit-conversion fix for the existing "metric" RDS.
##
## Bug: 30b_build_metric_patch.R used a DBH-based guard
##
##   if (dbh_max > 100) metric_already <- TRUE
##
## and short-circuited every HT / BA / BAL / CCFL / QMD / SICOND conversion
## whenever the upstream RDS already had DBH in cm. Diagnosis (2026-04-22)
## against calibration/data/conus_remeasurement_pairs_metric.rds:
##
##   DBH1  max = 191       (cm, correctly converted)
##   HT1   max = 309       (ft, NOT converted; 309 m is physically impossible)
##   BA1   max = 1289      (ft2/ac, NOT converted)
##   BAL1  max = 1288      (ft2/ac, NOT converted)
##   CCFL1 max = 1119      (ft2/ac, NOT converted)
##   QMD1  max = 66.4      (in, NOT converted)
##   TPH1  max = 37307     (trees/ha, correctly converted from TPA * 2.47)
##
## This script reads the hybrid RDS, converts only the variables that are
## still imperial (detected by per-variable sanity thresholds), recomputes
## the downstream interactions that consumed them (ba_x_rd, bal_x_rd,
## BAL_SW1, BAL_HW1), and writes a fixed RDS. rd_add itself is dimensionless
## and unchanged (computed from DBH_cm and TPH_unadj_ha only).
##
## Idempotent: runs its own per-variable sanity checks; if nothing is left
## to convert, exits cleanly.
##
## Usage:
##   Rscript calibration/R/30d_fix_units.R \
##       --in  calibration/data/conus_remeasurement_pairs_metric.rds \
##       --out calibration/data/conus_remeasurement_pairs_metric_v2.rds
##
## Author: A. Weiskittel, 2026-04-22
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--in",  type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds",
              dest    = "in_path"),
  make_option("--out", type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_v2.rds",
              dest    = "out_path")
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS 30d: surgical unit fix\n")
cat("Input :", opts$in_path, "\n")
cat("Output:", opts$out_path, "\n")
cat("==============================================================\n")

stopifnot(file.exists(opts$in_path))

d <- readRDS(opts$in_path)
setDT(d)
cat(sprintf("Loaded %s rows, %d cols\n",
            format(nrow(d), big.mark = ","), ncol(d)))

## ---- Conversion factors ----------------------------------------------------
IN_TO_CM     <- 2.54
FT_TO_M      <- 0.3048
FTAC_TO_M2HA <- 0.229568           # ft^2 / ac  ->  m^2 / ha

## ---- Per-variable state detection -----------------------------------------
## Thresholds rely on physical / biological limits. If the maximum observed
## value exceeds what is physically plausible in metric, the variable is
## still imperial.
##
##   HT    max plausible   : Hyperion at 115.6 m, use 150 as threshold
##   CCH   crown closure ht: same scale as HT
##   QMD   max plausible cm: ~200 (giant sequoia stand), use 150 threshold
##   BA    max plausible   : ~300 m^2/ha in redwood, use 500 threshold
##   BAL   same as BA
##   CCFL  max plausible   : ~800 as CCF can exceed 100, use 500 threshold
##   SICOND max plausible m: ~50 m site index, use 60 threshold
##
## If max exceeds threshold, apply conversion. NA-only columns are skipped.
convert_if_imperial <- function(dt, cols, factor, threshold, label) {
  for (v in cols) {
    if (!v %in% names(dt)) next
    mx <- suppressWarnings(max(dt[[v]], na.rm = TRUE))
    if (!is.finite(mx)) {
      cat(sprintf("  %-12s all NA, skipping\n", v))
      next
    }
    if (mx > threshold) {
      dt[, (v) := get(v) * factor]
      new_mx <- suppressWarnings(max(dt[[v]], na.rm = TRUE))
      cat(sprintf("  %-12s was %8.2f (%s), now %8.2f\n",
                  v, mx, label, new_mx))
    } else {
      cat(sprintf("  %-12s max %8.2f already metric, leaving\n", v, mx))
    }
  }
}

cat("\nStep 1. Height-like variables (ft -> m, threshold 150)\n")
convert_if_imperial(
  d,
  c("HT1", "HT2", "CCH1", "CCH2", "HT40_1", "HTlorey1"),
  FT_TO_M, threshold = 150, label = "ft"
)

cat("\nStep 2. Site-index height (ft -> m, threshold 60)\n")
convert_if_imperial(
  d,
  c("SICOND", "SICOND_FVS"),
  FT_TO_M, threshold = 60, label = "ft"
)

cat("\nStep 3. QMD (in -> cm, threshold 100)\n")
## QMD max in cm could be ~150 for giant sequoia, but p99 stays tight. Use
## a per-variable check: if max > 100 AND p99 > 40, likely cm (already
## metric). Our data shows max = 66, p99 = 16 for QMD1, which is clearly
## inches. Safer threshold: anything with max < 100 is likely still imperial.
for (v in c("QMD1", "QMD2")) {
  if (!v %in% names(d)) next
  mx <- suppressWarnings(max(d[[v]], na.rm = TRUE))
  p99 <- suppressWarnings(quantile(d[[v]], 0.99, na.rm = TRUE))
  if (!is.finite(mx)) next
  if (mx < 100 && p99 < 40) {
    d[, (v) := get(v) * IN_TO_CM]
    cat(sprintf("  %-12s was max=%6.2f p99=%5.2f (in), now max=%6.2f\n",
                v, mx, p99,
                suppressWarnings(max(d[[v]], na.rm = TRUE))))
  } else {
    cat(sprintf("  %-12s max %.2f already cm, leaving\n", v, mx))
  }
}

cat("\nStep 4. Basal-area-like variables (ft2/ac -> m2/ha, threshold 500)\n")
convert_if_imperial(
  d,
  c("BA1", "BA2", "BAL1", "BAL2", "CCFL1", "CCFL2"),
  FTAC_TO_M2HA, threshold = 500, label = "ft2/ac"
)

## ---- Recompute / fix dependents that consumed imperial BA / BAL -----------
##
## rd_add is DBH-and-TPH based, so it is correct and untouched.
## ba_x_rd, bal_x_rd, BAL_SW1, BAL_HW1 all used BA1 / BAL1 in their original
## (imperial) units. We rebuild them from the now-metric values.

cat("\nStep 5. Rebuild dependent interactions from metric BA / BAL\n")
if (all(c("BA1", "rd_add") %in% names(d))) {
  d[, ba_x_rd := BA1 * rd_add]
  cat(sprintf("  ba_x_rd   rebuilt: median=%.3f  p95=%.3f\n",
              median(d$ba_x_rd, na.rm = TRUE),
              quantile(d$ba_x_rd, 0.95, na.rm = TRUE)))
}
if (all(c("BAL1", "rd_add") %in% names(d))) {
  d[, bal_x_rd := BAL1 * rd_add]
  cat(sprintf("  bal_x_rd  rebuilt: median=%.3f  p95=%.3f\n",
              median(d$bal_x_rd, na.rm = TRUE),
              quantile(d$bal_x_rd, 0.95, na.rm = TRUE)))
}

## Rebuild BAL_SW1 / BAL_HW1 by re-applying the sw_share from step 7 of 30b.
## Softwood share is computed from BA contributions per species. The share
## itself is dimensionless, so whether BA is in ft2/ac or m2/ha does not
## change the ratio. We recompute anyway to be safe, then apply to metric
## BAL1.
if (all(c("BAL1", "BA1", "TPH1", "TPH_UNADJ1",
          "plot_key", "SPCD", "STATUS1") %in% names(d))) {
  cat("  Rebuilding BAL_SW1 / BAL_HW1 from softwood share ...\n")
  d[, is_sw := SPCD < 300]
  bawt <- d[STATUS1 == 1 & !is.na(TPH1) & TPH1 > 0,
            .(BA_sw = sum(fifelse(is_sw,  BA1 / TPH1 * TPH_UNADJ1, 0),
                           na.rm = TRUE),
              BA_hw = sum(fifelse(!is_sw, BA1 / TPH1 * TPH_UNADJ1, 0),
                           na.rm = TRUE)),
            by = plot_key]
  sw_share <- bawt[, .(plot_key,
                       sw_share_new = BA_sw / pmax(BA_sw + BA_hw, 1e-6))]
  d <- merge(d, sw_share, by = "plot_key", all.x = TRUE)
  d[, BAL_SW1 := BAL1 *      sw_share_new]
  d[, BAL_HW1 := BAL1 * (1 - sw_share_new)]
  d[, is_sw := NULL]
  d[, sw_share_new := NULL]
  cat(sprintf("  BAL_SW1 / BAL_HW1: mean=%.2f / %.2f m2/ha\n",
              mean(d$BAL_SW1, na.rm = TRUE),
              mean(d$BAL_HW1, na.rm = TRUE)))
}

## ---- Attach a units attribute so downstream code can assert it ------------
attr(d, "units") <- list(
  DBH   = "cm",
  HT    = "m",
  BA    = "m2/ha",
  BAL   = "m2/ha",
  CCFL  = "m2/ha",
  QMD   = "cm",
  TPH   = "trees/ha",
  SI    = "m",
  fixed_on = format(Sys.time(), "%Y-%m-%d %H:%M"),
  fixed_by = "30d_fix_units.R"
)

## ---- Sanity spot-check -----------------------------------------------------
cat("\nSanity spot-check (should all be metric now):\n")
pct <- function(x, p = c(0.01, 0.5, 0.95, 0.99, 1.0)) {
  q <- quantile(x, p, na.rm = TRUE)
  paste0(names(q), "=", round(q, 2), collapse = " ")
}
for (v in c("DBH1", "HT1", "BA1", "BAL1", "CCFL1", "QMD1",
            "SICOND", "TPH1", "rd_add", "ba_x_rd", "bal_x_rd",
            "BAL_SW1", "BAL_HW1")) {
  if (v %in% names(d)) cat(sprintf("  %-12s %s\n", v, pct(d[[v]])))
}

cat(sprintf("\nWriting fixed RDS to %s ...\n", opts$out_path))
saveRDS(d, opts$out_path)
cat(sprintf("Done. %s rows, %d cols\n",
            format(nrow(d), big.mark = ","), ncol(d)))
