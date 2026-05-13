## =============================================================================
## 40_build_cond_from_csv.R
##
## Build calibration/data/fia_cond.rds from the CONUS-scale ENTIRE_COND.csv
## file that feeds the ECOSUBCD join in 30_build_conus_dataset.R.
##
## Keeps the composite plot-measurement key (STATECD, UNITCD, COUNTYCD, PLOT,
## INVYR, CONDID) plus all Phase 2 modifier columns. Reduces a 775 MB CSV
## (2.3M rows) to a compact RDS ready for join against the pairs file.
##
## Default behavior: keep the DOMINANT forested condition per plot-
## measurement (largest CONDPROP_UNADJ among COND_STATUS_CD==1 rows).
## This matches the convention used elsewhere in the pipeline.
##
## Usage:
##   Rscript calibration/R/40_build_cond_from_csv.R \
##       --cond_csv ~/FIA/ENTIRE_COND.csv \
##       --out      calibration/data/fia_cond.rds
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--cond_csv", type = "character",
              default = "~/FIA/ENTIRE_COND.csv",
              dest    = "cond_csv"),
  make_option("--out",      type = "character",
              default = "calibration/data/fia_cond.rds",
              dest    = "out_path"),
  make_option("--keep_nonforest", action = "store_true", default = FALSE,
              help = "Keep COND_STATUS_CD != 1 rows (default drops them)")
)) |> parse_args()

cat("==============================================================\n")
cat("Build fia_cond.rds from ENTIRE_COND.csv\n")
cat("  cond_csv:", opts$cond_csv, "\n")
cat("  out     :", opts$out_path, "\n")
cat("==============================================================\n")

stopifnot(file.exists(opts$cond_csv))

keep_cols <- c(
  ## Composite key, what the pairs file can join on
  "STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR", "CONDID",
  ## Also keep CN / PLT_CN for any downstream work that has them
  "CN", "PLT_CN",
  ## Status and dominant-condition selector
  "COND_STATUS_CD", "CONDPROP_UNADJ",
  ## Classifications
  "FORTYPCD", "FLDTYPCD", "OWNCD", "OWNGRPCD",
  "STDAGE", "STDSZCD", "HABTYPCD1", "GSSTKCD", "ECOSUBCD",
  ## Phase 2 modifier codes
  "STDORGCD",
  "DSTRBCD1", "DSTRBYR1", "DSTRBCD2", "DSTRBYR2", "DSTRBCD3", "DSTRBYR3",
  "TRTCD1",   "TRTYR1",   "TRTCD2",   "TRTYR2",   "TRTCD3",   "TRTYR3",
  ## Carbon pools (nice to have)
  "CARBON_LITTER", "CARBON_SOIL_ORG"
)

cat("Reading ENTIRE_COND.csv ...\n")
t0 <- Sys.time()
hdr <- fread(opts$cond_csv, nrows = 0L)
select_cols <- intersect(keep_cols, names(hdr))
missing     <- setdiff(keep_cols, select_cols)
if (length(missing) > 0) {
  cat("  Note: missing from source (will skip): ",
      paste(missing, collapse = ", "), "\n", sep = "")
}
cond <- fread(opts$cond_csv, select = select_cols, showProgress = FALSE)
cat(sprintf("  Loaded: %s rows, %d cols in %.1f s\n",
            format(nrow(cond), big.mark = ","), ncol(cond),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (!opts$keep_nonforest) {
  before <- nrow(cond)
  cond <- cond[COND_STATUS_CD == 1]
  cat(sprintf("  Forested only (COND_STATUS_CD=1): %s rows (dropped %s)\n",
              format(nrow(cond), big.mark = ","),
              format(before - nrow(cond), big.mark = ",")))
}

## Dominant condition per plot-measurement by largest CONDPROP_UNADJ
cat("Picking dominant condition per plot-measurement ...\n")
setorder(cond, STATECD, UNITCD, COUNTYCD, PLOT, INVYR, -CONDPROP_UNADJ)
cond_dom <- cond[, .SD[1L], by = .(STATECD, UNITCD, COUNTYCD, PLOT, INVYR)]
cat(sprintf("  Dominant-condition rows: %s\n",
            format(nrow(cond_dom), big.mark = ",")))

cat("Summary of modifier coverage (% non-NA, on dominant conditions):\n")
for (v in c("STDORGCD","DSTRBCD1","DSTRBYR1","TRTCD1","TRTYR1",
            "FORTYPCD","STDAGE","OWNGRPCD","ECOSUBCD")) {
  if (v %in% names(cond_dom)) {
    pct <- 100 * mean(!is.na(cond_dom[[v]]))
    cat(sprintf("  %-14s %6.1f%%\n", v, pct))
  }
}

cat(sprintf("\nWriting %s ...\n", opts$out_path))
saveRDS(cond_dom, opts$out_path, compress = "xz")
sz <- file.info(opts$out_path)$size
cat(sprintf("Done. Output size: %.1f MB\n", sz / 1e6))
