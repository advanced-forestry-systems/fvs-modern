## =============================================================================
## 40_enrich_cond_modifiers.R
##
## FVS-CONUS Phase 2 staging: attach condition-level modifiers to the
## remeasurement pairs file.
##
## Joins FIA COND table on (PLT_CN, CONDID) to pull:
##
##   STDORGCD   stand origin (natural vs plantation)
##   DSTRBCD1-3 disturbance codes, first three agents observed
##   DSTRBYR1-3 year of disturbance for each
##   TRTCD1-3   treatment codes, first three treatments
##   TRTYR1-3   year of treatment for each
##   FORTYPCD   forest type (operational classification)
##   FLDTYPCD   field-called type (pre-algorithmic classification)
##   OWNCD      ownership class
##   OWNGRPCD   ownership group
##   STDAGE     stand age
##   STDSZCD    stand-size class
##   HABTYPCD1  habitat-type classification 1
##   GSSTKCD    growing-stock stocking
##   CARBON_LITTER
##   CARBON_SOIL_ORG
##
## Output: conus_remeasurement_pairs_metric_cond.rds. Non-destructive;
## leaves the canonical metric RDS untouched. Designed to run after
## Phase 1 fits have completed and Phase 2 planning begins.
##
## Staging notes:
##
## 1. STDORGCD as a plantation vs natural indicator lets DG, HG, mortality
##    models absorb plantation effects via an indicator term rather than
##    through a separate variant.
## 2. DSTRBCD modifiers (fire = 30, insect = 10-19, disease = 20-29, wind =
##    50, harvest = 80, other = 90+) allow the mortality model to separate
##    harvest from natural mortality with better resolution than the
##    current AGENTCD2 < 80 filter.
## 3. TRTCD modifiers (cutting = 10, natural regen = 20, artificial regen =
##    30, site preparation = 40, other = 50) pair with STDORGCD to let the
##    DG model respond to recent treatment.
##
## USAGE:
##   Rscript calibration/R/40_enrich_cond_modifiers.R \
##       --pairs calibration/data/conus_remeasurement_pairs_metric.rds \
##       --cond  calibration/data/fia_cond.rds \
##       --out   calibration/data/conus_remeasurement_pairs_metric_cond.rds
##
## The COND table can be compiled from the same state zips used by
## 30_build_conus_dataset.R; add a step that reads COND.csv from each state
## zip in parallel and saves a consolidated fia_cond.rds.
##
## Author: A. Weiskittel, 2026-04-22
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--pairs", type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds",
              dest    = "pairs_path"),
  make_option("--cond",  type = "character",
              default = "calibration/data/fia_cond.rds",
              dest    = "cond_path"),
  make_option("--out",   type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond.rds",
              dest    = "out_path")
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS Phase 2 COND-level enrichment\n")
cat("Pairs :", opts$pairs_path, "\n")
cat("COND  :", opts$cond_path, "\n")
cat("Output:", opts$out_path, "\n")
cat("==============================================================\n")

stopifnot(file.exists(opts$pairs_path))

pairs <- readRDS(opts$pairs_path)
setDT(pairs)
cat(sprintf("Pairs : %s rows, %d cols\n",
            format(nrow(pairs), big.mark = ","), ncol(pairs)))

## ---- 1. Require or build FIA COND ----------------------------------------

if (!file.exists(opts$cond_path)) {
  stop(sprintf(paste0(
    "COND RDS not found at %s.\n",
    "Build it first by reading COND.csv from each state zip under ",
    "calibration/data/fia_zips/.\n",
    "A template is in scripts/40_build_cond_from_zips.R (to be written)."
  ), opts$cond_path))
}

cond <- readRDS(opts$cond_path)
setDT(cond)
cat(sprintf("COND  : %s rows, %d cols\n",
            format(nrow(cond), big.mark = ","), ncol(cond)))

## ---- 2. Column selection --------------------------------------------------

keep_cols <- c(
  "PLT_CN", "CONDID",
  "STDORGCD",
  "DSTRBCD1", "DSTRBYR1", "DSTRBCD2", "DSTRBYR2", "DSTRBCD3", "DSTRBYR3",
  "TRTCD1",   "TRTYR1",   "TRTCD2",   "TRTYR2",   "TRTCD3",   "TRTYR3",
  "FORTYPCD", "FLDTYPCD",
  "OWNCD", "OWNGRPCD",
  "STDAGE", "STDSZCD",
  "HABTYPCD1",
  "GSSTKCD",
  "CARBON_LITTER", "CARBON_SOIL_ORG"
)
avail <- intersect(keep_cols, names(cond))
missing <- setdiff(keep_cols, avail)
if (length(missing) > 0) {
  cat("  COND missing columns (will skip):", paste(missing, collapse = ", "), "\n")
}
cond_sel <- cond[, ..avail]

## ---- 3. Build join keys matching pairs ------------------------------------
##
## Pairs file uses PLT_CN1 and CONDID1 for the t1 condition. We attach the
## t1-condition covariates as a baseline. For Phase 2 we can also attach the
## t2 condition to capture mid-interval events.

if (!"PLT_CN1" %in% names(pairs)) {
  stop("pairs RDS missing PLT_CN1 (t1 plot condition key).")
}

setnames(cond_sel, "PLT_CN", "PLT_CN1", skip_absent = TRUE)
if ("CONDID" %in% names(cond_sel)) {
  setnames(cond_sel, "CONDID", "CONDID1")
}

## Rename all non-key columns with a _cond1 suffix so they do not collide
## with any existing pair columns.
key_cols <- c("PLT_CN1", "CONDID1")
suffix_cols <- setdiff(names(cond_sel), key_cols)
new_names   <- paste0(suffix_cols, "_cond1")
setnames(cond_sel, suffix_cols, new_names)

cat(sprintf("Joining %d COND columns to pairs on (PLT_CN1, CONDID1) ...\n",
            length(new_names)))

if ("CONDID1" %in% names(pairs)) {
  pairs_out <- merge(pairs, cond_sel, by = key_cols, all.x = TRUE)
} else {
  ## Fallback if CONDID1 was not carried through. Join on PLT_CN1 alone and
  ## accept one-to-many row expansion, warning at the end.
  cat("  CONDID1 not in pairs; joining on PLT_CN1 only (may duplicate rows).\n")
  pairs_out <- merge(pairs, cond_sel[, ..c("PLT_CN1", new_names)],
                     by = "PLT_CN1", all.x = TRUE)
}

cat(sprintf("  Joined rows: %s (pre-join %s)\n",
            format(nrow(pairs_out), big.mark = ","),
            format(nrow(pairs),     big.mark = ",")))

## ---- 4. Derive compact modifier flags -------------------------------------
##
## Simple indicators for downstream Stan models. These sit alongside the
## raw codes; keep both.

if ("STDORGCD_cond1" %in% names(pairs_out)) {
  pairs_out[, is_plantation := fifelse(STDORGCD_cond1 == 1, 1L, 0L, na = NA_integer_)]
}
if ("DSTRBCD1_cond1" %in% names(pairs_out)) {
  pairs_out[, had_fire_t1     := as.integer(DSTRBCD1_cond1 == 30)]
  pairs_out[, had_insect_t1   := as.integer(DSTRBCD1_cond1 %in% 10:19)]
  pairs_out[, had_disease_t1  := as.integer(DSTRBCD1_cond1 %in% 20:29)]
  pairs_out[, had_wind_t1     := as.integer(DSTRBCD1_cond1 == 50)]
  pairs_out[, had_harvest_t1  := as.integer(DSTRBCD1_cond1 == 80)]
}
if ("TRTCD1_cond1" %in% names(pairs_out)) {
  pairs_out[, had_cutting_t1   := as.integer(TRTCD1_cond1 == 10)]
  pairs_out[, had_site_prep_t1 := as.integer(TRTCD1_cond1 == 40)]
}

cat("\nModifier coverage (% non-NA):\n")
for (v in c("is_plantation", "had_fire_t1", "had_insect_t1", "had_harvest_t1",
            "had_cutting_t1")) {
  if (v %in% names(pairs_out)) {
    pct <- 100 * mean(!is.na(pairs_out[[v]]))
    mn  <- mean(pairs_out[[v]], na.rm = TRUE)
    cat(sprintf("  %-16s coverage=%.1f%%  mean=%.4f\n", v, pct, mn))
  }
}

cat(sprintf("\nWriting enriched pairs to %s ...\n", opts$out_path))
saveRDS(pairs_out, opts$out_path)
cat(sprintf("Done. %s rows, %d cols\n",
            format(nrow(pairs_out), big.mark = ","), ncol(pairs_out)))
