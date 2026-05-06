## =============================================================================
## 42_build_species_traits.R
##
## Assemble a SPCD-keyed functional trait table for FVS-CONUS, for use both
## as covariates in HT-DBH and as a trait-informed prior on the species
## random intercept in DG, HG, HCB, and mortality:
##
##   z_sp ~ Normal( W . traits_sp, sigma_sp )
##
## Sources (in order of precedence when multiple exist):
##
## 1. FIA REF_SPECIES:
##      WOOD_SPGR_GREENVOL_DRYWT    -> wood_specific_gravity (WD / SG)
##      SHADE_TOLERANCE (1-5)        -> shade_tolerance_num
##      EVERGREEN (0/1)              -> evergreen
##      DECIDUOUS (0/1)              -> deciduous
##      GROWTH_HABIT_CD              -> growth_habit
##
## 2. leaf_longevity_by_spcd.csv (Hallik et al. 2009, Reich et al. 1999,
##    Niinemets & Lukjanova 2002, Sargent 1897, and others; compiled by AW):
##      LL_months                     -> leaf_longevity_months
##
## 3. Chave et al. 2009 global WD database:
##    species lookup keyed to genus+species; fills WD where FIA is NA.
##      wood_specific_gravity        <- Chave WD (g/cm^3)
##
## 4. TRY plant trait database export (user must provide a TRY TraitID slice
##    for SLA mm^2/mg and seed dry mass mg; see README in calibration/traits/):
##      SLA (TraitID 3117)            -> sla_mm2_mg
##      Seed mass (TraitID 26)        -> seed_mass_mg
##      Plant height at maturity (TraitID 3106) -> plant_height_tr_m
##
## 5. National Champion Tree Program (nationalchampiontree.org):
##    per-species maximum DBH and height among registered champions;
##    scraped or manually compiled to calibration/traits/champion_trees.csv.
##      max_ht_champ_m, max_dbh_champ_cm
##
## 6. FIA empirical maxima: 99th percentile HT1 and DBH1 per SPCD from
##    conus_remeasurement_pairs_metric.rds (computed on Cardinal).
##      max_ht_fia_p99_m, max_dbh_fia_p99_cm
##
## 7. Bechtold (2003) / Gonzalez-Benecke maximum crown-width allometry:
##    per-species max crown width at a reference DBH; fills
##      max_crown_width_m   (optional)
##
## Output: calibration/traits/species_traits.rds and .csv.
##
## Usage:
##   Rscript scripts/42_build_species_traits.R
##   Rscript scripts/42_build_species_traits.R --out calibration/traits/
##
## Author: A. Weiskittel, 2026-04-15
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--out",          type = "character",
              default = "calibration/traits"),
  make_option("--ref_species",  type = "character",
              default = "calibration/traits/REF_SPECIES.csv",
              help = "FIA REF_SPECIES CSV (from FIA DataMart)"),
  make_option("--leaf_long",    type = "character",
              default = "calibration/traits/leaf_longevity_by_spcd.csv"),
  make_option("--chave_wd",     type = "character",
              default = "calibration/traits/chave_wood_density.csv",
              help = "Optional Chave et al. 2009 global WD lookup"),
  make_option("--try_traits",   type = "character",
              default = "calibration/traits/try_trait_export.csv",
              help = "Optional TRY export with SPCD + TraitID + Value"),
  make_option("--champion",     type = "character",
              default = "calibration/traits/champion_trees.csv",
              help = "Optional National Champion Tree Program table"),
  make_option("--burns_honkala", type = "character",
              default = "calibration/traits/burns_honkala_shade_tolerance.csv",
              help = "Burns & Honkala (1990) shade tolerance fallback"),
  make_option("--niinemets",    type = "character",
              default = "calibration/traits/niinemets_shade_tolerance.csv",
              help = "Niinemets & Valladares (2006) continuous shade tolerance"),
  make_option("--potter_vcc",   type = "character",
              default = "calibration/traits/potter_vcc_scores.csv",
              help = "Potter et al. (2017) climate vulnerability scores"),
  make_option("--fia_pairs",    type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds",
              help = "For 99th percentile HT/DBH empirical maxima per SPCD"),
  make_option("--verbose",      action = "store_true", default = TRUE)
)) |> parse_args()

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path, label) {
  if (!file.exists(path)) {
    message("  [skip] ", label, " not found at ", path)
    return(NULL)
  }
  message("  [read] ", label, " from ", path)
  fread(path)
}

## -----------------------------------------------------------------------------
## 1. FIA REF_SPECIES
## -----------------------------------------------------------------------------
message("\n== FIA REF_SPECIES ==")
ref <- safe_read_csv(opts$ref_species, "REF_SPECIES")
if (!is.null(ref)) {
  setnames(ref, toupper(names(ref)))

  ## Core columns always present in FIADB REF_SPECIES exports
  ref <- ref[, .(
    SPCD,
    GENUS,
    SPECIES,
    COMMON_NAME,
    SFTWD_HRDWD           = fifelse(is.na(SFTWD_HRDWD), "U", SFTWD_HRDWD),
    wood_specific_gravity  = as.numeric(WOOD_SPGR_GREENVOL_DRYWT)
  )]

  ## SHADE_TOLERANCE may or may not be in the export (the ENTIRE_REF_SPECIES
  ## from FIA DataMart omits it). If present, use it; otherwise leave NA so
  ## the trait matrix builder median-imputes later.
  if ("SHADE_TOLERANCE" %in% names(fread(opts$ref_species, nrows = 0))) {
    st <- fread(opts$ref_species, select = c("SPCD", "SHADE_TOLERANCE"))
    ref <- merge(ref, st[, .(SPCD, shade_tolerance_num = as.numeric(SHADE_TOLERANCE))],
                 by = "SPCD", all.x = TRUE)
  } else {
    message("  SHADE_TOLERANCE column not in export; will be NA (median-imputed downstream)")
    ref[, shade_tolerance_num := NA_real_]
  }

  ## EVERGREEN / DECIDUOUS may also be absent. Derive from SFTWD_HRDWD as a
  ## reasonable approximation: S = softwood ~ evergreen, H = hardwood ~ deciduous.
  if ("EVERGREEN" %in% names(fread(opts$ref_species, nrows = 0))) {
    ev <- fread(opts$ref_species, select = c("SPCD", "EVERGREEN", "DECIDUOUS"))
    ref <- merge(ref, ev[, .(SPCD,
                              evergreen = as.integer(EVERGREEN),
                              deciduous = as.integer(DECIDUOUS))],
                 by = "SPCD", all.x = TRUE)
  } else {
    message("  EVERGREEN/DECIDUOUS columns not in export; deriving from SFTWD_HRDWD")
    ref[, evergreen := fifelse(SFTWD_HRDWD == "S", 1L, 0L)]
    ref[, deciduous := fifelse(SFTWD_HRDWD == "H", 1L, 0L)]
  }

  ## Binary softwood indicator (useful even when shade_tolerance is missing)
  ref[, softwood := fifelse(SFTWD_HRDWD == "S", 1L, 0L)]

  message("  rows: ", nrow(ref), "  with WD: ",
          sum(!is.na(ref$wood_specific_gravity)),
          "  shade_tol present: ",
          sum(!is.na(ref$shade_tolerance_num)))
} else {
  stop("REF_SPECIES is required; place it at ", opts$ref_species,
       " (FIADB data dictionary -> REF_SPECIES table).")
}

## -----------------------------------------------------------------------------
## 1b. Niinemets & Valladares (2006) shade tolerance (primary fallback)
## Continuous 1-5 scale from doi.org/10.6084/m9.figshare.c.3309258.
## Higher quality than Burns & Honkala due to continuous measurements.
## Precedence: REF_SPECIES > Niinemets > Burns & Honkala
## -----------------------------------------------------------------------------
message("\n== Niinemets & Valladares shade tolerance ==")
nv <- safe_read_csv(opts$niinemets, "Niinemets & Valladares shade tolerance")
if (!is.null(nv)) {
  setnames(nv, toupper(names(nv)))
  nv <- nv[, .(SPCD = as.integer(SPCD),
               shade_tol_nv = as.numeric(SHADE_TOLERANCE))]
  nv <- nv[!is.na(SPCD) & !is.na(shade_tol_nv)]
  nv <- nv[, .(shade_tol_nv = mean(shade_tol_nv)), by = SPCD]
  ref <- merge(ref, nv, by = "SPCD", all.x = TRUE)
  filled <- is.na(ref$shade_tolerance_num) & !is.na(ref$shade_tol_nv)
  ref[filled, shade_tolerance_num := shade_tol_nv]
  message("  filled shade_tolerance from Niinemets for ", sum(filled), " species")
  ref[, shade_tol_nv := NULL]
} else {
  message("  [skip] Niinemets lookup not found")
}

## Burns & Honkala (1990) fills any remaining gaps after Niinemets
message("\n== Burns & Honkala shade tolerance fallback ==")
bh <- safe_read_csv(opts$burns_honkala, "Burns & Honkala shade tolerance")
if (!is.null(bh)) {
  bh <- bh[, .(SPCD = as.integer(SPCD),
               shade_tol_bh = as.numeric(shade_tolerance_num))]
  bh <- bh[!is.na(SPCD) & !is.na(shade_tol_bh)]
  bh <- bh[, .(shade_tol_bh = mean(shade_tol_bh)), by = SPCD]
  ref <- merge(ref, bh, by = "SPCD", all.x = TRUE)
  filled <- is.na(ref$shade_tolerance_num) & !is.na(ref$shade_tol_bh)
  ref[filled, shade_tolerance_num := shade_tol_bh]
  message("  filled shade_tolerance from B&H for ", sum(filled), " additional species")
  ref[, shade_tol_bh := NULL]
} else {
  message("  [skip] Burns & Honkala lookup not found")
}

## -----------------------------------------------------------------------------
## 1c. Potter et al. (2017) climate vulnerability scores (CAPTURE framework)
## doi.org/10.1007/s11056-017-9569-5
## Composite vulnerability_score plus sub-components: climate_exposure,
## sensitivity, low_adaptive_capacity. Higher score = more vulnerable.
## -----------------------------------------------------------------------------
message("\n== Potter et al. (2017) climate vulnerability ==")
vcc <- safe_read_csv(opts$potter_vcc, "Potter VCC scores")
if (!is.null(vcc)) {
  setnames(vcc, toupper(names(vcc)))
  vcc <- vcc[, .(SPCD             = as.integer(SPCD),
                 vulnerability_score = as.numeric(VULNERABILITY_SCORE),
                 climate_exposure    = as.numeric(CLIMATE_EXPOSURE),
                 sensitivity         = as.numeric(SENSITIVITY),
                 low_adaptive_cap    = as.numeric(LOW_ADAPTIVE_CAPACITY))]
  vcc <- vcc[!is.na(SPCD)]
  ## Some species appear in multiple subspecies; take the mean per SPCD
  vcc <- vcc[, lapply(.SD, mean, na.rm = TRUE), by = SPCD]
  ref <- merge(ref, vcc, by = "SPCD", all.x = TRUE)
  message("  rows merged: ", sum(!is.na(ref$vulnerability_score)),
          " species with VCC scores")
} else {
  ref[, `:=`(vulnerability_score = NA_real_,
             climate_exposure    = NA_real_,
             sensitivity         = NA_real_,
             low_adaptive_cap    = NA_real_)]
}

## -----------------------------------------------------------------------------
## 2. Leaf longevity (Aaron's compiled table, Apr 2026)
## -----------------------------------------------------------------------------
message("\n== Leaf longevity ==")
ll <- safe_read_csv(opts$leaf_long, "leaf longevity")
if (!is.null(ll)) {
  setnames(ll, 1:2, c("SPCD", "leaf_longevity_months"))
  ll <- ll[, .(SPCD = as.integer(SPCD),
               leaf_longevity_months = as.numeric(leaf_longevity_months))]
  ll <- ll[!is.na(SPCD) & !is.na(leaf_longevity_months)]
  ll <- ll[, .(leaf_longevity_months = mean(leaf_longevity_months)), by = SPCD]
  message("  rows: ", nrow(ll))
  ref <- merge(ref, ll, by = "SPCD", all.x = TRUE)
} else {
  ref[, leaf_longevity_months := NA_real_]
}

## -----------------------------------------------------------------------------
## 3. Chave WD fallback
## -----------------------------------------------------------------------------
message("\n== Chave WD fallback ==")
chave <- safe_read_csv(opts$chave_wd, "Chave WD")
if (!is.null(chave)) {
  setnames(chave, toupper(names(chave)))
  chave <- chave[, .(GENUS,
                     SPECIES,
                     WD_chave = as.numeric(WOOD_DENSITY_G_CM3))]
  ref <- merge(ref, chave, by = c("GENUS", "SPECIES"), all.x = TRUE)
  filled <- is.na(ref$wood_specific_gravity) & !is.na(ref$WD_chave)
  ref[filled, wood_specific_gravity := WD_chave]
  message("  filled WD from Chave for ", sum(filled), " species")
  ref[, WD_chave := NULL]
}

## -----------------------------------------------------------------------------
## 4. TRY traits (SLA, seed mass, plant height)
## -----------------------------------------------------------------------------
message("\n== TRY traits ==")
tr <- safe_read_csv(opts$try_traits, "TRY")
if (!is.null(tr)) {
  ## Expected long format: SPCD, TraitID, TraitName, Value, StdValue
  setnames(tr, toupper(names(tr)))
  tr_w <- dcast(tr[, .(SPCD, TraitID, StdValue)],
                SPCD ~ TraitID,
                value.var = "StdValue",
                fun.aggregate = mean, na.rm = TRUE)
  setnames(tr_w,
           old = intersect(c("3117", "26", "3106"), names(tr_w)),
           new = c("sla_mm2_mg", "seed_mass_mg", "plant_height_tr_m")[
             match(intersect(c("3117","26","3106"), names(tr_w)),
                   c("3117","26","3106"))])
  ref <- merge(ref, tr_w, by = "SPCD", all.x = TRUE)
  message("  rows merged from TRY: ", nrow(tr_w))
} else {
  ref[, `:=`(sla_mm2_mg        = NA_real_,
             seed_mass_mg      = NA_real_,
             plant_height_tr_m = NA_real_)]
}

## -----------------------------------------------------------------------------
## 5. National Champion Tree Program
## -----------------------------------------------------------------------------
message("\n== Champion tree records ==")
ch <- safe_read_csv(opts$champion, "Champion trees")
if (!is.null(ch)) {
  setnames(ch, toupper(names(ch)))
  ch <- ch[, .(
    SPCD,
    max_ht_champ_m   = as.numeric(HEIGHT_M),
    max_dbh_champ_cm = as.numeric(DBH_CM),
    max_crown_champ_m = as.numeric(CROWN_WIDTH_M)
  )]
  ch <- ch[!is.na(SPCD)]
  ch <- ch[, lapply(.SD, max, na.rm = TRUE), by = SPCD]
  for (cc in setdiff(names(ch), "SPCD")) {
    ch[[cc]][is.infinite(ch[[cc]])] <- NA_real_
  }
  ref <- merge(ref, ch, by = "SPCD", all.x = TRUE)
  message("  rows: ", nrow(ch))
} else {
  ref[, `:=`(max_ht_champ_m    = NA_real_,
             max_dbh_champ_cm  = NA_real_,
             max_crown_champ_m = NA_real_)]
}

## -----------------------------------------------------------------------------
## 6. FIA empirical 99th percentile HT/DBH per SPCD
## -----------------------------------------------------------------------------
message("\n== FIA empirical maxima (99th percentile) ==")
if (file.exists(opts$fia_pairs)) {
  fia <- readRDS(opts$fia_pairs)
  setDT(fia)
  fia_max <- fia[!is.na(DBH1) & !is.na(HT1), .(
    max_dbh_fia_p99_cm = quantile(DBH1, 0.99, na.rm = TRUE),
    max_ht_fia_p99_m   = quantile(HT1,  0.99, na.rm = TRUE),
    n_fia              = .N
  ), by = SPCD]
  ref <- merge(ref, fia_max, by = "SPCD", all.x = TRUE)
  message("  SPCDs with FIA p99: ", nrow(fia_max))
} else {
  message("  [skip] ", opts$fia_pairs, " not found; run on Cardinal")
  ref[, `:=`(max_dbh_fia_p99_cm = NA_real_,
             max_ht_fia_p99_m   = NA_real_,
             n_fia              = NA_integer_)]
}

## -----------------------------------------------------------------------------
## 7. Harmonize a "max height" column (consensus across sources)
##    Precedence: champion tree max height > TRY plant_height > FIA p99
## -----------------------------------------------------------------------------
ref[, max_ht_m := pmax(max_ht_champ_m, plant_height_tr_m, max_ht_fia_p99_m,
                        na.rm = TRUE)]
ref[, max_dbh_cm := pmax(max_dbh_champ_cm, max_dbh_fia_p99_cm, na.rm = TRUE)]

## -----------------------------------------------------------------------------
## 8. Missingness summary
## -----------------------------------------------------------------------------
message("\n== Trait completeness (n / %) ==")
trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "sla_mm2_mg", "seed_mass_mg",
                "plant_height_tr_m", "max_ht_m", "max_dbh_cm",
                "max_crown_champ_m",
                "vulnerability_score", "climate_exposure",
                "sensitivity", "low_adaptive_cap")
for (cc in trait_cols) {
  have <- sum(!is.na(ref[[cc]]))
  message(sprintf("  %-28s %6d  (%4.1f%%)",
                  cc, have, 100 * have / nrow(ref)))
}

## -----------------------------------------------------------------------------
## 9. Write outputs
## -----------------------------------------------------------------------------
out_rds <- file.path(opts$out, "species_traits.rds")
out_csv <- file.path(opts$out, "species_traits.csv")
saveRDS(ref, out_rds)
fwrite(ref, out_csv)
message("\nWrote: ", out_rds)
message("Wrote: ", out_csv)
message("Total SPCDs: ", nrow(ref))
