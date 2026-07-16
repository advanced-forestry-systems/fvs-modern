##=============================================================================
## eval_climate_csv_explore.R
##
## Catalog and rank the 40 climate variables in ALL_SI_m.csv against DG.
##
## Steps:
##   1. Load ALL_SI_m.csv (362K rows, FIA + other sources)
##   2. Subset to FIA records (SOURCE == "FIA")
##   3. Join to matched-pairs dataset by SPCD + LAT + LON
##   4. Rank each climate variable by univariate R^2 against dg_obs_a
##   5. Rank conditional R^2 added on top of core covariates
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-15
##=============================================================================

library(data.table)

CSV_FILE  <- "/users/PUOM0008/crsfaaron/SiteIndex/ALL_SI_m.csv"
DATA_FILE <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
OUT_DIR   <- "calibration/output/conus/site_explore"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("== eval_climate_csv_explore.R ==\n\n")

cat("Loading climate CSV ..."); flush.console()
clim <- fread(CSV_FILE)
cat(" done. nrow =", nrow(clim), "ncol =", ncol(clim), "\n")
cat("Columns:", paste(names(clim), collapse=", "), "\n")
cat("SOURCE table:\n")
print(table(clim$SOURCE))

# Subset to FIA only
clim_fia <- clim[SOURCE == "FIA"]
cat("\nFIA rows in climate CSV:", nrow(clim_fia), "\n")
cat("Unique (LAT, LON):", uniqueN(clim_fia[, .(LAT, LON)]), "\n")

# Climate variable candidates (skip metadata + derived ratios for now)
clim_vars <- c("mat","map","gsp","mtcm","mmin","mtwm","mmax","sday","fday",
               "ffp","dd5","gsdd5","d100","dd0","mmindd0","tdiff","adi",
               "adimindd0","dd5mtcm","gspdd5","gspmtcm","gsptd","mapdd5",
               "mapmtcm","maptd","mtcmgsp","mtcmmap","pratio","prdd5",
               "prmtcm","sdi","sdimindd0","tdgsp","tdmap")
clim_vars <- intersect(clim_vars, names(clim_fia))
cat("\nClimate variables available:", length(clim_vars), "\n")

# Load matched-pairs
cat("\nLoading matched-pairs ..."); flush.console()
dat <- as.data.table(readRDS(DATA_FILE))
cat(" done. nrow =", nrow(dat), "\n")
cat("Has LAT/LON?", all(c("LAT","LON") %in% names(dat)), "\n")
cat("Has lat/lon?", all(c("lat","lon") %in% names(dat)), "\n")

# Standardize column names
if ("lat" %in% names(dat) && !"LAT" %in% names(dat)) setnames(dat, "lat", "LAT")
if ("lon" %in% names(dat) && !"LON" %in% names(dat)) setnames(dat, "lon", "LON")

if (!all(c("LAT","LON") %in% names(dat))) {
  cat("WARNING: LAT/LON not in matched-pairs; cannot join.\n")
  cat("Available coordinate-like cols:",
      paste(grep("lat|lon|LAT|LON|coord", names(dat), value=TRUE), collapse=", "), "\n")
  quit(status = 1)
}

# Round to 5 decimal places to allow exact match (FIA fuzzed coordinates)
dat[, lat_r := round(LAT, 5)]
dat[, lon_r := round(LON, 5)]
clim_fia[, lat_r := round(LAT, 5)]
clim_fia[, lon_r := round(LON, 5)]

# Try species-aware first (SPCD + lat + lon)
if ("SPCD" %in% names(clim_fia) && "SPCD" %in% names(dat)) {
  cat("\nJoining by SPCD + lat_r + lon_r ...\n")
  setkey(clim_fia, SPCD, lat_r, lon_r)
  setkey(dat, SPCD, lat_r, lon_r)
  joined <- clim_fia[dat, nomatch = NULL,
                     on = c("SPCD","lat_r","lon_r")]
  cat("  joined rows:", nrow(joined), "\n")
}

# Fallback: lat + lon only
if (!exists("joined") || nrow(joined) < 1000) {
  cat("\nFalling back to lat_r + lon_r join ...\n")
  # Pick one row per (lat, lon) from climate CSV
  clim_unique <- unique(clim_fia, by = c("lat_r","lon_r"))
  setkey(clim_unique, lat_r, lon_r)
  setkey(dat, lat_r, lon_r)
  joined <- clim_unique[, c("lat_r","lon_r", clim_vars), with=FALSE][dat,
                          nomatch = NULL, on = c("lat_r","lon_r")]
  cat("  joined rows:", nrow(joined), "\n")
}

if (nrow(joined) < 10000) {
  cat("Join failed. Aborting.\n"); quit(status = 1)
}

# Compute dg_obs_a + filter
joined[, dg_obs_a := (DBH2 - DBH1) / YEARS]
joined <- joined[is.finite(dg_obs_a) & dg_obs_a > -0.5 & dg_obs_a < 5 &
                  is.finite(DBH1) & DBH1 >= 2.54 & is.finite(CR1) & CR1 > 0]
cat("\nClean rows for analysis:", nrow(joined), "\n\n")

# Core covariates
joined[, ln_dbh        := log(DBH1)]
joined[, ln_cr_adj     := log((CR1 + 0.2) / 1.2)]
joined[, ln_bal_sw_adj := log(pmax(BAL_SW1, 0) + 0.01)]

# Univariate R^2 per climate variable
cat("=== Univariate R^2 for dg_obs_a ===\n")
results <- data.table()
for (v in clim_vars) {
  if (!v %in% names(joined)) next
  y <- joined$dg_obs_a
  x <- joined[[v]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5000) next
  fit <- lm(y[ok] ~ x[ok])
  r2 <- summary(fit)$r.squared
  results <- rbind(results,
                   data.table(var = v, n = sum(ok), R2 = r2,
                              sign = sign(coef(fit)[2])))
}
setorder(results, -R2)
print(results)
fwrite(results, file.path(OUT_DIR, "climate_univariate_R2.csv"))

# Conditional R^2 added after core
cat("\n=== Conditional R^2 added (over ln_dbh + ln_cr_adj + ln_bal_sw_adj) ===\n")
core_fit <- lm(dg_obs_a ~ ln_dbh + ln_cr_adj + ln_bal_sw_adj, data = joined)
core_r2 <- summary(core_fit)$r.squared
cat("Core R^2 =", round(core_r2, 4), "\n\n")

cond_results <- data.table()
for (v in clim_vars) {
  if (!v %in% names(joined)) next
  f <- as.formula(paste0("dg_obs_a ~ ln_dbh + ln_cr_adj + ln_bal_sw_adj + ", v))
  fit <- tryCatch(lm(f, data = joined), error = function(e) NULL)
  if (is.null(fit)) next
  r2 <- summary(fit)$r.squared
  cond_results <- rbind(cond_results,
                        data.table(var = v, R2_added = r2 - core_r2,
                                   coef = coef(fit)[v]))
}
setorder(cond_results, -R2_added)
print(cond_results)
fwrite(cond_results, file.path(OUT_DIR, "climate_conditional_R2.csv"))

cat("\nDone. Outputs in:", OUT_DIR, "\n")
