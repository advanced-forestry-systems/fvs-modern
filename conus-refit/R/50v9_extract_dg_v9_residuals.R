##=============================================================================
## 50v9_extract_dg_v9_residuals.R
##
## Modifier-residual extractor for DG_Kuehne v9 (and v8) base fits.
##
## Key fixes relative to legacy 50_extract_base_residuals.R:
##   1. Filters dg_obs_a > 0.01 BEFORE log() — avoids σ_resid = 21 bug
##      caused by log(pmax(., 0.001)) injecting -6.9 outliers.
##   2. Reads v9 meta format (sp_levels at top level, not under prep_meta).
##   3. Knows the v9 piecewise BGI / mapdd5 + FT random effect eta formula.
##   4. Honors the same filter set as the v9 base fit, so residuals
##      correspond to the actual fitted observations.
##
## CLI:
##   --fit       v9 fit RDS
##   --meta      v9 meta RDS (defaults to fit path with _fit -> _meta)
##   --variant   v8 (BGI piecewise) | v9 (mapdd5 piecewise)
##   --pairs     matched-pairs data
##   --traits    species traits
##   --out       residual bundle RDS
##   --subsample N
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-15
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

FIT_FILE  <- get_arg("fit")
META_FILE <- get_arg("meta", sub("_fit\\.rds$", "_meta.rds", FIT_FILE %||% ""))
VARIANT   <- get_arg("variant", "v8")
PAIRS_FILE  <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")
OUT_FILE  <- get_arg("out", paste0("calibration/output/conus/dg_kue/", VARIANT, "/dg_kuehne_", VARIANT, "_residuals.rds"))
SUBSAMPLE <- as.integer(get_arg("subsample", "0"))

stopifnot(!is.null(FIT_FILE), file.exists(FIT_FILE))
stopifnot(file.exists(META_FILE))
stopifnot(VARIANT %in% c("v8", "v9"))

cat("== 50v9_extract_dg_v9_residuals.R ==\n")
cat("  fit    :", FIT_FILE, "\n")
cat("  meta   :", META_FILE, "\n")
cat("  variant:", VARIANT, "\n")
cat("  out    :", OUT_FILE, "\n\n")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

cat("Loading fit + meta + pairs ..."); flush.console()
fit   <- readRDS(FIT_FILE)
meta  <- readRDS(META_FILE)
pairs <- as.data.table(readRDS(PAIRS_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done\n")

sp_levels <- meta$sp_levels
L1_levels <- meta$L1_levels
L2_levels <- meta$L2_levels
L3_levels <- meta$L3_levels
FT_levels <- meta$FT_levels

cat(sprintf("Levels: sp=%d L1=%d L2=%d L3=%d FT=%d\n",
            length(sp_levels), length(L1_levels), length(L2_levels),
            length(L3_levels), length(FT_levels)))

# Recreate derived columns and apply same filters as base driver
pairs[, dg_obs_a       := (DBH2 - DBH1) / YEARS]
pairs[, sqrt_years     := sqrt(YEARS)]
pairs[, ln_dbh         := log(DBH1)]
pairs[, ln_cr_adj      := log((CR1 + 0.2) / 1.2)]
pairs[, ln_bal_sw_adj  := log(BAL_SW1 + 0.01)]
pairs[, rd_additive    := sdi_additive1 / SDImax_brms]
pairs[, sdi_complexity := sdi_additive1 / pmax(SDI1, 1.0)]

# For v9 we also need mapdd5 — join from climate CSV
if (VARIANT == "v9") {
  clim_csv <- "/users/PUOM0008/crsfaaron/SiteIndex/ALL_SI_m.csv"
  cat("Loading climate CSV for mapdd5 join ..."); flush.console()
  clim <- fread(clim_csv)
  clim_fia <- clim[SOURCE == "FIA", .(SPCD, LAT, LON, mapdd5)]
  clim_fia[, lat_r := round(LAT, 5)]
  clim_fia[, lon_r := round(LON, 5)]
  clim_fia[, c("LAT","LON") := NULL]
  pairs[, lat_r := round(LAT, 5)]
  pairs[, lon_r := round(LON, 5)]
  setkey(clim_fia, SPCD, lat_r, lon_r)
  setkey(pairs, SPCD, lat_r, lon_r)
  pairs <- clim_fia[pairs, on = c("SPCD","lat_r","lon_r")]
  cat(" done\n")
}

# Apply v8/v9 filters (same as base driver) — **KEY FIX**: dg_obs_a > 0.01
filt <- is.finite(pairs$DBH1) & pairs$DBH1 >= 2.54 & is.finite(pairs$DBH2) &
        is.finite(pairs$CR1) & pairs$CR1 > 0 & pairs$CR1 <= 1.0 &
        is.finite(pairs$YEARS) & pairs$YEARS >= 1 & pairs$YEARS <= 20 &
        !is.na(pairs$EPA_L1_CODE) & !is.na(pairs$EPA_L2_CODE) & !is.na(pairs$EPA_L3_CODE) &
        pairs$EPA_L1_CODE != "" & pairs$EPA_L2_CODE != "" & pairs$EPA_L3_CODE != "" &
        pairs$TREESTATUS1 == 1 & pairs$TREESTATUS2 == 1 &
        is.finite(pairs$BAL_SW1) & pairs$BAL_SW1 >= 0 &
        is.finite(pairs$BAL_HW1) & pairs$BAL_HW1 >= 0 &
        is.finite(pairs$rd_additive) & pairs$rd_additive > 0 & pairs$rd_additive < 3.0 &
        is.finite(pairs$sdi_complexity) & pairs$sdi_complexity > 0 & pairs$sdi_complexity < 10 &
        is.finite(pairs$BA1) & pairs$BA1 >= 0 &
        !is.na(pairs$FORTYPCD_cond1) & pairs$FORTYPCD_cond1 > 0 &
        # KEY FIX: require positive growth strictly above floor
        pairs$dg_obs_a > 0.01 & pairs$dg_obs_a < 5.0
if (VARIANT == "v8") filt <- filt & is.finite(pairs$bgi)
if (VARIANT == "v9") filt <- filt & is.finite(pairs$mapdd5)

d <- pairs[filt]
cat("After filters:", nrow(d), "rows\n")

# Keep only species that match base fit
d <- d[SPCD %in% sp_levels]
cat("After species match:", nrow(d), "rows\n")

# Build indices that match base
d[, sp_idx := match(SPCD, sp_levels)]
d[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
d[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
d[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
d[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]

# Drop rows where any level is NA (level not present in base)
d <- d[!is.na(sp_idx) & !is.na(L1_idx) & !is.na(L2_idx) & !is.na(L3_idx) & !is.na(FT_idx)]
cat("After valid level match:", nrow(d), "rows\n")

if (SUBSAMPLE > 0 && SUBSAMPLE < nrow(d)) {
  set.seed(42)
  d <- d[sample(.N, SUBSAMPLE)]
  cat("Subsampled to:", nrow(d), "rows\n")
}

# Trait matrix and softwood (for v8/v9 with mapdd5 x softwood interaction)
trait_cols <- meta$trait_cols
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
softwood_by_sp <- traits_sub$softwood
softwood_by_sp[is.na(softwood_by_sp)] <- 0
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}
softwood_per_tree <- softwood_by_sp[d$sp_idx]
softwood_per_tree_c <- softwood_per_tree - mean(softwood_per_tree)

# Posterior-mean accessors
get_mean <- function(v) {
  out <- tryCatch(fit$draws(variables = v, format = "draws_matrix"),
                   error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  mean(as.numeric(out))
}
vec_mean <- function(base, n) {
  out <- tryCatch(fit$draws(variables = base, format = "draws_matrix"),
                   error = function(e) NULL)
  if (is.null(out)) return(rep(0, n))
  colMeans(as.matrix(out))
}

cat("Pulling posterior means ..."); flush.console()
b0  <- get_mean("b0")
b1  <- get_mean("b1"); b2 <- get_mean("b2"); b3 <- get_mean("b3")
b4  <- get_mean("b4"); b5 <- get_mean("b5")
b6  <- get_mean("b6")
b7  <- get_mean("b7"); b8 <- get_mean("b8")
b9a <- get_mean("b9a"); b9b <- get_mean("b9b")
b11 <- get_mean("b11")
b12 <- get_mean("b12"); b13 <- get_mean("b13"); b14 <- get_mean("b14"); b15 <- get_mean("b15")

trait_effect       <- vec_mean("trait_effect", length(sp_levels))
species_site_slope <- vec_mean("species_site_slope", length(sp_levels))
z_sp     <- vec_mean("z_sp", length(sp_levels))
z_L1     <- vec_mean("z_L1", length(L1_levels))
z_L2     <- vec_mean("z_L2", length(L2_levels))
z_L3     <- vec_mean("z_L3", length(L3_levels))
z_FT     <- vec_mean("z_FT", length(FT_levels))
if (VARIANT == "v8") {
  z_L1_bgi <- vec_mean("z_L1_bgi", length(L1_levels))
} else {
  z_L1_bgi <- vec_mean("z_L1_site", length(L1_levels))
}
cat(" done\n")

# Knots from meta
if (VARIANT == "v8") {
  climvar <- d$bgi
  knot1 <- meta$bgi_knots[1] %||% NA
  knot2 <- meta$bgi_knots[2] %||% NA
} else {
  # v9: mapdd5 standardized in driver
  climvar <- (d$mapdd5 - meta$mapdd5_mean) / meta$mapdd5_sd
  knot1 <- meta$mapdd5_knots[1]
  knot2 <- meta$mapdd5_knots[2]
}
clim_b1 <- climvar
clim_b2 <- pmax(climvar - knot1, 0)
clim_b3 <- pmax(climvar - knot2, 0)

b_site <- b6 + z_L1_bgi[d$L1_idx] + species_site_slope[d$sp_idx]

eta <- b0 +
  trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
  z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
  z_FT[d$FT_idx] +
  b1 * d$ln_dbh + b2 * d$DBH1 + b3 * d$ln_cr_adj +
  b4 * d$ln_bal_sw_adj + b5 * d$BAL_HW1 +
  b_site * clim_b1 + b9a * clim_b2 + b9b * clim_b3 +
  b7 * (d$BA1 * 0.2296 * d$rd_additive) +
  b8 * (d$BAL_SW1 * d$rd_additive) +
  b11 * d$sdi_complexity +
  b12 * (climvar * d$rd_additive) +
  b13 * (climvar * d$ln_dbh) +
  b14 * (climvar * softwood_per_tree_c) +
  b15 * (climvar * d$ln_cr_adj)

obs_raw  <- d$dg_obs_a
residual <- log(obs_raw) - eta    # NO pmax — we filtered above
weight   <- d$sqrt_years

cat(sprintf("Residual summary: n=%d  mean=%.4f  sd=%.4f  p01=%.3f  p99=%.3f\n",
            sum(is.finite(residual)),
            mean(residual, na.rm = TRUE),
            sd(residual, na.rm = TRUE),
            quantile(residual, 0.01, na.rm = TRUE),
            quantile(residual, 0.99, na.rm = TRUE)))

# Keep modifier-needed columns
keep <- intersect(c("SPCD","sp_idx","EPA_L1_CODE","L1_idx","YEARS",
                     "is_plantation",
                     "had_fire_t1","had_insect_t1","had_disease_t1",
                     "had_wind_t1","had_harvest_t1",
                     "had_cutting_t1","had_site_prep_t1",
                     "years_since_dstrb","years_since_trt",
                     "dstrb_decay_5yr","dstrb_decay_10yr","dstrb_decay_20yr",
                     "trt_decay_5yr","trt_decay_10yr","trt_decay_20yr"),
                   names(d))
out <- d[, ..keep]
out[, eta_base := eta]
out[, obs_raw := obs_raw]
out[, residual := residual]
out[, weight := weight]

# Coerce NAs to 0 on indicators
for (c in grep("^(is_plantation|had_)", names(out), value = TRUE)) {
  v <- out[[c]]; v[is.na(v)] <- 0L; out[[c]] <- as.integer(v)
}
for (c in grep("_decay_", names(out), value = TRUE)) {
  v <- out[[c]]; v[is.na(v)] <- 0.0; out[[c]] <- as.numeric(v)
}

bundle <- list(
  model = paste0("dg_kuehne_", VARIANT),
  family = "log",
  fit_path = FIT_FILE,
  data = out,
  L1_levels = L1_levels,
  n_rows = nrow(out),
  resid_sd = sd(residual, na.rm = TRUE)
)

saveRDS(bundle, OUT_FILE)
cat("\nResidual bundle saved to:", OUT_FILE, "\n")
cat("sigma_resid (empirical):", round(bundle$resid_sd, 3), "\n")
cat("Done.\n")
