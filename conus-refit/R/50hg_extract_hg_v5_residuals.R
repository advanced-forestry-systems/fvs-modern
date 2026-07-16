##=============================================================================
## 50hg_extract_hg_v5_residuals.R
##
## Modifier-residual extractor for HG_Organon v5 (BGI piecewise + FT RE).
## Same σ-resid-bug fix as 50v9_extract_dg_v9_residuals.R: filter
## hg_obs_a > 0.01 before log() residual.
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
PAIRS_FILE  <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")
OUT_FILE  <- get_arg("out")
SUBSAMPLE <- as.integer(get_arg("subsample", "0"))

stopifnot(!is.null(FIT_FILE), file.exists(FIT_FILE))
stopifnot(file.exists(META_FILE))
stopifnot(!is.null(OUT_FILE))

cat("== 50hg_extract_hg_v5_residuals.R ==\n")
cat("  fit:", FIT_FILE, "\n  out:", OUT_FILE, "\n\n")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

fit   <- readRDS(FIT_FILE)
meta  <- readRDS(META_FILE)
pairs <- as.data.table(readRDS(PAIRS_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))

sp_levels <- meta$sp_levels
L1_levels <- meta$L1_levels
L2_levels <- meta$L2_levels
L3_levels <- meta$L3_levels
FT_levels <- meta$FT_levels

# Recreate derived columns matching base driver
pairs[, hg_obs_a := (HT2 - HT1) / YEARS]
pairs[, sqrt_years := sqrt(YEARS)]
pairs[, ln_dbh := log(DBH1)]
pairs[, ln_ht := log(pmax(HT1, 1.5))]
pairs[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
pairs[, bal_log := log((BAL_SW1 + BAL_HW1) + 5)]
if (!"SLOPE" %in% names(pairs)) pairs[, SLOPE := 0]
if (!"ASPECT" %in% names(pairs)) pairs[, ASPECT := 0]
pairs[!is.finite(SLOPE), SLOPE := 0]
pairs[!is.finite(ASPECT), ASPECT := 0]
pairs[, slope_pct := as.numeric(SLOPE)]
pairs[, cos_aspect := cos(as.numeric(ASPECT) * pi / 180)]

# KEY FIX: hg_obs_a > 0.01 (positive growth only)
filt <- with(pairs,
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(HT1) & HT1 > 1.5 & is.finite(HT2) & HT2 > 1.5 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(bgi) &
  is.finite(BA1) & BA1 >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0 &
  hg_obs_a > 0.01 & hg_obs_a < 5.0)
d <- pairs[filt]
cat("After filters:", nrow(d), "rows\n")

d <- d[SPCD %in% sp_levels]
d[, sp_idx := match(SPCD, sp_levels)]
d[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
d[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
d[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
d[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]
d <- d[!is.na(sp_idx) & !is.na(L1_idx) & !is.na(L2_idx) &
        !is.na(L3_idx) & !is.na(FT_idx)]
cat("After level match:", nrow(d), "rows\n")

if (SUBSAMPLE > 0 && SUBSAMPLE < nrow(d)) {
  set.seed(42)
  d <- d[sample(.N, SUBSAMPLE)]
}

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

a0 <- get_mean("a0"); a1 <- get_mean("a1"); a2 <- get_mean("a2")
a3 <- get_mean("a3"); a4 <- get_mean("a4"); a5 <- get_mean("a5")
a6 <- get_mean("a6"); a7 <- get_mean("a7"); a8 <- get_mean("a8")
a9a <- get_mean("a9a"); a9b <- get_mean("a9b"); a10 <- get_mean("a10")

trait_effect       <- vec_mean("trait_effect", length(sp_levels))
species_site_slope <- vec_mean("species_site_slope", length(sp_levels))
z_L1     <- vec_mean("z_L1", length(L1_levels))
z_L2     <- vec_mean("z_L2", length(L2_levels))
z_L3     <- vec_mean("z_L3", length(L3_levels))
z_FT     <- vec_mean("z_FT", length(FT_levels))
z_L1_bgi <- vec_mean("z_L1_bgi", length(L1_levels))

# Apply BGI piecewise basis (knots from meta)
knot1 <- meta$bgi_knots[1]
knot2 <- meta$bgi_knots[2]
bgi <- d$bgi
bgi_b2 <- pmax(bgi - knot1, 0)
bgi_b3 <- pmax(bgi - knot2, 0)
b_site <- a4 + z_L1_bgi[d$L1_idx] + species_site_slope[d$sp_idx]

eta <- a0 +
  trait_effect[d$sp_idx] +
  z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
  a1 * d$ln_dbh + a2 * d$ln_ht + a3 * d$ln_cr_adj +
  b_site * bgi + a9a * bgi_b2 + a9b * bgi_b3 +
  a5 * d$bal_log + a6 * (d$BA1 * 0.2296) +
  a7 * d$slope_pct + a8 * d$cos_aspect +
  a10 * (bgi * d$bal_log)

obs_raw  <- d$hg_obs_a
residual <- log(obs_raw) - eta
weight   <- d$sqrt_years

cat(sprintf("Residual: n=%d mean=%.4f sd=%.4f p01=%.3f p99=%.3f\n",
            sum(is.finite(residual)),
            mean(residual, na.rm = TRUE),
            sd(residual, na.rm = TRUE),
            quantile(residual, 0.01, na.rm = TRUE),
            quantile(residual, 0.99, na.rm = TRUE)))

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

for (c in grep("^(is_plantation|had_)", names(out), value = TRUE)) {
  v <- out[[c]]; v[is.na(v)] <- 0L; out[[c]] <- as.integer(v)
}
for (c in grep("_decay_", names(out), value = TRUE)) {
  v <- out[[c]]; v[is.na(v)] <- 0.0; out[[c]] <- as.numeric(v)
}

bundle <- list(
  model = "hg_organon_v5",
  family = "log",
  fit_path = FIT_FILE,
  data = out,
  sp_levels = sp_levels,
  L1_levels = L1_levels,
  n_rows = nrow(out),
  resid_sd = sd(residual, na.rm = TRUE)
)
saveRDS(bundle, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\nsigma_resid:", round(bundle$resid_sd, 3), "\n")
