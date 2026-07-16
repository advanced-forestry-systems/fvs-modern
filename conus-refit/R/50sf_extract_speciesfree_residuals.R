##=============================================================================
## 50sf_extract_speciesfree_residuals.R
##
## Modifier-residual extractor for CR / HCB / Mortality species-free base fits.
## Knows the new v5+ architectures (trait_effect + L1/L2/L3 + FT random effect).
##
## CLI:
##   --fit       species-free fit RDS
##   --meta      meta RDS (defaults to fit -> meta)
##   --component cr | hcb | mort
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
COMPONENT <- get_arg("component")  # cr | hcb | mort
PAIRS_FILE  <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")
OUT_FILE  <- get_arg("out")
SUBSAMPLE <- as.integer(get_arg("subsample", "0"))

stopifnot(!is.null(FIT_FILE), file.exists(FIT_FILE))
stopifnot(file.exists(META_FILE))
stopifnot(COMPONENT %in% c("cr","hcb","mort"))
stopifnot(!is.null(OUT_FILE))

cat("== 50sf_extract_speciesfree_residuals.R ==\n")
cat("  fit      :", FIT_FILE, "\n")
cat("  component:", COMPONENT, "\n")
cat("  out      :", OUT_FILE, "\n\n")

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

cat(sprintf("Levels: sp=%d L1=%d L2=%d L3=%d FT=%d\n",
            length(sp_levels), length(L1_levels), length(L2_levels),
            length(L3_levels), length(FT_levels)))

# ln_csi or similar climate
if ("climate_si" %in% names(pairs)) {
  med <- median(pairs$climate_si, na.rm = TRUE)
  pairs[!is.finite(climate_si), climate_si := med]
  pairs[, ln_csi := log(pmax(climate_si, 0.1))]
} else {
  pairs[, ln_csi := 0]
}
pairs[!is.finite(ln_csi), ln_csi := 0]

# Component-specific derived columns + filters
if (COMPONENT == "cr") {
  pairs[, delta_CR_a := (CR2 - CR1) / YEARS]
  filt <- with(pairs,
    is.finite(DBH1) & DBH1 >= 2.54 &
    is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
    is.finite(CR2) & CR2 > 0 & CR2 <= 1.0 &
    is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
    !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
    !is.na(EPA_L2_CODE) & EPA_L2_CODE != "" &
    !is.na(EPA_L3_CODE) & EPA_L3_CODE != "" &
    TREESTATUS1 == 1 & TREESTATUS2 == 1 &
    is.finite(BA1) & BA1 >= 0 &
    is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
    !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0 &
    is.finite(delta_CR_a) & delta_CR_a > -0.5 & delta_CR_a < 0.5
  )
  obs_col <- "delta_CR_a"
  family <- "identity"
  weight_n <- 1
} else if (COMPONENT == "hcb") {
  pairs[, ratio := 1 - CR1]
  pairs[, ln_ht := log(pmax(HT1, 1.5))]
  pairs[, ln_dbh := log(DBH1)]
  pairs[, bal_over_ht := (BAL_SW1 + BAL_HW1) / (HT1 + 1)]
  pairs[, sqrt_ba := sqrt(BA1 * 0.2296)]
  pairs[, ln_cspi_shift := ln_csi]
  filt <- with(pairs,
    is.finite(DBH1) & DBH1 >= 2.54 &
    is.finite(HT1) & HT1 > 1 &
    is.finite(ratio) & ratio > 0.01 & ratio < 0.99 &
    !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
    !is.na(EPA_L2_CODE) & EPA_L2_CODE != "" &
    !is.na(EPA_L3_CODE) & EPA_L3_CODE != "" &
    is.finite(BA1) & BA1 >= 0 &
    is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
    !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
  )
  obs_col <- "ratio"
  family <- "logit"
  weight_n <- 1
} else if (COMPONENT == "mort") {
  pairs[, alive := as.integer(TREESTATUS2 == 1)]
  pairs[, rd_ratio := sdi_additive1 / SDImax_brms]
  pairs[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]
  filt <- with(pairs,
    TREESTATUS1 == 1 & !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1,2) &
    is.finite(DBH1) & DBH1 >= 2.54 &
    is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
    is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
    !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
    !is.na(EPA_L2_CODE) & EPA_L2_CODE != "" &
    !is.na(EPA_L3_CODE) & EPA_L3_CODE != "" &
    is.finite(BA1) & BA1 >= 0 &
    is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
    is.finite(rd_ratio) & rd_ratio >= 0 &
    !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
  )
  obs_col <- "alive"
  family <- "cloglog"
  weight_n <- 1
}

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
  cat("Subsampled to:", nrow(d), "rows\n")
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

cat("Pulling posterior means ...\n")

trait_effect <- vec_mean("trait_effect", length(sp_levels))
z_L1 <- vec_mean("z_L1", length(L1_levels))
z_L2 <- vec_mean("z_L2", length(L2_levels))
z_L3 <- vec_mean("z_L3", length(L3_levels))
z_FT <- vec_mean("z_FT", length(FT_levels))

if (COMPONENT == "cr") {
  b0 <- get_mean("b0"); b1 <- get_mean("b1"); b2 <- get_mean("b2")
  b3 <- get_mean("b3"); b4 <- get_mean("b4"); b5 <- get_mean("b5"); b6 <- get_mean("b6")
  eta <- b0 + trait_effect[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
    b1 * d$DBH1 + b2 * d$DBH1^2 +
    b3 * (d$BA1 * 0.2296) + b4 * (d$BAL_SW1 + d$BAL_HW1) +
    b5 * d$CR1 + b6 * d$ln_csi
  obs_raw <- d$delta_CR_a
  residual <- obs_raw - eta   # identity link
} else if (COMPONENT == "hcb") {
  h0 <- get_mean("h0"); h1 <- get_mean("h1"); h2 <- get_mean("h2")
  h3 <- get_mean("h3"); h4 <- get_mean("h4"); h5 <- get_mean("h5")
  eta <- h0 + trait_effect[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
    h1 * d$ln_ht + h2 * d$ln_dbh + h3 * d$bal_over_ht +
    h4 * d$sqrt_ba + h5 * d$ln_cspi_shift
  obs_raw <- d$ratio
  # Beta likelihood; residual on logit scale
  residual <- qlogis(pmin(pmax(obs_raw, 1e-4), 1 - 1e-4)) - eta
} else if (COMPONENT == "mort") {
  b0 <- get_mean("b0"); b1 <- get_mean("b1"); b2 <- get_mean("b2")
  b3 <- get_mean("b3"); b4 <- get_mean("b4"); b5 <- get_mean("b5"); b6 <- get_mean("b6")
  eta <- b0 + trait_effect[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
    b1 * d$DBH1 + b2 * d$DBH1^2 +
    b3 * d$CR1 + b4 * d$ln_csi +
    b5 * (d$BAL_SW1 + d$BAL_HW1) + b6 * d$sqrt_ba_rd
  # P_surv_annual = exp(-exp(eta))
  # residual on cloglog scale: cloglog(p_mort_annual) - eta
  # P_mort_T = 1 - exp(-exp(eta) * T)
  # cloglog(p_mort_T) - log(T) = eta
  p_mort_obs <- 1 - d$alive
  # Empirical per-tree mortality probability over T years
  # For each row: residual = cloglog(p_mort_T) - eta - log(T)
  # Avoid log(0) when alive = 1
  T <- d$YEARS
  # Clamp p_mort to (eps, 1-eps) so cloglog is finite
  eps <- 1e-3
  p_clamped <- pmin(pmax(p_mort_obs, eps), 1 - eps)
  residual <- log(-log(1 - p_clamped)) - log(T) - eta
  obs_raw <- d$alive
}

weight <- rep(weight_n, nrow(d))

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
  model = paste0(COMPONENT, "_speciesfree"),
  family = family,
  fit_path = FIT_FILE,
  data = out,
  L1_levels = L1_levels,
  n_rows = nrow(out),
  resid_sd = sd(residual, na.rm = TRUE)
)

saveRDS(bundle, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")
cat("sigma_resid (empirical):", round(bundle$resid_sd, 3), "\n")
cat("Done.\n")
