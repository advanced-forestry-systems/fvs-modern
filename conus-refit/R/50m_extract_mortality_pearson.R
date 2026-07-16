##=============================================================================
## 50m_extract_mortality_pearson.R
##
## Mortality residual extractor using Pearson residual approach.
## Bypasses the cloglog-of-Bernoulli infinity issue by using:
##
##   p_mort_T_hat = 1 - exp(-exp(eta_safe) * T)
##   p_surv_T_hat = exp(-exp(eta_safe) * T)
##   residual    = (alive_obs - p_surv_T_hat) / sqrt(p_surv_T_hat * p_mort_T_hat)
##
## Pearson residual is approximately N(0, 1) at convergence for a well-
## specified model, so modifier_common's normal-likelihood machinery
## works without modification. Weight stays 1 (no annualization needed
## since exposure is in eta).
##
## CLI: same as 50sf with --component=mort_pearson
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

cat("== 50m_extract_mortality_pearson.R ==\n")
cat("  fit:", FIT_FILE, "\n")
cat("  out:", OUT_FILE, "\n\n")

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

pairs[, alive := as.integer(TREESTATUS2 == 1)]
pairs[, rd_ratio := sdi_additive1 / SDImax_brms]
pairs[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]
if ("climate_si" %in% names(pairs)) {
  med <- median(pairs$climate_si, na.rm = TRUE)
  pairs[!is.finite(climate_si), climate_si := med]
  pairs[, ln_csi := log(pmax(climate_si, 0.1))]
} else {
  pairs[, ln_csi := 0]
}
pairs[!is.finite(ln_csi), ln_csi := 0]

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

trait_effect <- vec_mean("trait_effect", length(sp_levels))
z_L1 <- vec_mean("z_L1", length(L1_levels))
z_L2 <- vec_mean("z_L2", length(L2_levels))
z_L3 <- vec_mean("z_L3", length(L3_levels))
z_FT <- vec_mean("z_FT", length(FT_levels))

b0 <- get_mean("b0"); b1 <- get_mean("b1"); b2 <- get_mean("b2")
b3 <- get_mean("b3"); b4 <- get_mean("b4"); b5 <- get_mean("b5"); b6 <- get_mean("b6")

eta <- b0 + trait_effect[d$sp_idx] +
  z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
  b1 * d$DBH1 + b2 * d$DBH1^2 +
  b3 * d$CR1 + b4 * d$ln_csi +
  b5 * (d$BAL_SW1 + d$BAL_HW1) + b6 * d$sqrt_ba_rd

# Pearson residual
eta_safe <- pmin(pmax(eta, -20), 5)
p_mort_T <- 1 - exp(-exp(eta_safe) * d$YEARS)
p_mort_T <- pmin(pmax(p_mort_T, 1e-4), 1 - 1e-4)
p_surv_T <- 1 - p_mort_T

residual <- (d$alive - p_surv_T) / sqrt(p_surv_T * p_mort_T)
weight   <- rep(1, nrow(d))

cat(sprintf("Pearson residual: n=%d mean=%.4f sd=%.4f p01=%.3f p99=%.3f\n",
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
out[, p_surv_T := p_surv_T]
out[, p_mort_T := p_mort_T]
out[, obs_raw := d$alive]
out[, residual := residual]
out[, weight := weight]

for (c in grep("^(is_plantation|had_)", names(out), value = TRUE)) {
  v <- out[[c]]; v[is.na(v)] <- 0L; out[[c]] <- as.integer(v)
}
for (c in grep("_decay_", names(out), value = TRUE)) {
  v <- out[[c]]; v[is.na(v)] <- 0.0; out[[c]] <- as.numeric(v)
}

bundle <- list(
  model = "mort_speciesfree_pearson",
  family = "pearson",
  fit_path = FIT_FILE,
  data = out,
  sp_levels = sp_levels,
  L1_levels = L1_levels,
  n_rows = nrow(out),
  resid_sd = sd(residual, na.rm = TRUE)
)

saveRDS(bundle, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")
cat("sigma_resid (empirical Pearson):", round(bundle$resid_sd, 3), "\n")
cat("Done.\n")
