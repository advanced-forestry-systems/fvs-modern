##=============================================================================
## 32_fit_hg_speciesfree_v8.R
## HG v8 driver: v7 + CCH crown closure at tree height
##=============================================================================

library(data.table)
library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/hg_unified_v8_bgi.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/hg/v8")
OUT_NAME  <- get_arg("outname",   "hg_sf_v8")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 32_fit_hg_speciesfree_v8.R ==\n")
cat("Stan:", STAN_FILE, "\n\n")

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done. Rows:", nrow(dat), "\n\n")

MIN_OBS_SPECIES <- 5000

dat[, hg_obs_a := (HT2 - HT1) / YEARS]
dat[, sqrt_years := sqrt(YEARS)]
dat[, ln_dbh := log(DBH1)]
dat[, ln_ht := log(pmax(HT1, 1.5))]
dat[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
dat[, bal_log := log((BAL_SW1 + BAL_HW1) + 5)]

if (!"SLOPE" %in% names(dat)) dat[, SLOPE := 0]
if (!"ASPECT" %in% names(dat)) dat[, ASPECT := 0]
dat[!is.finite(SLOPE), SLOPE := 0]
dat[!is.finite(ASPECT), ASPECT := 0]
dat[, slope_pct := as.numeric(SLOPE)]
dat[, cos_aspect := cos(as.numeric(ASPECT) * pi / 180)]
dat[!is.finite(slope_pct), slope_pct := 0]
dat[!is.finite(cos_aspect), cos_aspect := 0]

dat <- dat[
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(HT1) & HT1 > 1.5 & is.finite(HT2) & HT2 > 1.5 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 &
  is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(bgi) &
  is.finite(BA1) & BA1 >= 0 &
  is.finite(CCH1) & CCH1 >= 0 &                # v8: NEW filter
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0 &
  hg_obs_a > -0.5 & hg_obs_a < 5.0
]
cat("After column filters:", nrow(dat), "rows\n")

# Variety splits
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  pre_df_n <- sum(dat$SPCD == 202L)
  dat[SPCD == 202L & as.character(EPA_L1_CODE) == "7", SPCD := 2020L]
  dat[SPCD == 202L, SPCD := 2021L]
  cat(sprintf("DF variety split: %d -> %d coastal + %d rocky\n",
              pre_df_n, sum(dat$SPCD == 2020L), sum(dat$SPCD == 2021L)))
}
if (any(traits$SPCD == 1080L) && any(traits$SPCD == 1081L)) {
  pre_lp_n <- sum(dat$SPCD == 108L)
  dat[SPCD == 108L & as.character(EPA_L1_CODE) == "7", SPCD := 1080L]
  dat[SPCD == 108L, SPCD := 1081L]
  cat(sprintf("LP variety split: %d -> %d shore + %d rocky\n",
              pre_lp_n, sum(dat$SPCD == 1080L), sum(dat$SPCD == 1081L)))
}

sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("After species filter:", nrow(dat), "rows;", nrow(sp_counts), "species\n\n")

sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))
FT_levels <- sort(unique(as.integer(dat$FORTYPCD_cond1)))

dat[, sp_idx := match(SPCD, sp_levels)]
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
dat[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]

use_v3_traits <- all(c("climate_exposure", "low_adaptive_cap") %in% names(traits))
if (use_v3_traits) {
  cat("[traits] v3 layout\n")
  trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                  "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                  "climate_exposure", "sensitivity", "low_adaptive_cap")
} else {
  trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                  "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                  "vulnerability_score", "sensitivity")
}
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  cat("Subsampled to:", nrow(dat), "rows\n\n")
}

bgi_knots <- quantile(dat$bgi, c(0.33, 0.67), na.rm = TRUE)
cat("BGI knots:", round(bgi_knots, 3), "\n")
cat("CCH1: range", round(min(dat$CCH1), 3), "to", round(max(dat$CCH1), 3),
    "median", round(median(dat$CCH1), 3), "\n\n")

stan_data <- list(
  N_obs = nrow(dat),
  N_sp = length(sp_levels),
  N_L1 = length(L1_levels),
  N_L2 = length(L2_levels),
  N_L3 = length(L3_levels),
  N_FT = length(FT_levels),
  P_trait = ncol(W),

  hg_obs_a = dat$hg_obs_a,
  sqrt_years = dat$sqrt_years,
  ln_dbh = dat$ln_dbh,
  ln_ht = dat$ln_ht,
  ln_cr_adj = dat$ln_cr_adj,
  bal_log = dat$bal_log,
  bgi = dat$bgi,
  ba_metric = dat$BA1 * 0.2296,
  slope_pct = dat$slope_pct,
  cos_aspect = dat$cos_aspect,
  cch = dat$CCH1,                            # v8: NEW

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,
  FT_idx = dat$FT_idx,
  W = W,
  bgi_knot1 = unname(bgi_knots[1]),
  bgi_knot2 = unname(bgi_knots[2])
)

cat("=== Stan data ready ===\nN_obs =", stan_data$N_obs, "\n\n")
mod <- cmdstan_model(STAN_FILE)

if (SMOKE) {
  iter_warmup <- 50; iter_sampling <- 50; chains <- 2
} else {
  iter_warmup <- 1000; iter_sampling <- 1000; chains <- 4
}

t_start <- Sys.time()
fit <- mod$sample(
  data = stan_data, chains = chains, parallel_chains = chains,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  seed = 42, adapt_delta = 0.9, max_treedepth = 10, refresh = 100
)
wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("\nWall:", round(wall_min, 1), "min\n\n")

fit_path  <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

if (has_flag("compact")) {
  ll <- fit$draws("log_lik", format = "draws_matrix"); loo_res <- loo::loo(ll)
  saveRDS(loo_res, file.path(OUT_DIR, paste0(OUT_NAME, "_loo.rds"))); rm(ll); gc()
  cat("Saved compact LOO.\n")
} else {
  fit$save_object(fit_path)
}

vars <- c("a0", "a1", "a2", "a2_quad", "a3", "a4", "a5", "a6", "a7", "a8",
          "a9a", "a9b", "a10", "a11", "a11b",
          paste0("gamma[", seq_len(ncol(W)), "]"),
          paste0("gamma_site[", seq_len(ncol(W)), "]"),
          "sigma_L1", "sigma_L2", "sigma_L3", "sigma_FT",
          "sigma_sp", "sigma_L1_bgi", "sigma")
summary_df <- fit$summary(variables = vars, "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail")
names(summary_df)[names(summary_df) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summary_df, summ_path)
cat("=== KEY: CCH effect ===\n")
print(summary_df[summary_df$variable %in% c("a11", "a11b"), ])
cat("\nDone.\n")
