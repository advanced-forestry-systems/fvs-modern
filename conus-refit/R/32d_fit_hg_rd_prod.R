##=============================================================================
## 32_fit_hg_speciesfree_v5.R
## Height growth species-free B1 v5: BGI piecewise + FT random effect.
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

STAN_FILE <- get_arg("stan_file", "calibration/stan/hg_organon_speciesfree_v5_bgi.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/hg/v5")
OUT_NAME  <- get_arg("outname",   "hg_sf_v5")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")
HOLDOUT_SPCD_FILE <- get_arg("holdout_spcd_file", NA_character_)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 32_fit_hg_speciesfree_v5.R ==\n")

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

# Slope/aspect: coalesce NA to 0
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
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0 &
  hg_obs_a > -0.5 & hg_obs_a < 5.0
]
cat("After column filters:", nrow(dat), "rows\n")

if (!is.na(HOLDOUT_SPCD_FILE)) {
  holdout_spcd <- as.integer(readLines(HOLDOUT_SPCD_FILE))
  cat(sprintf("Excluding %d holdout species: %s\n",
              length(holdout_spcd),
              paste(holdout_spcd, collapse = ",")))
  pre_n <- nrow(dat)
  dat <- dat[!SPCD %in% holdout_spcd]
  cat(sprintf("  rows before: %d  after: %d  dropped: %d\n",
              pre_n, nrow(dat), pre_n - nrow(dat)))
}

# DF variety split: SPCD 202 -> 2020 (Coastal, EPA_L1=7) or 2021 (Rocky Mt, else).
# Activates only when species_traits_v2.rds is in use (has SPCD 2020 and 2021 rows).
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  pre_df_n <- sum(dat$SPCD == 202L)
  dat[SPCD == 202L & as.character(EPA_L1_CODE) == "7", SPCD := 2020L]
  dat[SPCD == 202L, SPCD := 2021L]
  cat(sprintf("DF variety split: %d SPCD 202 records -> %d coastal (2020) + %d rocky (2021)\n",
              pre_df_n, sum(dat$SPCD == 2020L), sum(dat$SPCD == 2021L)))
}

# Lodgepole pine variety split: SPCD 108 -> 1080 (Shore pine, EPA_L1=7) or 1081 (Rocky Mt, else)
if (any(traits$SPCD == 1080L) && any(traits$SPCD == 1081L)) {
  pre_lp_n <- sum(dat$SPCD == 108L)
  dat[SPCD == 108L & as.character(EPA_L1_CODE) == "7", SPCD := 1080L]
  dat[SPCD == 108L, SPCD := 1081L]
  cat(sprintf("Lodgepole variety split: %d SPCD 108 records -> %d shore (1080) + %d rocky (1081)\n",
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

# Auto-detect traits_v3 (decomposed Potter VCC: CE + S + LAC) vs v2 (composite vuln_score + S)
use_v3_traits <- all(c("climate_exposure", "low_adaptive_cap") %in% names(traits))
if (use_v3_traits) {
  cat("[traits] detected v3 layout: using decomposed Potter components (CE+S+LAC)\n")
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

knots <- quantile(dat$bgi, c(0.33, 0.67), na.rm = TRUE)
cat("BGI knots (33/67%):", round(knots, 3), "\n\n")

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

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,
  FT_idx = dat$FT_idx,
  W = W,

  bgi_knot1 = unname(knots[1]),
  bgi_knot2 = unname(knots[2])
)
## v6_relht-only data: raw HT and species max height
if (grepl("v6_relht|relht", STAN_FILE)) {
  stan_data$ht_obs_m   <- as.numeric(dat$HT1)
  stan_data$max_ht_sp_m <- as.numeric(traits[match(sp_levels, SPCD), max_ht_m])
  if (any(!is.finite(stan_data$max_ht_sp_m))) {
    med <- median(stan_data$max_ht_sp_m, na.rm=TRUE)
    stan_data$max_ht_sp_m[!is.finite(stan_data$max_ht_sp_m)] <- med
    cat("v6_relht: median-filled missing max_ht_m\n")
  }
  cat("v6_relht: added ht_obs_m + max_ht_sp_m\n")
}
## v6_climate_relht: additionally pass climate_si z-scored per observation
if (grepl("climate_relht", STAN_FILE)) {
  cs <- as.numeric(dat$climate_si)
  med <- median(cs, na.rm=TRUE); sdv <- sd(cs, na.rm=TRUE)
  cs[!is.finite(cs)] <- med
  stan_data$climate_si_z <- (cs - med) / pmax(sdv, 0.1)
  cat("v6_climate_relht: added climate_si_z (n=", length(stan_data$climate_si_z), ", range", round(range(stan_data$climate_si_z),2), ")\n")
}

cat("=== Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n\n")

if (grepl("v8rd|v8_rd", STAN_FILE)) { rd <- dat$sdi_additive1 / dat$SDImax_brms; rd[!is.finite(rd)] <- median(rd[is.finite(rd)]); stan_data$rd_additive <- rd; br <- dat$BAL_SW1 + dat$BAL_HW1; br[!is.finite(br)] <- 0; stan_data$bal_raw <- br; message("v8rd: rd+bal_raw added") }
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
meta_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

if (has_flag("compact")) {
  ll <- fit$draws("log_lik", format = "draws_matrix"); loo_res <- loo::loo(ll)
  saveRDS(loo_res, file.path(OUT_DIR, paste0(OUT_NAME, "_loo.rds"))); rm(ll); gc()
  cat("Saved compact LOO (skipped multi-GB save_object).\n")
} else {
  fit$save_object(fit_path)
}
vars <- c("a0", paste0("a", 1:8), "a9a", "a9b", "a10",
          paste0("gamma[", seq_len(ncol(W)), "]"),
          paste0("gamma_site[", seq_len(ncol(W)), "]"),
          "sigma_L1","sigma_L2","sigma_L3","sigma_FT","sigma_L1_bgi","sigma")
if ("sigma_sp" %in% fit$metadata()$stan_variables) vars <- c(vars, "sigma_sp")
vars <- intersect(unique(c(vars,"a_bard","a_blrd")), fit$metadata()$stan_variables)
summary_df <- fit$summary(variables = vars, "mean","median","sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat","ess_bulk","ess_tail")
names(summary_df)[names(summary_df) %in% c("5%","95%")] <- c("q5","q95")
data.table::fwrite(summary_df, summ_path)
saveRDS(list(form = "hg_speciesfree_v5", trait_cols = trait_cols,
             stan_file = STAN_FILE,
             bgi_knots = knots,
             sp_levels = sp_levels, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels,
             FT_levels = FT_levels,
             summary = summary_df, n_obs = stan_data$N_obs,
             wall_min = wall_min), meta_path)

cat("=== a coefs ===\n")
print(summary_df[grepl("^a[0-9]", summary_df$variable), ])
cat("\n=== sigmas ===\n")
print(summary_df[grepl("^sigma", summary_df$variable), ])
cat("\nDone.\n")
