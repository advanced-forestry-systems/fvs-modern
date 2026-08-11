# Fit species-free SURVIVAL-framed model (Greg's framing, our architecture)
#
# Models the probability of survival directly (higher eta = higher survival)
# in our species-free hierarchical structure, with Greg's predictors
# (nonlinear crown ratio + crown closure at tree tip). Mathematically
# equivalent in fit to the mortality-framed model but coefficients read on
# the survival scale.
#
# Run from Cardinal (project root /users/PUOM0008/crsfaaron/fvs-conus):
#   module load gdal/3.7.3 gcc/12.3.0 geos/3.12.0 proj/9.2.1 R/4.4.0
#   Rscript R/34_fit_survival_speciesfree.R --subsample=30000   # smoke
#   Rscript R/34_fit_survival_speciesfree.R                     # full
#
# Author: Aaron Weiskittel
# Date: 2026-05-20

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/gompit_survival_speciesfree.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/mort/survival_speciesfree")
OUT_NAME  <- get_arg("outname",   "survival_speciesfree")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")
USE_CCH   <- !has_flag("no_cch")     # allow disabling cch if not in data yet

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")

dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))

MIN_OBS_SPECIES <- 5000

# ---- Crown closure at tree tip (cch) resolution -------------------------
cch_candidates <- c("CCH", "cch", "CCH1", "CCH_TT", "crown_closure_tip",
                    "CCFL_tip", "CC_tip", "cch_init")
cch_col <- intersect(cch_candidates, names(dat))
if (length(cch_col) == 0) {
  if (USE_CCH) {
    message("WARNING: crown closure at tree tip not found in data. ",
            "Falling back to cch = 0 (cch terms will be inert). ",
            "Provide cch in the data or pass --no_cch to silence this. ",
            "See GREG_MORTALITY_RECONCILIATION for the data dependency.")
  }
  dat[, cch := 0]
} else {
  cch_col <- cch_col[1]
  message("Using crown closure column: ", cch_col)
  dat[, cch := as.numeric(get(cch_col))]
  dat[!is.finite(cch), cch := 0]
}

# ---- Derive response and interval (mirror v1 driver) --------------------
dat[, alive := as.integer(TREESTATUS2 == 1)]       # survived if status2 == 1
dat[, T_years := YEARS]                             # measurement interval

# ---- Standard covariates (mirror v1 driver) -----------------------------
if ("climate_si" %in% names(dat)) {
  med <- median(dat$climate_si, na.rm = TRUE)
  dat[!is.finite(climate_si), climate_si := med]
  dat[, ln_csi := log(pmax(climate_si, 0.1))]
} else {
  dat[, ln_csi := 0]
}
dat[!is.finite(ln_csi), ln_csi := 0]
dat[, rd_ratio := sdi_additive1 / SDImax_brms]      # relative density (mirrors v1 driver)
dat[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]

dat <- dat[
  TREESTATUS1 == 1 &
  !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1, 2) &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(cch) &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  is.finite(BA1) & BA1 >= 0 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(rd_ratio) & rd_ratio >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
]

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

sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))
FT_levels <- sort(unique(as.integer(dat$FORTYPCD_cond1)))

trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(2026)
  dat <- dat[sort(sample.int(nrow(dat), SUBSAMPLE))]
}

stan_data <- list(
  N_obs = nrow(dat), N_sp = length(sp_levels),
  N_L1 = length(L1_levels), N_L2 = length(L2_levels),
  N_L3 = length(L3_levels), N_FT = length(FT_levels),
  P_trait = ncol(W),
  alive = as.integer(dat$alive), T_years = dat$T_years,
  dbh = dat$DBH1, dbh_sq = dat$DBH1^2,
  cr_init = dat$CR1, ln_csi = dat$ln_csi,
  bal_metric = (dat$BAL_SW1 + dat$BAL_HW1),
  sqrt_ba_rd = dat$sqrt_ba_rd, cch = dat$cch,
  sp_idx = match(dat$SPCD, sp_levels),
  L1_idx = match(as.character(dat$EPA_L1_CODE), L1_levels),
  L2_idx = match(as.character(dat$EPA_L2_CODE), L2_levels),
  L3_idx = match(as.character(dat$EPA_L3_CODE), L3_levels),
  FT_idx = match(as.integer(dat$FORTYPCD_cond1), FT_levels),
  W = W
)

cat("=== Survival species-free Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("Survival rate:", round(mean(stan_data$alive), 4), "\n\n")

mod <- cmdstan_model(STAN_FILE)
if (SMOKE) { iw <- 50; is_ <- 50; ch <- 2 } else { iw <- 1000; is_ <- 1000; ch <- 4 }

t0 <- Sys.time()
fit <- mod$sample(data = stan_data, chains = ch, parallel_chains = ch,
                  iter_warmup = iw, iter_sampling = is_,
                  adapt_delta = 0.9, max_treedepth = 11, seed = 2026, refresh = 100)
wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# Parameter summary first (fast: restricted to the model coefficients).
vars <- c("b0", paste0("b", 1:6), "b3b", "b7", "b7b",
          paste0("gamma[", 1:ncol(W), "]"),
          "sigma_L1", "sigma_L2", "sigma_L3", "sigma_FT")
summ <- fit$summary(variables = vars, "mean", "median", "sd",
                    ~quantile(.x, probs = c(0.05, 0.95)), "rhat", "ess_bulk", "ess_tail")
names(summ)[names(summ) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv")))

# Compute LOO directly and save the compact loo object (a few MB), instead of
# serializing the full multi-GB log_lik matrix. Skip for smoke runs.
if (!SMOKE) {
  ll <- fit$draws("log_lik", format = "draws_matrix")
  loo_res <- loo::loo(ll)
  saveRDS(loo_res, file.path(OUT_DIR, paste0(OUT_NAME, "_loo.rds")))
  cat("\nLOO ELPD (survival):",
      sprintf("%.1f (SE %.1f)\n",
              loo_res$estimates["elpd_loo", "Estimate"],
              loo_res$estimates["elpd_loo", "SE"]))
  rm(ll); gc()
  # Small parameter draws for downstream coefficient work.
  pdraws <- fit$draws(variables = vars, format = "draws_matrix")
  saveRDS(pdraws, file.path(OUT_DIR, paste0(OUT_NAME, "_param_draws.rds")))

  # Save the ecoregion / forest-type random-effect summary tables, keyed by the
  # real codes, in the same format 61b emits. These are small (the z_* vectors,
  # not the 100K log_lik), so this avoids the save_object hang on the full fit
  # while still giving the species-free bundle everything it needs for the
  # mortality (survival) component to be extracted, benchmarked, and integrated.
  save_re <- function(zname, levels_vec, tag) {
    s <- tryCatch(fit$summary(variables = zname, "mean", "sd"),
                  error = function(e) NULL)
    if (is.null(s) || nrow(s) != length(levels_vec)) return(invisible(NULL))
    data.table::fwrite(
      data.frame(level = levels_vec, mean = s$mean, sd = s$sd),
      file.path(OUT_DIR, paste0(OUT_NAME, "_re_", tag, ".csv")))
  }
  save_re("z_L1", L1_levels, "L1"); save_re("z_L2", L2_levels, "L2")
  save_re("z_L3", L3_levels, "L3"); save_re("z_FT", FT_levels, "FT")
  cat("Saved RE summary tables (L1/L2/L3/FT) for species-free bundle assembly.\n")
}

saveRDS(list(form = "survival_speciesfree", trait_cols = trait_cols,
             cch_used = (length(cch_col) > 0),
             summary = summ, n_obs = stan_data$N_obs, wall_min = wall_min),
        file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds")))

cat("\n=== survival fit complete (wall", round(wall_min, 1), "min) ===\n")
cat("b3 (CR survival effect, expect positive):",
    round(summ[summ$variable == "b3", "mean"][[1]], 4), "\n")
cat("b5 (BAL survival effect, expect negative):",
    round(summ[summ$variable == "b5", "mean"][[1]], 4), "\n")
cat("\nNext: run compare_survival_vs_mortality.R to confirm LOO equivalence ",
    "with the mortality-framed fit.\n")
