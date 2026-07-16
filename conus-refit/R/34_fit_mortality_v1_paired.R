# Fit species-free mortality v1 (NO cch) on the v2 cch seed-2026 100K data
#
# Purpose: provide a PAIRED baseline for the crown-closure gain. The existing
# production v1 fit (mort_sf_100k_prod) used set.seed(42) and has no cch column
# handling, so its 100K subset differs from the v2 cch / survival fits (seed
# 2026). That makes a paired loo_compare invalid and forces a per-observation
# ELPD contrast. This driver removes that limitation: it reproduces the v2 cch
# data prep EXACTLY (same filter incl. is.finite(cch), same seed 2026, same
# subsample) but fits the v1 Stan model (no cch, no nonlinear CR). The result
# is on the identical 100K rows in identical order as v2 cch, so
#   loo_compare(v2_cch_loo, v1_paired_loo)
# is a valid paired test of the pure gain from crown closure at tree tip plus
# quadratic crown ratio (Greg Johnson's key covariates).
#
# Run from Cardinal (project root /users/PUOM0008/crsfaaron/fvs-conus):
#   module load gcc/12.3.0 R/4.4.0
#   Rscript R/34_fit_mortality_v1_paired.R --subsample=30000   # smoke
#   Rscript R/34_fit_mortality_v1_paired.R                     # full
#
# Author: Aaron Weiskittel
# Date: 2026-05-21

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

# v1 (no cch) Stan model, but data prep mirrors the v2 cch driver
STAN_FILE <- get_arg("stan_file", "calibration/stan/gompit_mortality_speciesfree.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/mort/v1_paired_100k_prod")
OUT_NAME  <- get_arg("outname",   "mort_v1_paired_100k")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- "calibration/traits/species_traits.rds"

dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))

MIN_OBS_SPECIES <- 5000

# ---- Crown closure at tree tip (cch) resolution -------------------------
# Reproduced from the v2 cch driver so the row set after filtering and
# subsampling is IDENTICAL. cch itself is NOT passed to the v1 Stan model;
# it is only computed here to keep the is.finite(cch) filter in play (which
# drops zero rows because non-finite cch is coalesced to 0).
cch_candidates <- c("CCH", "cch", "CCH1", "CCH_TT", "crown_closure_tip",
                    "CCFL_tip", "CC_tip", "cch_init")
cch_col <- intersect(cch_candidates, names(dat))
if (length(cch_col) == 0) {
  stop("Crown closure column not found; v2 cch run used CCH1. Cannot guarantee ",
       "paired data without it. Provide cch or align the filter manually.")
}
cch_col <- cch_col[1]
cat("Using crown closure column (for paired filtering only):", cch_col, "\n")
dat[, cch := as.numeric(get(cch_col))]
dat[!is.finite(cch), cch := 0]

# ---- Derive response and interval (mirror v2 cch driver) ----------------
dat[, alive := as.integer(TREESTATUS2 == 1)]
dat[, T_years := YEARS]

# ---- Standard covariate construction (mirror v2 cch driver) -------------
if ("climate_si" %in% names(dat)) {
  med <- median(dat$climate_si, na.rm = TRUE)
  dat[!is.finite(climate_si), climate_si := med]
  dat[, ln_csi := log(pmax(climate_si, 0.1))]
} else {
  dat[, ln_csi := 0]
}
dat[!is.finite(ln_csi), ln_csi := 0]

dat[, rd_ratio := sdi_additive1 / SDImax_brms]
dat[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]

# Filter is IDENTICAL to the v2 cch driver, including is.finite(cch)
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

# IDENTICAL subsample to v2 cch / survival: set.seed(2026)
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(2026)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
}

# v1 Stan data block: NO cch (everything else mirrors the v2 cch driver)
stan_data <- list(
  N_obs = nrow(dat),
  N_sp  = length(sp_levels),
  N_L1  = length(L1_levels),
  N_L2  = length(L2_levels),
  N_L3  = length(L3_levels),
  N_FT  = length(FT_levels),
  P_trait = ncol(W),
  alive   = as.integer(dat$alive),
  T_years = dat$T_years,
  dbh     = dat$DBH1,
  dbh_sq  = dat$DBH1^2,
  cr_init = dat$CR1,
  ln_csi  = dat$ln_csi,
  bal_metric = (dat$BAL_SW1 + dat$BAL_HW1),
  sqrt_ba_rd = dat$sqrt_ba_rd,
  sp_idx  = match(dat$SPCD, sp_levels),
  L1_idx  = match(as.character(dat$EPA_L1_CODE), L1_levels),
  L2_idx  = match(as.character(dat$EPA_L2_CODE), L2_levels),
  L3_idx  = match(as.character(dat$EPA_L3_CODE), L3_levels),
  FT_idx  = match(as.integer(dat$FORTYPCD_cond1), FT_levels),
  W       = W
)

cat("=== Stan data ready (v1 paired, no cch) ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("Survival rate:", round(mean(stan_data$alive), 4), "\n\n")

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
  adapt_delta = 0.9, max_treedepth = 11, seed = 2026, refresh = 100
)
wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

meta_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

# Compute LOO directly and save the compact loo object (skip the giant matrix).
if (!SMOKE) {
  ll <- fit$draws("log_lik", format = "draws_matrix")
  loo_res <- loo::loo(ll)
  saveRDS(loo_res, file.path(OUT_DIR, paste0(OUT_NAME, "_loo.rds")))
  cat("\nLOO ELPD (mortality v1 paired):",
      sprintf("%.1f (SE %.1f)\n",
              loo_res$estimates["elpd_loo", "Estimate"],
              loo_res$estimates["elpd_loo", "SE"]))
  rm(ll); gc()
}

# v1 parameter set: NO b3b / b7 / b7b
vars <- c("b0", paste0("b", 1:6),
          paste0("gamma[", 1:ncol(W), "]"),
          "sigma_L1", "sigma_L2", "sigma_L3", "sigma_FT")
summary_df <- fit$summary(variables = vars, "mean", "median", "sd",
                          ~quantile(.x, probs = c(0.05, 0.95)), "rhat", "ess_bulk", "ess_tail")
names(summary_df)[names(summary_df) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summary_df, summ_path)

saveRDS(list(form = "mort_v1_paired", trait_cols = trait_cols,
             sp_levels = sp_levels, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels, FT_levels = FT_levels,
             summary = summary_df, n_obs = stan_data$N_obs,
             seed = 2026, wall_min = wall_min),
        meta_path)

cat("\n=== v1 paired fit complete ===\n")
cat("Wall minutes:", round(wall_min, 1), "\n")
cat("This fit is on the SAME seed-2026 100K rows as v2 cch and survival,\n")
cat("so loo_compare(v2_cch_loo, v1_paired_loo) is a valid paired test of\n")
cat("the crown-closure + nonlinear-CR gain.\n")
