##=============================================================================
## 32_fit_dg_kuehne_v7.R
##
## DG Kuehne v7 driver: species-aware (with traits) + forest type + ecoregion
## random effects, BGI as climate-derived site variable, continuous RD
## interactions, sdi_complexity main effect.
##
## Aaron's biometric directives (2026-05-15):
##   - Drop SICOND (circular)
##   - Drop ELEV (not trusted)
##   - Use BGI directly as site variable (non-circular)
##   - Species, forest type, ecoregion as random effects
##   - Plantation back to modifier (not in this base model)
##   - Continuous RD modifiers (BA*RD, BAL*RD, BGI*RD)
##
## CLI:
##   --stan_file=PATH (default v7 stan)
##   --outdir=PATH
##   --outname=NAME
##   --subsample=N
##   --smoke
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-15
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

STAN_FILE <- get_arg("stan_file", "calibration/stan/dg_kuehne2022_v7_climate_hierarchy.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/dg_kue/v7")
OUT_NAME  <- get_arg("outname",   "dg_kuehne_v7")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 32_fit_dg_kuehne_v7.R ==\n")
cat("Stan:", STAN_FILE, "\n\n")

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- "calibration/traits/species_traits.rds"

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done. Rows:", nrow(dat), "\n\n")

## 1. Kuehne-specific data prep --------------------------------------------
MIN_OBS_SPECIES <- 5000

dat[, dg_obs_a       := (DBH2 - DBH1) / YEARS]
dat[, sqrt_years     := sqrt(YEARS)]
dat[, ln_dbh         := log(DBH1)]
dat[, ln_cr_adj      := log((CR1 + 0.2) / 1.2)]
dat[, ln_bal_sw_adj  := log(BAL_SW1 + 0.01)]

# Derived per Aaron's biometric guidance
dat[, rd_additive    := sdi_additive1 / SDImax_brms]
dat[, sdi_complexity := sdi_additive1 / pmax(SDI1, 1.0)]

# Climate variable: BGI (non-circular, climate-derived)
# Used directly (not log-transformed) due to negative values (76% are negative)

dat <- dat[
  is.finite(DBH1) & DBH1 >= 2.54 & is.finite(DBH2) &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 &
  is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(bgi) &
  is.finite(rd_additive) & rd_additive > 0 & rd_additive < 3.0 &
  is.finite(sdi_complexity) & sdi_complexity > 0 & sdi_complexity < 10 &
  is.finite(BA1) & BA1 >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0 &
  dg_obs_a > -0.5 & dg_obs_a < 5.0
]
cat("After column filters:", nrow(dat), "rows\n")

sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("After species filter:", nrow(dat), "rows;", nrow(sp_counts), "species\n\n")

## 2. Build indices (species, ecoregion, forest type) ----------------------
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

cat("N_sp/N_L1/N_L2/N_L3/N_FT =",
    length(sp_levels), "/", length(L1_levels), "/",
    length(L2_levels), "/", length(L3_levels), "/",
    length(FT_levels), "\n\n")

## 3. Trait matrix W -------------------------------------------------------
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

## 4. Subsample if requested -----------------------------------------------
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  cat("Subsampled to:", nrow(dat), "rows\n\n")
}

## 5. Stan data ------------------------------------------------------------
stan_data <- list(
  N_obs = nrow(dat),
  N_sp = length(sp_levels),
  N_L1 = length(L1_levels),
  N_L2 = length(L2_levels),
  N_L3 = length(L3_levels),
  N_FT = length(FT_levels),
  P_trait = ncol(W),

  dg_obs_a = dat$dg_obs_a,
  sqrt_years = dat$sqrt_years,
  ln_dbh = dat$ln_dbh,
  dbh = dat$DBH1,
  ln_cr_adj = dat$ln_cr_adj,
  ln_bal_sw_adj = dat$ln_bal_sw_adj,
  bal_hw = dat$BAL_HW1,

  bgi = dat$bgi,
  ba_metric = dat$BA1 * 0.2296,
  bal_sw_metric = dat$BAL_SW1,
  rd_additive = dat$rd_additive,
  sdi_complexity = dat$sdi_complexity,

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,
  FT_idx = dat$FT_idx,
  W = W
)

cat("=== Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("bgi range:", round(quantile(stan_data$bgi, c(0,0.5,1)), 2), "\n")
cat("rd_additive range:", round(quantile(stan_data$rd_additive, c(0,0.5,1)), 3), "\n\n")

## 6. Compile + sample -----------------------------------------------------
cat("Compiling Stan model:", STAN_FILE, "\n"); flush.console()
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

## 7. Save -----------------------------------------------------------------
fit_path  <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
meta_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

fit$save_object(fit_path)
vars <- c(paste0("b", c(0:9, 11:12)),
          paste0("gamma[", seq_len(ncol(W)), "]"),
          paste0("gamma_site[", seq_len(ncol(W)), "]"),
          "sigma_sp", "sigma_L1", "sigma_L2", "sigma_L3",
          "sigma_FT", "sigma_L1_bgi", "sigma")
summary_df <- fit$summary(variables = vars, "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail")
names(summary_df)[names(summary_df) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summary_df, summ_path)
saveRDS(list(form = "kuehne_v7", trait_cols = trait_cols,
             stan_file = STAN_FILE,
             sp_levels = sp_levels, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels,
             FT_levels = FT_levels,
             summary = summary_df, n_obs = stan_data$N_obs,
             wall_min = wall_min), meta_path)

cat("=== b coefficients ===\n")
print(summary_df[grepl("^b[0-9]", summary_df$variable), ])
cat("\n=== sigmas + sigma_sp + sigma_FT ===\n")
print(summary_df[grepl("^sigma", summary_df$variable), ])
cat("\nDone.\n")
