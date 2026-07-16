##=============================================================================
## 32_fit_dg_kuehne_v8c.R
##
## DG Kuehne v8c driver: v8 base + additive CSPI v4 productivity term.
## Mirrors the CSPI ingestion pattern from 36_fit_htdbh_speciesfree.R
##
## CLI (in addition to v8):
##   --cspi_file=PATH  CSV with PLT_CN_cond1, cspi_v4 columns
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-25
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

STAN_FILE <- get_arg("stan_file", "calibration/stan/dg_kuehne2022_v8c_bgi_cspi.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/dg_kue/v8c")
OUT_NAME  <- get_arg("outname",   "dg_kuehne_v8c")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
BGI_Q1    <- as.numeric(get_arg("bgi_q1", "0.33"))
BGI_Q2    <- as.numeric(get_arg("bgi_q2", "0.67"))
SMOKE     <- has_flag("smoke")
HOLDOUT_SPCD_FILE <- get_arg("holdout_spcd_file", NA_character_)
CSPI_FILE <- get_arg("cspi_file", NA_character_)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 32_fit_dg_kuehne_v8c.R ==\n")
cat("Stan:", STAN_FILE, "\n")
cat("CSPI file:", if (!is.na(CSPI_FILE)) CSPI_FILE else "(NONE - aborting)", "\n\n")

if (is.na(CSPI_FILE) || !file.exists(CSPI_FILE)) {
  stop("v8c requires --cspi_file=... pointing to cspi_v4_at_calib_plots.csv")
}

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done. Rows:", nrow(dat), "\n\n")

## ---- CSPI v4 join (mirrors 36_fit_htdbh_speciesfree.R) -----------------
suppressMessages(library(bit64))
cspi_lk <- data.table::fread(CSPI_FILE)
dat[, .pid := as.character(PLT_CN_cond1)]
cspi_lk[, .pid := as.character(PLT_CN_cond1)]
dat[cspi_lk, cspi_ext := i.cspi_v4, on = ".pid"]
n_cspi_found <- sum(is.finite(dat$cspi_ext))
.med_cspi <- median(dat$cspi_ext, na.rm = TRUE)
dat[, cspi_v4 := data.table::fifelse(is.finite(cspi_ext), cspi_ext, .med_cspi)]
dat[, ln_cspi_shift := log(pmax(cspi_v4, 0.1))]
dat[!is.finite(ln_cspi_shift), ln_cspi_shift := log(.med_cspi)]
cat(sprintf("CSPI v4 join: matched %d of %d rows (%.1f%%); median-filled rest at %.2f m\n",
            n_cspi_found, nrow(dat), 100 * n_cspi_found / nrow(dat), .med_cspi))

## 1. Kuehne data prep -----------------------------------------------------
MIN_OBS_SPECIES <- 5000

dat[, dg_obs_a       := (DBH2 - DBH1) / YEARS]
dat[, sqrt_years     := sqrt(YEARS)]
dat[, ln_dbh         := log(DBH1)]
dat[, ln_cr_adj      := log((CR1 + 0.2) / 1.2)]
dat[, ln_bal_sw_adj  := log(BAL_SW1 + 0.01)]

dat[, rd_additive    := sdi_additive1 / SDImax_brms]
dat[, sdi_complexity := sdi_additive1 / pmax(SDI1, 1.0)]

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
  is.finite(ln_cspi_shift) &
  dg_obs_a > -0.5 & dg_obs_a < 5.0
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

# Variety splits (same as v8)
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  pre_df_n <- sum(dat$SPCD == 202L)
  dat[SPCD == 202L & as.character(EPA_L1_CODE) == "7", SPCD := 2020L]
  dat[SPCD == 202L, SPCD := 2021L]
  cat(sprintf("DF variety split: %d SPCD 202 records -> %d coastal (2020) + %d rocky (2021)\n",
              pre_df_n, sum(dat$SPCD == 2020L), sum(dat$SPCD == 2021L)))
}
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

## 2. Indices --------------------------------------------------------------
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

## 3. Trait matrix ---------------------------------------------------------
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

softwood_by_sp <- traits_sub$softwood
softwood_by_sp[is.na(softwood_by_sp)] <- 0

for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}

softwood_per_tree <- softwood_by_sp[dat$sp_idx]
softwood_per_tree <- softwood_per_tree - mean(softwood_per_tree)

## 4. Subsample ------------------------------------------------------------
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  softwood_per_tree <- softwood_per_tree[idx]
  cat("Subsampled to:", nrow(dat), "rows\n\n")
}

## 5. BGI knots -----------------------------------------------------------
bgi_knots <- quantile(dat$bgi, c(BGI_Q1, BGI_Q2), na.rm = TRUE)
cat("BGI knots placed at q", BGI_Q1, "/", BGI_Q2, " = ",
    round(bgi_knots[1], 3), " / ", round(bgi_knots[2], 3), "\n", sep="")
cat("ln_cspi_shift summary: range",
    round(min(dat$ln_cspi_shift), 3), "to", round(max(dat$ln_cspi_shift), 3),
    "median", round(median(dat$ln_cspi_shift), 3), "\n\n")

## 6. Stan data -----------------------------------------------------------
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
  softwood = softwood_per_tree,

  ln_cspi_shift = dat$ln_cspi_shift,

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,
  FT_idx = dat$FT_idx,
  W = W,

  bgi_knot1 = unname(bgi_knots[1]),
  bgi_knot2 = unname(bgi_knots[2])
)

cat("=== Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n\n")

## 7. Compile + sample ---------------------------------------------------
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

## 8. Save ----------------------------------------------------------------
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

vars <- c("b0", paste0("b", 1:8), "b9a", "b9b",
          paste0("b", 11:15), "a_cspi",
          paste0("gamma[", seq_len(ncol(W)), "]"),
          paste0("gamma_site[", seq_len(ncol(W)), "]"),
          "sigma_sp", "sigma_L1", "sigma_L2", "sigma_L3",
          "sigma_FT", "sigma_L1_bgi", "sigma")
summary_df <- fit$summary(variables = vars, "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail")
names(summary_df)[names(summary_df) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summary_df, summ_path)
saveRDS(list(form = "kuehne_v8c", trait_cols = trait_cols,
             stan_file = STAN_FILE,
             sp_levels = sp_levels, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels,
             FT_levels = FT_levels,
             bgi_knots = c(stan_data$bgi_knot1, stan_data$bgi_knot2),
             cspi_file = CSPI_FILE,
             summary = summary_df, n_obs = stan_data$N_obs,
             wall_min = wall_min), meta_path)

cat("=== KEY COEFFICIENT (v8c) ===\n")
print(summary_df[summary_df$variable == "a_cspi", ])
cat("\n=== BGI block ===\n")
print(summary_df[summary_df$variable %in%
                  c("b6","b9a","b9b","b12","b13","b14","b15"), ])
cat("\n=== sigmas ===\n")
print(summary_df[grepl("^sigma", summary_df$variable), ])
cat("\nDone.\n")
