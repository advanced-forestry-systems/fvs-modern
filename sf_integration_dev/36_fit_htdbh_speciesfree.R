##=============================================================================
## 36_fit_htdbh_speciesfree.R
##
## Species-free B1 height-diameter (Wykoff log-normal) fit driver.
## Architecture: trait_effect (W*gamma) + L1/L2/L3 ecoregion REs + forest type
## RE, NO species random intercept. Mirrors the data convention of the other
## B1 species-free components (CR, HCB, mortality): v2 data, BAL_SW1+BAL_HW1,
## sdi_additive1/SDImax_brms relative density, climate_si, and the
## Douglas-fir / lodgepole variety split when v2 traits are in use.
##
## CLI:
##   --stan_file=PATH  --traits=PATH  --outdir=PATH  --outname=NAME
##   --subsample=N  --smoke
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-21
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

STAN_FILE <- get_arg("stan_file", "calibration/stan/ht_dbh_wykoff_speciesfree.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/ht_dbh/speciesfree")
OUT_NAME  <- get_arg("outname",   "htdbh_speciesfree")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 36_fit_htdbh_speciesfree.R ==\n")
cat("Stan:", STAN_FILE, "\n\n")

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))

# --- optional external CSPI productivity swap (overwrites climate_si) ---
.CSPI_FILE <- get_arg("cspi_file", NA_character_)
if (!is.na(.CSPI_FILE) && file.exists(.CSPI_FILE)) {
  suppressMessages(library(bit64))
  .lk <- data.table::fread(.CSPI_FILE)
  dat[, .pid := as.character(PLT_CN_cond1)]
  .lk[, .pid := as.character(PLT_CN_cond1)]
  dat[.lk, cspi_ext := i.cspi_v4, on = ".pid"]
  .med <- median(dat$cspi_ext, na.rm = TRUE)
  dat[, climate_si := data.table::fifelse(is.finite(cspi_ext), cspi_ext, .med)]
  cat(sprintf("CSPI swap: external cspi_v4 for %d of %d rows (median fill rest)\n",
              sum(is.finite(dat$cspi_ext)), nrow(dat)))
}

traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done. Rows:", nrow(dat), "\n\n")

MIN_OBS_SPECIES <- as.integer(get_arg("min_sp", "5000"))  # lower for full species coverage under unified z_sp

# ln(climate site index), shifted, coalesce NA to median
if ("climate_si" %in% names(dat)) {
  med <- median(dat$climate_si, na.rm = TRUE)
  dat[!is.finite(climate_si), climate_si := med]
  dat[, ln_cspi_shift := log(pmax(climate_si, 0.1))]
} else if ("CSPI" %in% names(dat)) {
  med <- median(dat$CSPI, na.rm = TRUE)
  dat[!is.finite(CSPI), CSPI := med]
  dat[, ln_cspi_shift := log(pmax(CSPI, 0.1))]
} else {
  dat[, ln_cspi_shift := 0]
}
dat[!is.finite(ln_cspi_shift), ln_cspi_shift := 0]

# Relative density and density interactions (mirror mortality/CR/HCB convention)
dat[, rd_ratio := sdi_additive1 / SDImax_brms]
dat[, bal_metric := (BAL_SW1 + BAL_HW1)]
dat[, ba_metric  := BA1 * 0.2296]
dat[, sqrt_ba    := sqrt(pmax(ba_metric, 0))]
dat[, ba_x_rd    := ba_metric  * rd_ratio]
dat[, bal_x_rd   := bal_metric * rd_ratio]

dat <- dat[
  TREESTATUS1 == 1 &
  is.finite(HT1) & HT1 >= 1.37 & HT1 < 85 &
  is.finite(DBH1) & DBH1 >= 2.54 & DBH1 < 250 &
  is.finite(BA1) & BA1 > 0 &
  is.finite(bal_metric) & bal_metric >= 0 &
  is.finite(rd_ratio) & rd_ratio > 0 & rd_ratio < 2 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
]
cat("After column filters:", nrow(dat), "rows\n")

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

cat("N_sp/N_L1/N_L2/N_L3/N_FT =",
    length(sp_levels), "/", length(L1_levels), "/",
    length(L2_levels), "/", length(L3_levels), "/",
    length(FT_levels), "\n\n")

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
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  cat("Subsampled to:", nrow(dat), "rows\n\n")
}

stan_data <- list(
  N_obs = nrow(dat),
  N_sp = length(sp_levels),
  N_L1 = length(L1_levels),
  N_L2 = length(L2_levels),
  N_L3 = length(L3_levels),
  N_FT = length(FT_levels),
  P_trait = ncol(W),

  ht_obs = dat$HT1,
  dbh = dat$DBH1,
  bal = dat$bal_metric,
  sqrt_ba = dat$sqrt_ba,
  ln_cspi_shift = dat$ln_cspi_shift,
  ba_x_rd = dat$ba_x_rd,
  bal_x_rd = dat$bal_x_rd,

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,
  FT_idx = dat$FT_idx,
  W = W
)

cat("=== Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("HT1 range:", round(quantile(stan_data$ht_obs, c(0,0.5,1)), 2), "\n\n")

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
vars <- c("b0", "b1", "a_bal", "a_ba", "a_cspi", "a_bard", "a_blrd",
          paste0("gamma[", seq_len(ncol(W)), "]"),
          "sigma_L1","sigma_L2","sigma_L3","sigma_FT","sigma")
if ("sigma_sp" %in% fit$metadata()$stan_variables) vars <- c(vars, "sigma_sp")
summary_df <- fit$summary(variables = vars, "mean","median","sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat","ess_bulk","ess_tail")
names(summary_df)[names(summary_df) %in% c("5%","95%")] <- c("q5","q95")
data.table::fwrite(summary_df, summ_path)
saveRDS(list(form = "htdbh_speciesfree", trait_cols = trait_cols,
             stan_file = STAN_FILE,
             sp_levels = sp_levels, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels,
             FT_levels = FT_levels,
             summary = summary_df, n_obs = stan_data$N_obs,
             wall_min = wall_min), meta_path)

cat("=== b / a coefs ===\n")
print(summary_df[grepl("^b[0-9]|^a_", summary_df$variable), ])
cat("\n=== sigmas ===\n")
print(summary_df[grepl("^sigma", summary_df$variable), ])
cat("\nDone.\n")
