##=============================================================================
## 35_fit_ingrowth_negbinom_v4.R
##
## Ingrowth driver, v4 architecture:
##   - drops stand_age (unreliable in FIA; many -999 values)
##   - adds HT40_1 (top height) as the stand structure metric
##   - switches from RD1 (Curtis) to rd_sdimax (SDI / SDImax_brms)
##   - keeps ln_BA, ln_BAL, ln_CSI, clim_pca1
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

STAN_FILE <- get_arg("stan_file", "calibration/stan/ingrowth_negbinom_v2.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/ingrowth")
OUT_NAME  <- get_arg("outname",   "ingrowth_negbinom_v4")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")
STATES_FILTER <- get_arg("states", NA_character_)  # comma-separated, e.g., "AR,CA"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 35_fit_ingrowth_negbinom_v4.R ==\n")
cat("Stan:", STAN_FILE, "\n\n")

## 1. Load GRM, identify INGROWTH events ------------------------------------
GRM_DIR   <- "calibration/data/raw_fia"
grm_files <- list.files(GRM_DIR, pattern = "_TREE_GRM_COMPONENT\\.csv$",
                        full.names = TRUE)
if (!is.na(STATES_FILTER)) {
  states <- trimws(strsplit(STATES_FILTER, ",")[[1]])
  grm_files <- grm_files[basename(grm_files) %in% paste0(states, "_TREE_GRM_COMPONENT.csv")]
}
cat("GRM files:", length(grm_files), "\n")

gr_cols <- c("TRE_CN", "PLT_CN", "STATECD", "SUBP_COMPONENT_AL_FOREST")
cat("Loading GRM ..."); flush.console()
grm <- rbindlist(lapply(grm_files, fread, select = gr_cols, showProgress = FALSE))
cat(" done. Rows:", nrow(grm), "\n")

ingrowth <- grm[SUBP_COMPONENT_AL_FOREST == "INGROWTH"]
plot_counts <- ingrowth[, .(n_recruits = .N), by = PLT_CN]
all_plots <- unique(grm[, .(PLT_CN, STATECD)])
all_plots <- merge(all_plots, plot_counts, by = "PLT_CN", all.x = TRUE)
all_plots[is.na(n_recruits), n_recruits := 0L]
rm(grm, ingrowth); gc()
cat("Plot-level: total =", nrow(all_plots),
    " with ingrowth =", sum(all_plots$n_recruits > 0), "\n\n")

## 2. Load plot covariates from matched-pairs file --------------------------
PAIRS_PATH <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
cat("Loading pairs ..."); flush.console()
trees <- as.data.table(readRDS(PAIRS_PATH))
cat(" done\n")

# Plot-level aggregation. Use HT40_1 / HTlorey1 for top height, rd_sdimax for
# relative density. Drop STDAGE (unreliable, -999 values).
plot_covs <- trees[, .(
  ba_t1       = first(BA1) * 0.2296,
  bal_mean_t1 = mean(BAL1, na.rm = TRUE) * 0.2296,
  ht40_t1     = first(HT40_1) * 0.3048,        # feet to meters
  htlorey_t1  = first(HTlorey1) * 0.3048,
  rd_sdimax   = first(rd_sdimax),
  rd_add      = first(rd_add),
  years       = first(YEARS),
  cspi        = first(cspi),
  bgi         = first(bgi),
  climate_si  = first(climate_si),
  clim_pca1   = first(clim_pca1),
  EPA_L1_CODE = first(EPA_L1_CODE),
  EPA_L2_CODE = first(EPA_L2_CODE),
  EPA_L3_CODE = first(EPA_L3_CODE),
  dom_spcd    = SPCD[which.max(DBH1 * (TREESTATUS1 == 1))]
), by = PLT_CN_cond1]
setnames(plot_covs, "PLT_CN_cond1", "PLT_CN")
plot_covs[, PLT_CN := as.numeric(PLT_CN)]
all_plots[, PLT_CN := as.numeric(PLT_CN)]

dat <- merge(all_plots, plot_covs, by = "PLT_CN")
cat("After inner join:", nrow(dat), "plots\n")

## 3. Filter to plots with complete data ----------------------------------
dat <- dat[
  is.finite(years) & years >= 1 & years <= 20 &
  is.finite(ba_t1) & ba_t1 >= 0 &
  is.finite(rd_sdimax) & rd_sdimax > 0 & rd_sdimax < 3.0 &
  is.finite(ht40_t1) & ht40_t1 > 0 &
  is.finite(climate_si) & climate_si > 0 &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
  !is.na(EPA_L2_CODE) & EPA_L2_CODE != "" &
  !is.na(EPA_L3_CODE) & EPA_L3_CODE != "" &
  !is.na(dom_spcd)
]
dat[is.na(bal_mean_t1), bal_mean_t1 := 0]
dat[is.na(clim_pca1), clim_pca1 := 0]

dat[, ln_ba    := log(ba_t1 + 1)]
## DECORRELATE: log(BAL_mean+1) is collinear with ln_ba (offsetting +4.81/-4.39 in the prior fit,
## giving a backwards positive density->recruitment response). Replace with the BAL/BA ratio
## (relative share of BA in larger trees), which is not collinear with total BA. Stan model
## unchanged: b2 now multiplies this ratio (still passed in the ln_bal data slot).
dat[, ln_bal   := bal_mean_t1 / pmax(ba_t1, 0.01)]
dat[, ln_csi   := log(pmax(climate_si, 0.01))]
dat[, ln_ht40  := log(ht40_t1 + 1)]
dat[, log_years := log(years)]

cat("After filtering:", nrow(dat), "plots\n")
cat("Mean ingrowth/plot:", round(mean(dat$n_recruits), 2),
    " | rd_sdimax median:", round(median(dat$rd_sdimax), 2),
    " | ht40 median:", round(median(dat$ht40_t1), 1), "m\n\n")

## 4. Build trait matrix W_dom + indices -----------------------------------
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")
traits <- as.data.table(readRDS(TRAITS_FILE))
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
for (col in trait_cols) {
  med <- median(traits[[col]], na.rm = TRUE)
  traits[is.na(get(col)), (col) := med]
}
trait_means <- sapply(traits[, trait_cols, with = FALSE], mean)
trait_sds   <- sapply(traits[, trait_cols, with = FALSE], sd)

# DF variety split: dom_spcd 202 -> 2020 (Coastal, EPA_L1=7) or 2021 (Rocky Mt, else)
# Activates only when species_traits_v2.rds is in use.
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  pre_df_n <- sum(dat$dom_spcd == 202L)
  dat[dom_spcd == 202L & as.character(EPA_L1_CODE) == "7", dom_spcd := 2020L]
  dat[dom_spcd == 202L, dom_spcd := 2021L]
  cat(sprintf("DF variety split (plot dom_spcd): %d -> %d coastal (2020) + %d rocky (2021)\n",
              pre_df_n, sum(dat$dom_spcd == 2020L), sum(dat$dom_spcd == 2021L)))
}

# Lodgepole pine variety split on dom_spcd
if (any(traits$SPCD == 1080L) && any(traits$SPCD == 1081L)) {
  pre_lp_n <- sum(dat$dom_spcd == 108L)
  dat[dom_spcd == 108L & as.character(EPA_L1_CODE) == "7", dom_spcd := 1080L]
  dat[dom_spcd == 108L, dom_spcd := 1081L]
  cat(sprintf("Lodgepole variety split (plot dom_spcd): %d -> %d shore (1080) + %d rocky (1081)\n",
              pre_lp_n, sum(dat$dom_spcd == 1080L), sum(dat$dom_spcd == 1081L)))
}

dat <- merge(dat, traits[, c("SPCD", trait_cols), with = FALSE],
             by.x = "dom_spcd", by.y = "SPCD", all.x = TRUE)
for (col in trait_cols) {
  med <- median(dat[[col]], na.rm = TRUE)
  dat[is.na(get(col)), (col) := med]
}
W_dom <- as.matrix(dat[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W_dom))) {
  W_dom[, j] <- (W_dom[, j] - trait_means[j]) / trait_sds[j]
}

L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]; W_dom <- W_dom[idx, , drop = FALSE]
  cat("Subsampled to:", nrow(dat), "plots\n\n")
}

## 5. Stan data -----------------------------------------------------------
stan_data <- list(
  N_plots   = nrow(dat),
  N_L1      = length(L1_levels),
  N_L2      = length(L2_levels),
  N_L3      = length(L3_levels),
  P_trait   = ncol(W_dom),
  n_recruits = as.integer(dat$n_recruits),
  log_years  = dat$log_years,
  ln_ba      = dat$ln_ba,
  ln_bal     = dat$ln_bal,
  rd_sdimax  = dat$rd_sdimax,
  ln_csi     = dat$ln_csi,
  ln_ht40    = dat$ln_ht40,
  clim_pca1  = dat$clim_pca1,
  L1_idx     = dat$L1_idx,
  L2_idx     = dat$L2_idx,
  L3_idx     = dat$L3_idx,
  W_dom      = W_dom
)
cat("=== Stan data ===\n")
cat("N_plots   =", stan_data$N_plots, "\n")
cat("N_L1/L2/L3=", stan_data$N_L1, "/", stan_data$N_L2, "/", stan_data$N_L3, "\n")
cat("Mean recruits =", round(mean(stan_data$n_recruits), 2), "\n")
cat("Var/Mean overdisp =", round(var(stan_data$n_recruits) / mean(stan_data$n_recruits), 2), "\n\n")

## 6. Compile + sample ---------------------------------------------------
cat("Compiling:", STAN_FILE, "\n"); flush.console()
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

## 7. Save ----------------------------------------------------------------
fit_path  <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
meta_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

fit$save_object(fit_path)
vars <- c(paste0("b", 0:6), paste0("gamma[", seq_len(ncol(W_dom)), "]"),
          "sigma_L1", "sigma_L2", "sigma_L3", "phi")
summary_df <- fit$summary(variables = vars, "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail")
names(summary_df)[names(summary_df) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summary_df, summ_path)
saveRDS(list(form = "negbinom_v2", trait_cols = trait_cols,
             trait_means = trait_means, trait_sds = trait_sds,  # saved for out-of-sample prediction
             stan_file = STAN_FILE, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels,
             summary = summary_df, n_plots = stan_data$N_plots,
             wall_min = wall_min), meta_path)

cat("=== b coefficients ===\n")
print(summary_df[grepl("^b[0-9]", summary_df$variable), ])
cat("\n=== gamma trait coefs ===\n")
print(summary_df[grepl("^gamma", summary_df$variable), ])
cat("\nDone.\n")
