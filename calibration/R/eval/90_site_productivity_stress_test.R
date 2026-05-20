# =============================================================================
# Title: Site productivity stress test for HG and DG_Kuehne base models
# Author: A. Weiskittel
# Date: 2026-05-13
# Description: Refits HG (or DG_Kuehne) under three different site productivity
#              variables: composite site productivity index (cspi), climate
#              site index (climate_si), and biomass growth index (bgi). All
#              other model architecture is identical. Compares sigma, gamma,
#              site coefficient, and posterior predictive accuracy across the
#              three site metrics.
#
# This is the first systematic site-variable comparison for FVS-CONUS. All
# eight production fits to date have used cspi. The manuscript can claim
# robustness to site metric choice if and only if the three fits produce
# comparable sigma and consistent direction on disturbance modifiers.
#
# CLI flags:
#   --model={hg, dg_kue}        target base model (default: hg)
#   --site_var={cspi, bgi, climate_si}  the site productivity variable
#   --subsample=N               default 500K for HG, 200K for DG_Kue
#   --outdir=PATH               default calibration/output/conus/<model>/site_stress
#
# Run on Cardinal (one fit per site variable; submit 3 jobs):
#   for SV in cspi bgi climate_si; do
#     sbatch calibration/slurm/submit_hg_site_stress.sh $SV
#   done
#
# Outputs:
#   <outdir>/<model>_<site_var>_fit.rds
#   <outdir>/<model>_<site_var>_meta.rds
#   <outdir>/<model>_<site_var>_summary.csv
#
# Then 90b_site_stress_compare.R reads the three banked fits and writes a
# comparison table.
# =============================================================================

library(data.table)
library(cmdstanr)
library(posterior)

# --- CLI parsing -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

MODEL    <- get_arg("model", "hg")
VARIANT  <- get_arg("variant", "b1")
SITE_VAR <- get_arg("site_var", "cspi")
stopifnot(SITE_VAR %in% c("cspi", "bgi", "climate_si"))
stopifnot(MODEL %in% c("hg", "dg_kue"))
stopifnot(VARIANT %in% c("b1", "b2"))

PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR   <- get_arg("outdir",
                     file.path(PROJ_ROOT,
                               paste0("calibration/output/conus/", MODEL,
                                      "/site_stress")))
SUBSAMPLE <- as.integer(get_arg("subsample",
                                if (MODEL == "hg") "500000" else "200000"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Site productivity stress test\n"))
cat(sprintf("  model:     %s\n", MODEL))
cat(sprintf("  variant:   %s\n", VARIANT))
cat(sprintf("  site_var:  %s\n", SITE_VAR))
cat(sprintf("  outdir:    %s\n", OUT_DIR))
cat(sprintf("  subsample: %d\n\n", SUBSAMPLE))

# --- Load data --------------------------------------------------------------

DATA_FILE   <- file.path(PROJ_ROOT,
  "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- file.path(PROJ_ROOT, "calibration/traits/species_traits.rds")

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done\n\n")

# --- Apply filters ---------------------------------------------------------
SHIFT_CSPI <- 1.0
SHIFT_BGI  <- 1.0  # bgi is non-negative
SHIFT_CSI  <- 0    # climate_si is positive

# Pick the site variable and apply the right shift before log
site_col <- SITE_VAR
shift    <- switch(SITE_VAR, cspi = SHIFT_CSPI, bgi = SHIFT_BGI, climate_si = SHIFT_CSI)

dat[, site_value := .SD[[1]], .SDcols = site_col]
dat[, ln_site := log(pmax(site_value + shift, 0.01))]

dat <- dat[
  is.finite(DBH1) & DBH1 >= 2.54 & is.finite(DBH2) &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(site_value) & site_value > -shift
]

# Model-specific filters
if (MODEL == "dg_kue") {
  MIN_OBS_SPECIES <- 5000
  dat <- dat[is.finite(BAL_SW1) & BAL_SW1 >= 0 &
             is.finite(BAL_HW1) & BAL_HW1 >= 0 &
             is.finite(ba_x_rd) & is.finite(bal_x_rd)]
  dat[, dg_obs_a := (DBH2 - DBH1) / YEARS]
  dat[, sqrt_years := sqrt(YEARS)]
  dat[, ln_dbh := log(DBH1)]
  dat[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
  dat[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
  dat <- dat[dg_obs_a > -0.5 & dg_obs_a < 5.0]
} else if (MODEL == "hg") {
  MIN_OBS_SPECIES <- 5000
  # HG uses simpler covariates: DBH, CR, BAL, site
  dat <- dat[is.finite(HT1) & is.finite(HT2) & HT2 > HT1 &
             is.finite(BAL1) & BAL1 >= 0]
  dat[, hg_obs_a := (HT2 - HT1) / YEARS]
  dat[, sqrt_years := sqrt(YEARS)]
  dat[, ln_dbh := log(DBH1)]
  dat[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
  dat <- dat[hg_obs_a > 0.0 & hg_obs_a < 2.0]
}
cat("After column filters:", nrow(dat), "rows\n")

# Species filter
sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("After species filter (n >=", MIN_OBS_SPECIES, "):", nrow(dat), "rows\n")
cat("Retained species:", nrow(sp_counts), "\n")

# Subsample
set.seed(42)
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  cat("Subsampled to:", nrow(dat), "rows\n\n")
}

# --- HG covariate prep (parallel to Kuehne, different covariate set) -------
if (MODEL == "hg") {
  # HG ORGANON form uses: ln_dbh, dbh, ln_cr_adj, bal, site, plus interactions
  MIN_OBS_SPECIES <- 5000
  dat <- dat[is.finite(HT1) & is.finite(HT2) & HT2 > HT1 &
             is.finite(BAL1) & BAL1 >= 0]
  dat[, hg_obs_a   := (HT2 - HT1) / YEARS]
  dat[, sqrt_years := sqrt(YEARS)]
  dat[, ln_dbh     := log(DBH1)]
  dat[, ln_cr_adj  := log((CR1 + 0.2) / 1.2)]
  dat[, bal_m      := BAL1 * 0.2296]   # ft2/ac to m2/ha
  dat[, bal_log    := log(bal_m + 5.0)]
  dat <- dat[hg_obs_a > 0.0 & hg_obs_a < 2.0]
  cat("After HG column filters:", nrow(dat), "rows\n")

  sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
  dat <- dat[SPCD %in% sp_counts$SPCD]
  cat("After HG species filter:", nrow(dat), "rows;",
      nrow(sp_counts), "species\n")
}

# --- Build trait matrix W ---------------------------------------------------
trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
sp_levels <- sort(unique(dat$SPCD))
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols),
                     with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}

# --- Indices ----------------------------------------------------------------
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))

dat[, sp_idx := match(SPCD, sp_levels)]
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]

# --- Stan data + Stan file --------------------------------------------------
if (MODEL == "dg_kue") {
  STAN_FILE <- file.path(PROJ_ROOT,
                         "calibration/stan/dg_kuehne2022_speciesfree.stan")
  stan_data <- list(
    N_obs = nrow(dat), N_sp = length(sp_levels),
    N_L1 = length(L1_levels), N_L2 = length(L2_levels), N_L3 = length(L3_levels),
    P_trait = ncol(W),
    dg_obs_a = dat$dg_obs_a, sqrt_years = dat$sqrt_years,
    ln_dbh = dat$ln_dbh, dbh = dat$DBH1, ln_cr_adj = dat$ln_cr_adj,
    ln_bal_sw_adj = dat$ln_bal_sw_adj, bal_hw = dat$BAL_HW1,
    ln_csi = dat$ln_site,
    ba_x_rd = dat$ba_x_rd, bal_x_rd = dat$bal_x_rd,
    sp_idx = dat$sp_idx, L1_idx = dat$L1_idx,
    L2_idx = dat$L2_idx, L3_idx = dat$L3_idx,
    W = W
  )
} else if (MODEL == "hg") {
  # HG B1 uses the pre-compiled hg_organon_speciesfree binary; the source was
  # lost in git clean. Use the binary directly via cmdstanr::cmdstan_model
  # with the exe_file argument when ready_to_compile = FALSE. For now, fall
  # back to recompiling from a copy of the dg_kuehne speciesfree stan if the
  # HG source is unavailable.
  STAN_FILE <- if (VARIANT == "b1") {
    file.path(PROJ_ROOT, "calibration/stan/hg_organon_speciesfree.stan")
  } else {
    file.path(PROJ_ROOT, "calibration/stan/height_increment.stan")
  }
  if (!file.exists(STAN_FILE)) {
    stop(sprintf("HG Stan file not found: %s", STAN_FILE))
  }

  # The HG ORGANON form linear predictor uses:
  #   eta = a0 + trait_effect[sp] + a1 ln(dbh) + a2 dbh + a3 ln_cr_adj
  #       + a4 bal_log + a5 ln_csi + a6 ln(dbh)^2 + a7 dbh*ln_csi + a8 ln_csi^2
  # but the exact mapping depends on which version of the Stan file is used.
  # Build the data block conservatively with the same names the Stan model
  # expects. If the Stan file uses different field names, this call will fail
  # at sampling with a clear error and you can edit the data block.
  stan_data <- list(
    N_obs = nrow(dat), N_sp = length(sp_levels),
    N_L1 = length(L1_levels), N_L2 = length(L2_levels), N_L3 = length(L3_levels),
    P_trait = ncol(W),
    hg_obs_a = dat$hg_obs_a, sqrt_years = dat$sqrt_years,
    ln_dbh = dat$ln_dbh, dbh = dat$DBH1, ln_cr_adj = dat$ln_cr_adj,
    bal_log = dat$bal_log, ln_csi = dat$ln_site,
    sp_idx = dat$sp_idx, L1_idx = dat$L1_idx,
    L2_idx = dat$L2_idx, L3_idx = dat$L3_idx,
    W = W
  )
}

cat("=== Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("N_sp    =", stan_data$N_sp, "\n")
cat("ln_site mean / sd =", round(mean(stan_data$ln_csi), 3),
    "/", round(sd(stan_data$ln_csi), 3), "\n\n")

# --- Compile and sample -----------------------------------------------------

cat("Compiling Stan model:", STAN_FILE, "\n"); flush.console()
mod <- cmdstan_model(STAN_FILE)

t_start <- Sys.time()
fit <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  adapt_delta     = 0.9,
  max_treedepth   = 10,
  refresh         = 100
)
wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("\nSampling wall time:", round(wall_min, 1), "min\n\n")

# --- Save -------------------------------------------------------------------
OUT_NAME <- sprintf("%s_%s_b1", MODEL, SITE_VAR)
fit_path     <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
meta_path    <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summary_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

cat("Saving fit ..."); flush.console()
fit$save_object(fit_path)
cat(" done\n")

vars <- c(paste0("b", 0:8), paste0("gamma[", seq_len(ncol(W)), "]"),
          "sigma_L1", "sigma_L2", "sigma_L3", "sigma")
summary_df <- fit$summary(variables = vars,
                          "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail") |>
  data.table::as.data.table() |>
  data.table::setnames(c("5%", "95%"), c("q5", "q95"), skip_absent = TRUE)
data.table::fwrite(summary_df, summary_path)

meta <- list(
  form       = MODEL,
  variant    = "b1_speciesfree",
  site_var   = SITE_VAR,
  trait_cols = trait_cols,
  stan_file  = STAN_FILE,
  prep_meta  = list(sp = sp_levels, L1 = L1_levels, L2 = L2_levels,
                    L3 = L3_levels, site_shift = shift),
  summary    = summary_df,
  n_obs      = stan_data$N_obs,
  wall_min   = wall_min
)
saveRDS(meta, meta_path)
cat("Saved meta and summary\n\n")

# --- Sigma callout ---------------------------------------------------------
sigma_summary <- fit$summary("sigma", "mean", "sd",
                             ~quantile(.x, c(0.05, 0.95)), "rhat")
cat(sprintf("=== %s, site_var=%s sigma posterior ===\n", MODEL, SITE_VAR))
print(sigma_summary)
cat("\nDone.\n")
