##=============================================================================
## 32_fit_dg_kuehne_speciesfree.R
## Kuehne et al. (2022) diameter growth model, species-free (B1) variant
##
## Mirrors the May 7 banked v2 fit but drops z_sp from the Stan linear
## predictor. The trait fixed-effect block (W * gamma) is the only species
## signal. If the fit converges with sigma matching v2's banked value, the
## species-free hypothesis is validated on a second base model.
##
## Data shape (matches May 7 prep_meta exactly):
##   - 160K observations after Kuehne filters
##   - 75 species, 11 L1, 21 L2, 84 L3 ecodivision levels
##   - 8 trait columns: wood_specific_gravity, shade_tolerance_num, softwood,
##     leaf_longevity_months, max_ht_m, max_dbh_cm, vulnerability_score,
##     sensitivity
##
## Usage on Cardinal:
##   module load gcc/12.3.0 R/4.4.0
##   cd ~/fvs-modern
##   Rscript --vanilla calibration/R/32_fit_dg_kuehne_speciesfree.R
##
## CLI flags:
##   --stan_file=PATH       override default speciesfree.stan path
##   --outdir=PATH          override default output directory
##   --outname=NAME         override default output basename
##   --subsample=N          fit on a random N-row subsample (default: full)
##   --smoke                run a tiny smoke fit (50 iter warmup + 50 sampling)
##
## Author: Aaron Weiskittel + Claude
## Date: 2026-05-11
##=============================================================================

library(data.table)
library(cmdstanr)
library(posterior)

## CLI parsing ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/dg_kuehne2022_speciesfree.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/dg/speciesfree_pilot")
OUT_NAME  <- get_arg("outname",   "dg_kuehne_cspi_traits1_b1")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")
HOLDOUT_FILE <- get_arg("holdout_spcd_file", NA_character_)

cat("== 32_fit_dg_kuehne_speciesfree.R ==\n")
cat("Stan file: ", STAN_FILE, "\n")
cat("Output:    ", OUT_DIR, "/", OUT_NAME, "_*\n", sep = "")
cat("Subsample: ", if (is.na(SUBSAMPLE)) "(none, full data)" else SUBSAMPLE, "\n")
cat("Smoke:     ", SMOKE, "\n\n")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## 1. Load data --------------------------------------------------------------

DATA_FILE <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- "calibration/traits/species_traits.rds"

stopifnot(file.exists(DATA_FILE), file.exists(STAN_FILE), file.exists(TRAITS_FILE))

cat("Loading", DATA_FILE, "...\n"); flush.console()
dat <- as.data.table(readRDS(DATA_FILE))
cat("  loaded", nrow(dat), "rows x", ncol(dat), "cols\n")

cat("Loading", TRAITS_FILE, "...\n"); flush.console()
traits <- as.data.table(readRDS(TRAITS_FILE))
cat("  loaded", nrow(traits), "species rows\n\n")

## 2. Kuehne data prep -------------------------------------------------------
##    Mirrors the v2 fit. DBH in cm (data is metric), CR fraction in (0, 1].

MIN_OBS_SPECIES <- 5000
CSPI_SHIFT      <- 1.0

trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")

cat("Filtering data for Kuehne fit ...\n"); flush.console()

## Compute annual diameter growth and ancillary transforms.
dat[, dg_obs_periodic := DBH2 - DBH1]
dat[, dg_obs_a := dg_obs_periodic / YEARS]
dat[, sqrt_years := sqrt(YEARS)]
dat[, ln_dbh := log(DBH1)]
dat[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
dat[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
dat[, ln_csi := log(pmax(cspi + CSPI_SHIFT, 0.01))]

## Filters: positive growth tolerance, reasonable values, required columns
## non-missing. dg_obs_a > 0.001 will be enforced inside Stan via lognormal
## skip; we keep small negatives here so that gamma posteriors are not biased.
dat <- dat[
  is.finite(dg_obs_a) &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 &
  is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(cspi) & cspi > -CSPI_SHIFT &
  is.finite(ba_x_rd) & is.finite(bal_x_rd) &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  dg_obs_a > -0.5 & dg_obs_a < 5.0  ## annual growth in cm/yr, reasonable range
]
cat("  after column filters:", nrow(dat), "rows\n")

## Optional held-out species filter (for the generalization validation fit)
if (!is.na(HOLDOUT_FILE) && file.exists(HOLDOUT_FILE)) {
  holdout_spcd <- as.integer(readLines(HOLDOUT_FILE))
  before_n <- nrow(dat)
  dat <- dat[!SPCD %in% holdout_spcd]
  cat("  held-out species removed:", length(holdout_spcd),
      "(dropped", before_n - nrow(dat), "rows)\n")
}

## Species filter: minimum observation threshold
sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("  after species filter (n >=", MIN_OBS_SPECIES, "):", nrow(dat), "rows\n")
cat("  retained species:", nrow(sp_counts), "\n")

## Build integer indices (sorted SPCD, sorted unique EPA codes)
sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))

dat[, sp_idx := match(SPCD, sp_levels)]
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]

cat("  N_sp =", length(sp_levels), " N_L1 =", length(L1_levels),
    " N_L2 =", length(L2_levels), " N_L3 =", length(L3_levels), "\n\n")

## Build trait matrix W (N_sp x P_trait). Standardize each column (z-score).
cat("Building trait matrix W (75 x 8 expected) ...\n"); flush.console()
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
stopifnot(identical(traits_sub$SPCD, sp_levels))

W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
## Impute trait NA with column median (rare but happens for some species).
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) {
    med <- median(W[!na, j], na.rm = TRUE)
    W[na, j] <- med
  }
  ## Standardize
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}
cat("  W dim:", paste(dim(W), collapse = " x "), "\n")
cat("  any remaining NA in W:", any(is.na(W)), "\n\n")

## Optional subsample for the pilot ------------------------------------------
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  cat("  subsampled to", nrow(dat), "rows\n\n")
}

## 3. Build Stan data list ---------------------------------------------------

stan_data <- list(
  N_obs   = nrow(dat),
  N_sp    = length(sp_levels),
  N_L1    = length(L1_levels),
  N_L2    = length(L2_levels),
  N_L3    = length(L3_levels),
  P_trait = ncol(W),

  dg_obs_a       = dat$dg_obs_a,
  sqrt_years     = dat$sqrt_years,

  ln_dbh         = dat$ln_dbh,
  dbh            = dat$DBH1,
  ln_cr_adj      = dat$ln_cr_adj,
  ln_bal_sw_adj  = dat$ln_bal_sw_adj,
  bal_hw         = dat$BAL_HW1,
  ln_csi         = dat$ln_csi,
  ba_x_rd        = dat$ba_x_rd,
  bal_x_rd       = dat$bal_x_rd,

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,

  W = W
)

cat("=== Stan data assembled ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("N_sp    =", stan_data$N_sp, "\n")
cat("N_L1    =", stan_data$N_L1, "\n")
cat("N_L2    =", stan_data$N_L2, "\n")
cat("N_L3    =", stan_data$N_L3, "\n")
cat("P_trait =", stan_data$P_trait, "\n\n")

## 4. Compile + sample -------------------------------------------------------

cat("Compiling Stan model:", STAN_FILE, "\n"); flush.console()
mod <- cmdstan_model(STAN_FILE)

if (SMOKE) {
  iter_warmup   <- 50
  iter_sampling <- 50
  chains        <- 2
  cat("SMOKE TEST: 50 warmup + 50 sampling x 2 chains\n\n")
} else {
  iter_warmup   <- 1000
  iter_sampling <- 1000
  chains        <- 4
  cat("PRODUCTION: 1000 warmup + 1000 sampling x 4 chains\n\n")
}

t_start <- Sys.time()
fit <- mod$sample(
  data            = stan_data,
  chains          = chains,
  parallel_chains = chains,
  iter_warmup     = iter_warmup,
  iter_sampling   = iter_sampling,
  seed            = 42,
  adapt_delta     = 0.9,
  max_treedepth   = 10,
  refresh         = 100
)
t_end <- Sys.time()
wall_min <- as.numeric(difftime(t_end, t_start, units = "mins"))
cat("\nSampling wall time:", round(wall_min, 1), "min\n\n")

## 5. Save fit + meta + summary ----------------------------------------------

fit_path     <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
meta_path    <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summary_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

cat("Saving fit to:", fit_path, "\n"); flush.console()
fit$save_object(fit_path)

cat("Computing summary draws ...\n"); flush.console()
vars_to_summarize <- c(paste0("b", 0:8),
                       paste0("gamma[", seq_len(ncol(W)), "]"),
                       "sigma_L1", "sigma_L2", "sigma_L3", "sigma",
                       paste0("trait_effect[", seq_len(length(sp_levels)), "]"))
summary_df <- fit$summary(variables = vars_to_summarize,
                          "mean", "median", "sd", "mad",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail")
data.table::fwrite(summary_df, summary_path)

prep_meta <- list(
  sp         = sp_levels,
  L1         = L1_levels,
  L2         = L2_levels,
  L3         = L3_levels,
  cspi_shift = CSPI_SHIFT
)
meta <- list(
  form       = "kuehne_speciesfree_b1",
  site_var   = "cspi",
  traits     = TRUE,
  trait_cols = trait_cols,
  stan_file  = STAN_FILE,
  prep_meta  = prep_meta,
  summary    = summary_df,
  n_obs      = stan_data$N_obs,
  wall_min   = wall_min
)
saveRDS(meta, meta_path)
cat("Saved meta to:    ", meta_path, "\n")
cat("Saved summary to: ", summary_path, "\n\n")

## 6. Sigma callout (the manuscript decision number) ------------------------

sigma_summary <- fit$summary("sigma", "mean", "sd",
                             ~quantile(.x, c(0.05, 0.95)), "rhat")
cat("=== DG_Kuehne B1 sigma posterior ===\n")
print(sigma_summary)
cat("\nFor manuscript comparison, the v2 banked sigma was approximately 21.3\n")
cat("(see calibration/output/conus/dg/dg_kuehne_cspi_traits1_summary.csv).\n")
cat("If B1 sigma is within +/- 5 percent of v2, species-free is validated\n")
cat("on the second base model and the trait-driven story is robust.\n")

cat("\n=== Done ===\n")
