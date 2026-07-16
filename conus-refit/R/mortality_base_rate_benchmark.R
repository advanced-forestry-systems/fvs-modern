# Base-rate benchmark for our hierarchical mortality model
#
# Matches Greg Johnson's Kahneman-Tversky base-rate framing so our species-free
# hierarchical model can be compared apples-to-apples against (a) the species
# base rate and (b) Greg's per-species cr+cch Gompit model.
#
# The base rate for each species is its average annual mortality rate, ignoring
# all other covariates. Any useful model must beat the base rate. We compute
# per-species negative log likelihood for:
#   1. Base rate (species mean annual mortality, period-length adjusted)
#   2. Our hierarchical model (posterior mean prediction per observation)
# and report the NLL reduction per species, mirroring Greg's Figure 2.
#
# Run from Cardinal:
#   module load gdal/3.7.3 gcc/12.3.0 geos/3.12.0 proj/9.2.1 R/4.4.0
#   Rscript R/mortality_base_rate_benchmark.R \
#       --fit=calibration/output/conus/mort/mort_prod/mort_speciesfree_fit.rds
#
# Author: Aaron Weiskittel
# Date: 2026-05-20

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

FIT_FILE  <- get_arg("fit", "calibration/output/conus/mort/mort_prod/mort_speciesfree_fit.rds")
DATA_FILE <- get_arg("data", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_DIR   <- get_arg("outdir", "calibration/output/conus/mort/base_rate_benchmark")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- as.data.table(readRDS(DATA_FILE))

# ---- 1. Species base rate (period-length adjusted) ----------------------
# Base rate annual survival per species: geometric over observed intervals.
# annual_surv_sp = (sum of survivals weighted by exposure) approach:
# We estimate the annual mortality hazard h_sp such that the per-observation
# survival exp(-h_sp * T) best matches the observed alive/dead pattern.
# Closed form MLE for constant hazard: h = deaths / total_exposure_years.
base_rate <- dat[, .(
  deaths        = sum(alive == 0),
  exposure_yrs  = sum(T_years),
  n             = .N
), by = SPCD]
base_rate[, hazard_sp := deaths / exposure_yrs]
base_rate[, annual_surv_sp := exp(-hazard_sp)]

dat <- merge(dat, base_rate[, .(SPCD, hazard_sp)], by = "SPCD", all.x = TRUE)

# Base-rate per-observation log likelihood
dat[, log_p_surv_base := -hazard_sp * T_years]
dat[, ll_base := fifelse(alive == 1, log_p_surv_base,
                         log(1 - exp(log_p_surv_base)))]

# ---- 2. Hierarchical model per-observation log likelihood ----------------
# The fit already stores log_lik per observation in generated quantities.
# Pull the posterior mean log_lik per observation.
fit <- readRDS(FIT_FILE)
ll_draws <- fit$draws("log_lik", format = "draws_matrix")
# Posterior mean log pointwise predictive density per observation
dat[, ll_model := colMeans(ll_draws)[seq_len(nrow(dat))]]

# ---- 3. Per-species NLL comparison --------------------------------------
per_species <- dat[, .(
  n            = .N,
  nll_base     = -sum(ll_base,  na.rm = TRUE),
  nll_model    = -sum(ll_model, na.rm = TRUE)
), by = SPCD]
per_species[, nll_reduction := nll_base - nll_model]
per_species[, pct_reduction := 100 * nll_reduction / nll_base]
per_species[, model_wins := nll_model < nll_base]
setorder(per_species, -nll_reduction)

fwrite(per_species, file.path(OUT_DIR, "per_species_nll_comparison.csv"))

# ---- 4. Headline summary -------------------------------------------------
n_species   <- nrow(per_species)
n_wins      <- sum(per_species$model_wins)
total_base  <- sum(per_species$nll_base)
total_model <- sum(per_species$nll_model)
overall_pct <- 100 * (total_base - total_model) / total_base

cat("=== Mortality base-rate benchmark ===\n")
cat(sprintf("Species evaluated:           %d\n", n_species))
cat(sprintf("Species where model beats base rate: %d (%.1f%%)\n",
            n_wins, 100 * n_wins / n_species))
cat(sprintf("Total NLL base rate:         %.0f\n", total_base))
cat(sprintf("Total NLL hierarchical model:%.0f\n", total_model))
cat(sprintf("Overall NLL reduction:       %.1f%%\n", overall_pct))

# ---- 5. Species where model loses (compare to Greg's Table 2) ------------
losers <- per_species[model_wins == FALSE]
fwrite(losers, file.path(OUT_DIR, "species_model_loses.csv"))
cat(sprintf("\nSpecies where hierarchical model does NOT beat base rate: %d\n", nrow(losers)))
if (nrow(losers) > 0) {
  cat("These should be cross-checked against Greg's Table 2 (his 19 degraded species).\n")
  cat("If our hierarchical model beats base rate on species where Greg's per-species\n")
  cat("model degraded, that demonstrates the species-free thesis: borrowing strength\n")
  cat("through traits + ecoregion handles low-data species better than per-species fits.\n")
  print(head(losers[order(SPCD)], 25))
}

# ---- 6. Greg comparison hook --------------------------------------------
# If Greg's per-species NLL table is available, join and compare directly.
GREG_NLL <- get_arg("greg_nll", NA_character_)
if (!is.na(GREG_NLL) && file.exists(GREG_NLL)) {
  greg <- as.data.table(read.csv(GREG_NLL))
  # Expect columns SPCD and nll_revised (Greg's revised model per-species NLL)
  cmp <- merge(per_species[, .(SPCD, nll_model)],
               greg[, .(SPCD, nll_greg = nll_revised)], by = "SPCD")
  cmp[, ours_better := nll_model < nll_greg]
  fwrite(cmp, file.path(OUT_DIR, "ours_vs_greg_per_species.csv"))
  cat(sprintf("\nHead-to-head vs Greg: ours better on %d of %d species (%.1f%%)\n",
              sum(cmp$ours_better), nrow(cmp),
              100 * mean(cmp$ours_better)))
} else {
  cat("\nGreg per-species NLL table not provided (--greg_nll=...).\n")
  cat("Once available, this script will produce the direct head-to-head comparison.\n")
}

cat("\nOutput dir: ", OUT_DIR, "\n")
