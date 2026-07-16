# Finalize mortality/survival comparisons once the production fits land
#
# Runs automatically via SLURM dependency after the survival (survival_100k),
# mortality v2 cch (mort_v2_cch_100k), and v1 paired baseline
# (mort_v1_paired_100k) production fits complete. Produces:
#   1. Survival vs mortality-v2-cch LOO comparison (framing equivalence check;
#      both fit on identical seed-2026 100K data, so loo_compare is valid and
#      the two should be near-identical, confirming the survival framing is a
#      correct sign-flip of the mortality framing).
#   2. PAIRED cch gain: v2-cch vs v1-paired (both seed-2026, identical 100K
#      rows in identical order) via loo_compare. This is the rigorous test of
#      the crown-closure + nonlinear-CR gain (Greg's key covariates).
#   3. Secondary per-observation cch gain vs the older v1 production fit
#      (seed 42), run only if the paired v1 baseline is unavailable.
#   4. Base-rate benchmark hook for the survival model.
#
# Run from Cardinal (fvs-modern so calibration/ resolves):
#   module load gcc/12.3.0 R/4.4.0
#   Rscript calibration/R/finalize_mortality_comparisons.R
#
# Author: Aaron Weiskittel
# Date: 2026-05-21

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(loo)
})

ROOT <- "/users/PUOM0008/crsfaaron/fvs-conus"
MORT <- file.path(ROOT, "output", "conus", "mort")
OUT  <- file.path(MORT, "framing_and_cch_comparison")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load the loo objects the new drivers saved -----------------------
surv_loo_path     <- file.path(MORT, "survival_100k_prod",   "survival_100k_prod_loo.rds")
v2_loo_path       <- file.path(MORT, "v2_cch_100k_prod",     "mort_v2_cch_100k_prod_loo.rds")
v1paired_loo_path <- file.path(MORT, "v1_paired_100k_prod",  "mort_v1_paired_100k_loo.rds")

stopifnot(file.exists(surv_loo_path), file.exists(v2_loo_path))
surv_loo <- readRDS(surv_loo_path)
v2_loo   <- readRDS(v2_loo_path)

n_surv <- nrow(surv_loo$pointwise)
n_v2   <- nrow(v2_loo$pointwise)

cat("=== 1. Survival vs v2 cch (both seed-2026 100K, identical data) ===\n")
cat(sprintf("Survival ELPD: %.1f (SE %.1f), n = %d, per-obs %.4f\n",
            surv_loo$estimates["elpd_loo","Estimate"],
            surv_loo$estimates["elpd_loo","SE"], n_surv,
            surv_loo$estimates["elpd_loo","Estimate"] / n_surv))
cat(sprintf("v2 cch  ELPD: %.1f (SE %.1f), n = %d, per-obs %.4f\n",
            v2_loo$estimates["elpd_loo","Estimate"],
            v2_loo$estimates["elpd_loo","SE"], n_v2,
            v2_loo$estimates["elpd_loo","Estimate"] / n_v2))

if (n_surv == n_v2) {
  cmp_framing <- loo_compare(v2_loo, surv_loo)
  print(cmp_framing)
  saveRDS(cmp_framing, file.path(OUT, "survival_vs_v2cch_loo_compare.rds"))
  cat("\nInterpretation: these should be near-identical (|delta ELPD| within a\n")
  cat("few units), confirming the survival framing is a correct reparameterization\n")
  cat("of the mortality framing with the same predictors (cch + nonlinear CR).\n")
} else {
  cat("\nNOTE: n differs between survival and v2 cch; framing loo_compare skipped.\n")
}

# ---- 2. PAIRED cch gain: v2 cch vs v1 paired (identical seed-2026 data) ---
cat("\n=== 2. PAIRED cch gain: v2 cch vs v1 paired baseline ===\n")
if (file.exists(v1paired_loo_path)) {
  v1p_loo <- readRDS(v1paired_loo_path)
  n_v1p   <- nrow(v1p_loo$pointwise)
  cat(sprintf("v1 paired (no cch) ELPD: %.1f (SE %.1f), n = %d\n",
              v1p_loo$estimates["elpd_loo","Estimate"],
              v1p_loo$estimates["elpd_loo","SE"], n_v1p))
  cat(sprintf("v2 cch            ELPD: %.1f (SE %.1f), n = %d\n",
              v2_loo$estimates["elpd_loo","Estimate"],
              v2_loo$estimates["elpd_loo","SE"], n_v2))
  if (n_v1p == n_v2) {
    cmp_cch <- loo_compare(v1p_loo, v2_loo)
    print(cmp_cch)
    saveRDS(cmp_cch, file.path(OUT, "v2cch_vs_v1paired_loo_compare.rds"))
    delta <- v2_loo$estimates["elpd_loo","Estimate"] -
             v1p_loo$estimates["elpd_loo","Estimate"]
    cat(sprintf("\nPaired delta ELPD (v2 cch - v1 paired): %+.1f\n", delta))
    cat("Positive and large relative to its SE means crown closure at tree tip\n")
    cat("plus quadratic crown ratio improve predictive performance on identical\n")
    cat("data (Greg's key covariate finding confirmed, paired and rigorous).\n")
    paired_done <- TRUE
  } else {
    cat("\nNOTE: n differs (v1 paired ", n_v1p, " vs v2 ", n_v2,
        "); paired loo_compare invalid.\n", sep = "")
    paired_done <- FALSE
  }
} else {
  cat("v1 paired loo not found at ", v1paired_loo_path, "\n", sep = "")
  cat("(Run 34_fit_mortality_v1_paired.R.) Falling back to per-obs contrast.\n")
  paired_done <- FALSE
}

# ---- 3. Secondary per-obs cch gain vs old v1 (seed 42) -------------------
# Only if the paired comparison was not available, since this loads the ~6GB
# v1 production fit and recomputes its loo (different 100K subset, seed 42).
v1_fit_path <- file.path(MORT, "speciesfree", "mort_sf_100k_prod_fit.rds")
if (!paired_done && file.exists(v1_fit_path)) {
  cat("\n=== 3. Secondary per-obs cch gain: v2 cch vs old v1 (seed 42) ===\n")
  cat("Loading old v1 fit (large object, ~6GB)...\n")
  v1_fit <- readRDS(v1_fit_path)
  v1_ll  <- v1_fit$draws("log_lik", format = "draws_matrix")
  v1_loo <- loo(v1_ll)
  saveRDS(v1_loo, file.path(OUT, "v1_mortality_oldseed_loo.rds"))
  rm(v1_fit, v1_ll); gc()

  n_v1 <- nrow(v1_loo$pointwise)
  per_v1 <- v1_loo$estimates["elpd_loo","Estimate"] / n_v1
  per_v2 <- v2_loo$estimates["elpd_loo","Estimate"] / n_v2
  cat(sprintf("old v1 (no cch) per-obs ELPD: %.4f (n = %d, seed 42)\n", per_v1, n_v1))
  cat(sprintf("v2 (cch + nonlinear CR) per-obs ELPD: %.4f (n = %d)\n", per_v2, n_v2))
  cat(sprintf("Per-observation ELPD gain from cch + nonlinear CR: %+.5f\n",
              per_v2 - per_v1))
  cat("Approximate (different subsets); the paired comparison in section 2 is\n")
  cat("the rigorous version when available.\n")

  gain_tbl <- data.frame(
    model = c("old v1 (no cch, seed 42)", "v2 (cch + nonlinear CR)"),
    elpd_loo = c(v1_loo$estimates["elpd_loo","Estimate"],
                 v2_loo$estimates["elpd_loo","Estimate"]),
    n_obs = c(n_v1, n_v2),
    elpd_per_obs = c(per_v1, per_v2)
  )
  fwrite(gain_tbl, file.path(OUT, "cch_gain_per_obs_oldseed.csv"))
} else if (paired_done) {
  cat("\n=== 3. Secondary per-obs cch gain ===\n")
  cat("Skipped: paired comparison in section 2 is available and rigorous.\n")
}

# ---- 4. Base-rate benchmark for the survival model -----------------------
cat("\n=== 4. Base-rate benchmark ===\n")
cat("Run separately: Rscript calibration/R/mortality_base_rate_benchmark.R\n")
cat("The survival driver saves a compact loo plus parameter draws, not the full\n")
cat("fit; the base-rate benchmark should be run from a fit that retains\n")
cat("p_surv_annual, or recomputed from the parameter draws. Flagging as follow-up.\n")

cat("\nDone. Output dir: ", OUT, "\n")
