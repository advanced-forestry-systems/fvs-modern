# Compare survival-framed vs mortality-framed species-free models
#
# Confirms empirically that the two framings give identical LOO ELPD and
# per-observation predictions (they share the same likelihood with the
# exposure offset). Any difference beyond Monte Carlo noise indicates a
# parameterization or prior issue worth investigating.
#
# Run from Cardinal:
#   module load gdal/3.7.3 gcc/12.3.0 geos/3.12.0 proj/9.2.1 R/4.4.0
#   Rscript R/compare_survival_vs_mortality.R \
#       --surv=calibration/output/conus/mort/survival_speciesfree/survival_speciesfree_fit.rds \
#       --mort=calibration/output/conus/mort/mort_prod/mort_speciesfree_fit.rds
#
# Author: Aaron Weiskittel
# Date: 2026-05-20

suppressPackageStartupMessages({
  library(cmdstanr)
  library(loo)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

SURV_FIT <- get_arg("surv", "calibration/output/conus/mort/survival_speciesfree/survival_speciesfree_fit.rds")
MORT_FIT <- get_arg("mort", "calibration/output/conus/mort/mort_prod/mort_speciesfree_fit.rds")
OUT_DIR  <- get_arg("outdir", "calibration/output/conus/mort/framing_comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

surv <- readRDS(SURV_FIT)
mort <- readRDS(MORT_FIT)

# ---- 1. LOO ELPD for each ------------------------------------------------
loo_surv <- loo(surv$draws("log_lik"))
loo_mort <- loo(mort$draws("log_lik"))

cat("=== LOO ELPD ===\n")
cat(sprintf("Survival framing : %.1f (SE %.1f)\n",
            loo_surv$estimates["elpd_loo", "Estimate"],
            loo_surv$estimates["elpd_loo", "SE"]))
cat(sprintf("Mortality framing: %.1f (SE %.1f)\n",
            loo_mort$estimates["elpd_loo", "Estimate"],
            loo_mort$estimates["elpd_loo", "SE"]))

cmp <- loo_compare(loo_mort, loo_surv)
print(cmp)
saveRDS(list(loo_surv = loo_surv, loo_mort = loo_mort, compare = cmp),
        file.path(OUT_DIR, "framing_loo_compare.rds"))

# ---- 2. Per-observation predicted survival agreement --------------------
# Only meaningful if both fits used the same data and observation order.
if ("p_surv_annual" %in% surv$metadata()$stan_variables &&
    "p_mort_annual" %in% mort$metadata()$stan_variables) {
  ps <- colMeans(surv$draws("p_surv_annual", format = "draws_matrix"))
  pm <- colMeans(mort$draws("p_mort_annual", format = "draws_matrix"))
  n <- min(length(ps), length(pm))
  # survival should equal 1 - mortality if framings are equivalent and aligned
  agreement <- data.frame(
    p_surv_from_survival_model   = ps[seq_len(n)],
    p_surv_from_mortality_model  = 1 - pm[seq_len(n)]
  )
  agreement$abs_diff <- abs(agreement$p_surv_from_survival_model -
                            agreement$p_surv_from_mortality_model)
  fwrite(agreement, file.path(OUT_DIR, "per_obs_survival_agreement.csv"))
  cat(sprintf("\nPer-observation annual survival agreement:\n"))
  cat(sprintf("  mean abs diff = %.5f, max abs diff = %.5f\n",
              mean(agreement$abs_diff), max(agreement$abs_diff)))
  cat("  (Differences should be near zero if data order matches; large diffs\n")
  cat("   usually mean the two fits used different subsamples, not a model issue.)\n")
}

# ---- 3. Coefficient sign check ------------------------------------------
# Pull crown ratio and competition coefficients from each, confirm signs flip.
get_coef <- function(fit, var) {
  s <- fit$summary(variables = var, "mean")
  s$mean[1]
}
cat("\n=== Coefficient sign check (should be opposite signs) ===\n")
for (v in c("b3", "b5", "b6")) {
  cat(sprintf("  %s: survival = %+.4f, mortality = %+.4f\n",
              v, get_coef(surv, v), get_coef(mort, v)))
}

cat("\nConclusion guidance:\n")
cat("  If LOO ELPD matches within a fraction of an SE and coefficient signs are\n")
cat("  cleanly opposite, the two framings are equivalent. Adopt the survival\n")
cat("  framing for the manuscript (more interpretable, matches Greg's framing).\n")
cat("\nOutput dir: ", OUT_DIR, "\n")
