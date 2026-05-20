#!/usr/bin/env Rscript
## compare_b1_b2_dg_kue.R
##
## Decision-quality comparison between the species-free B1 DG_Kuehne pilot
## (dg_kuehne2022_speciesfree.stan, 200K subsample) and the production
## species-aware B2 DG_Kuehne fit (dg_kuehne2022_lognormal.stan, banked
## May 7). Mirrors compare_b1_b2_hg.R structure.
##
## Outputs:
##   1. Shared-parameter posterior comparison (b0..b8, gamma, sigma).
##   2. sigma comparison (the manuscript decision number).
##   3. sigma_sp recap from B2 (the residual species effect B1 dropped).
##   4. gamma posterior comparison.
##
## Reference v2 banked values (from May 7 summary CSV):
##   sigma:     1.913 [1.907, 1.918]   <- comparison target for B1
##   sigma_sp:  0.207 [0.179, 0.241]   <- residual species effect (should be small)
##
## Usage on Cardinal (login node is fine, both fits fit in 32G):
##   module load gcc/12.3.0 R/4.4.0
##   cd ~/fvs-modern
##   Rscript calibration/R/compare_b1_b2_dg_kue.R

suppressPackageStartupMessages({
  library(data.table)
  library(posterior)
})

b1_path <- "~/fvs-modern/calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds"
b2_path <- "~/fvs-modern/calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds"
b1_path <- path.expand(b1_path); b2_path <- path.expand(b2_path)

stopifnot(file.exists(b1_path), file.exists(b2_path))

cat("FVS-CONUS B1 (species-free) vs B2 (species-aware) DG_Kuehne comparison\n")
cat("======================================================================\n\n")

read_fit <- function(p, label) {
  cat(sprintf("Reading %s (%s) ... ", label, basename(p))); flush.console()
  fit <- readRDS(p)
  draws <- tryCatch(fit$draws(), error = function(e) as_draws_array(fit))
  cat(sprintf("draws array %s\n", paste(dim(draws), collapse = " x ")))
  list(label = label, fit = fit, draws = draws)
}

B1 <- read_fit(b1_path, "B1")
B2 <- read_fit(b2_path, "B2")

vars_B1 <- variables(B1$draws)
vars_B2 <- variables(B2$draws)

## 1. Shared fixed coefficients ----------------------------------------------

shared <- intersect(vars_B1, vars_B2)
focal <- c(paste0("b", 0:8), grep("^gamma", vars_B1, value = TRUE),
           "sigma", "sigma_L1", "sigma_L2", "sigma_L3")
focal <- intersect(focal, shared)

if (length(focal) > 0) {
  cat("\n=== Shared parameters: B1 vs B2 ===\n\n")
  s1 <- summarize_draws(subset_draws(B1$draws, variable = focal),
                        "mean", "sd", ~quantile(.x, c(0.05, 0.95)),
                        "rhat", "ess_bulk")
  s2 <- summarize_draws(subset_draws(B2$draws, variable = focal),
                        "mean", "sd", ~quantile(.x, c(0.05, 0.95)),
                        "rhat", "ess_bulk")
  cmp <- data.table(
    param   = s1$variable,
    B1_mean = round(s1$mean, 4),  B1_sd = round(s1$sd, 4),
    B1_q5   = round(s1$`5%`, 4),  B1_q95 = round(s1$`95%`, 4),
    B1_rhat = round(s1$rhat, 3),
    B2_mean = round(s2$mean, 4),  B2_sd = round(s2$sd, 4),
    B2_q5   = round(s2$`5%`, 4),  B2_q95 = round(s2$`95%`, 4),
    B2_rhat = round(s2$rhat, 3)
  )
  cmp[, mean_diff := round(B1_mean - B2_mean, 4)]
  cmp[, sd_ratio  := round(B1_sd / B2_sd, 3)]
  print(cmp, nrows = 50)

  fwrite(cmp, "~/fvs-modern/MEMORY/b1_b2_dg_kue_param_comparison.csv")
  cat("\nWrote ~/fvs-modern/MEMORY/b1_b2_dg_kue_param_comparison.csv\n")
}

## 2. sigma_sp recap from B2 -------------------------------------------------

cat("\n=== sigma_sp in B2 (the residual species effect B1 dropped) ===\n")
if ("sigma_sp" %in% vars_B2) {
  s_sp <- summarize_draws(subset_draws(B2$draws, variable = "sigma_sp"),
                          "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk")
  print(s_sp)
  cat(sprintf("\nB2 sigma_sp posterior mean = %.3f vs prior Normal+(0, 0.15)\n",
              s_sp$mean))
  cat("Per the HG B1 result, traits should explain most species variance.\n")
}

## 3. Sigma comparison: THE manuscript decision number -----------------------

cat("\n=== sigma posterior comparison (the manuscript decision) ===\n")
if ("sigma" %in% vars_B1 && "sigma" %in% vars_B2) {
  d1 <- as.numeric(subset_draws(B1$draws, variable = "sigma"))
  d2 <- as.numeric(subset_draws(B2$draws, variable = "sigma"))
  cat(sprintf("B1 sigma: mean = %.4f  sd = %.4f  q5..q95 = [%.4f, %.4f]\n",
              mean(d1), sd(d1), quantile(d1, 0.05), quantile(d1, 0.95)))
  cat(sprintf("B2 sigma: mean = %.4f  sd = %.4f  q5..q95 = [%.4f, %.4f]\n",
              mean(d2), sd(d2), quantile(d2, 0.05), quantile(d2, 0.95)))
  diff_pct <- 100 * (mean(d1) - mean(d2)) / mean(d2)
  cat(sprintf("B1 sigma is %+0.2f%% relative to B2.\n", diff_pct))
  if (abs(diff_pct) < 5) {
    cat("VERDICT: B1 ~= B2 (within 5%). Species-free works on DG_Kue.\n")
    cat("The trait-driven story is validated on a second base model.\n")
  } else if (diff_pct > 5) {
    cat("VERDICT: B1 > B2 by more than 5%. Species RE was carrying signal\n")
    cat("traits could not fully replace for DG_Kue. The species-free pivot\n")
    cat("works for HG but is base-model dependent.\n")
  } else {
    cat("VERDICT: B1 < B2. Surprising. Worth a sanity check.\n")
  }
}

## 4. gamma comparison (8 trait coefficients) --------------------------------

gammas <- grep("^gamma", shared, value = TRUE)
if (length(gammas) > 0) {
  cat("\n=== gamma (trait loading) comparison ===\n")
  cat(sprintf("Number of trait coefficients: %d\n", length(gammas)))
  s1g <- summarize_draws(subset_draws(B1$draws, variable = gammas),
                         "mean", "sd", ~quantile(.x, c(0.05, 0.95)))
  s2g <- summarize_draws(subset_draws(B2$draws, variable = gammas),
                         "mean", "sd", ~quantile(.x, c(0.05, 0.95)))
  gcmp <- data.table(
    trait_idx = s1g$variable,
    B1_mean = round(s1g$mean, 3),
    B1_ci   = sprintf("[%+0.3f, %+0.3f]", s1g$`5%`, s1g$`95%`),
    B2_mean = round(s2g$mean, 3),
    B2_ci   = sprintf("[%+0.3f, %+0.3f]", s2g$`5%`, s2g$`95%`),
    abs_mean_diff    = round(abs(s1g$mean - s2g$mean), 3),
    B1_excludes_zero = (s1g$`5%` * s1g$`95%`) > 0,
    B2_excludes_zero = (s2g$`5%` * s2g$`95%`) > 0
  )
  print(gcmp, nrows = 50)

  cat(sprintf("\nB1 gammas with 90%% CI excluding zero: %d / %d\n",
              sum(gcmp$B1_excludes_zero), nrow(gcmp)))
  cat(sprintf("B2 gammas with 90%% CI excluding zero: %d / %d\n",
              sum(gcmp$B2_excludes_zero), nrow(gcmp)))

  fwrite(gcmp, "~/fvs-modern/MEMORY/b1_b2_dg_kue_gamma_comparison.csv")
  cat("\nWrote ~/fvs-modern/MEMORY/b1_b2_dg_kue_gamma_comparison.csv\n")
}

cat("\n=== Done ===\n")
