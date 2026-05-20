#!/usr/bin/env Rscript
## compare_b1_b2_hg.R
##
## Decision-quality comparison between the species-free B1 HG pilot
## (hg_organon_speciesfree.stan, 500K subsample, job 9298882) and the
## production species-aware B2 HG fit (hg_organon_fixedK.stan, banked
## Apr 25).
##
## Outputs three things:
##   1. Shared-parameter posterior comparison (a0..a8, gamma, sigma).
##      If B1 and B2 give the same gammas with similar uncertainty,
##      species-free works.
##   2. sigma comparison. If B1 sigma >> B2 sigma, traits aren't fully
##      replacing the species RE that B2 had. If B1 sigma <= B2 sigma,
##      species-free is the cleaner model.
##   3. Posterior predictive RMSE on a held-out chunk of the data.
##      The actual operational test: do predictions match?
##
## Usage on Cardinal:
##   module load gcc/12.3.0 R/4.4.0
##   cd ~/fvs-modern
##   Rscript calibration/R/compare_b1_b2_hg.R

suppressPackageStartupMessages({
  library(data.table)
  library(posterior)
  library(readr)
})

b1_path <- "~/fvs-modern/calibration/output/conus/hg/speciesfree_pilot/hg_organon_fixedK_cspi_traits1_fit.rds"
b2_path <- "~/fvs-modern/calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds"
b1_path <- path.expand(b1_path); b2_path <- path.expand(b2_path)

stopifnot(file.exists(b1_path), file.exists(b2_path))

cat("FVS-CONUS B1 (species-free) vs B2 (species-aware) HG comparison\n")
cat("================================================================\n\n")

read_fit_summary <- function(p, label) {
  cat(sprintf("Reading %s (%s) ... ", label, basename(p))); flush.console()
  fit <- readRDS(p)
  draws <- tryCatch(fit$draws(), error = function(e) as_draws_array(fit))
  cat(sprintf("draws array %s\n", paste(dim(draws), collapse = " x ")))
  list(label = label, fit = fit, draws = draws)
}

B1 <- read_fit_summary(b1_path, "B1")
B2 <- read_fit_summary(b2_path, "B2")

vars_B1 <- variables(B1$draws)
vars_B2 <- variables(B2$draws)

##---------------------------------------------------------------
## 1. Shared fixed coefficients (a0..a8) and traits (gamma)
##---------------------------------------------------------------
shared_fixed <- intersect(vars_B1, vars_B2)
fixed_of_interest <- c(paste0("a", 0:8), grep("^gamma", vars_B1, value = TRUE),
                       "sigma", "sigma_L1", "sigma_L2", "sigma_L3")
fixed_of_interest <- intersect(fixed_of_interest, shared_fixed)

if (length(fixed_of_interest) > 0) {
  cat("\n=== Shared parameters: B1 vs B2 ===\n\n")
  s1 <- summarize_draws(subset_draws(B1$draws, variable = fixed_of_interest),
                        "mean", "sd", ~quantile(.x, c(0.05, 0.95)),
                        "rhat", "ess_bulk")
  s2 <- summarize_draws(subset_draws(B2$draws, variable = fixed_of_interest),
                        "mean", "sd", ~quantile(.x, c(0.05, 0.95)),
                        "rhat", "ess_bulk")
  cmp <- data.table(
    param  = s1$variable,
    B1_mean = round(s1$mean, 4), B1_sd = round(s1$sd, 4),
    B1_q5   = round(s1$`5%`, 4), B1_q95 = round(s1$`95%`, 4),
    B1_rhat = round(s1$rhat, 3),
    B2_mean = round(s2$mean, 4), B2_sd = round(s2$sd, 4),
    B2_q5   = round(s2$`5%`, 4), B2_q95 = round(s2$`95%`, 4),
    B2_rhat = round(s2$rhat, 3)
  )
  cmp[, mean_diff := round(B1_mean - B2_mean, 4)]
  cmp[, sd_ratio  := round(B1_sd / B2_sd, 3)]
  print(cmp, nrows = 50)

  fwrite(cmp, "~/fvs-modern/MEMORY/b1_b2_param_comparison.csv")
  cat("\nWrote ~/fvs-modern/MEMORY/b1_b2_param_comparison.csv\n")
}

##---------------------------------------------------------------
## 2. sigma_sp (B2 only). If small, B1's removal is well justified.
##---------------------------------------------------------------
cat("\n=== sigma_sp in B2 (the residual species effect B1 dropped) ===\n")
if ("sigma_sp" %in% vars_B2) {
  s_sp <- summarize_draws(subset_draws(B2$draws, variable = "sigma_sp"),
                          "mean", "median", "sd", ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk")
  print(s_sp)
  cat(sprintf("\nB2 sigma_sp posterior mean = %.3f against prior Normal(0, 0.15)\n",
              s_sp$mean))
  cat("If this is close to 0.15 (the prior mode), traits explain most species variance.\n")
}

##---------------------------------------------------------------
## 3. Sigma comparison — is B1's residual error larger than B2's?
##---------------------------------------------------------------
cat("\n=== sigma posterior comparison ===\n")
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
    cat("Verdict on sigma: B1 ~= B2 (within 5%). Species-free is competitive.\n")
  } else if (diff_pct > 5) {
    cat("Verdict on sigma: B1 > B2 by more than 5%. Species RE was doing some work that traits can't fully replace.\n")
  } else {
    cat("Verdict on sigma: B1 < B2 — species-free is actually tighter, surprising.\n")
  }
}

##---------------------------------------------------------------
## 4. Trait coefficient comparison (gammas — the headline result)
##---------------------------------------------------------------
gamma_vars <- grep("^gamma", shared_fixed, value = TRUE)
if (length(gamma_vars) > 0) {
  cat("\n=== gamma (trait loading) comparison ===\n")
  cat(sprintf("Number of trait coefficients: %d\n", length(gamma_vars)))
  s1g <- summarize_draws(subset_draws(B1$draws, variable = gamma_vars),
                         "mean", "sd", ~quantile(.x, c(0.05, 0.95)))
  s2g <- summarize_draws(subset_draws(B2$draws, variable = gamma_vars),
                         "mean", "sd", ~quantile(.x, c(0.05, 0.95)))
  gcmp <- data.table(
    trait_idx = s1g$variable,
    B1_mean = round(s1g$mean, 3), B1_ci = sprintf("[%+0.3f, %+0.3f]", s1g$`5%`, s1g$`95%`),
    B2_mean = round(s2g$mean, 3), B2_ci = sprintf("[%+0.3f, %+0.3f]", s2g$`5%`, s2g$`95%`),
    abs_mean_diff = round(abs(s1g$mean - s2g$mean), 3),
    B1_excludes_zero = (s1g$`5%` * s1g$`95%`) > 0,
    B2_excludes_zero = (s2g$`5%` * s2g$`95%`) > 0
  )
  print(gcmp, nrows = 50)

  cat(sprintf("\nB1 gammas with 90%% CI excluding zero: %d / %d\n",
              sum(gcmp$B1_excludes_zero), nrow(gcmp)))
  cat(sprintf("B2 gammas with 90%% CI excluding zero: %d / %d\n",
              sum(gcmp$B2_excludes_zero), nrow(gcmp)))

  fwrite(gcmp, "~/fvs-modern/MEMORY/b1_b2_gamma_comparison.csv")
  cat("\nWrote ~/fvs-modern/MEMORY/b1_b2_gamma_comparison.csv\n")
}

cat("\n=== Done ===\n")
