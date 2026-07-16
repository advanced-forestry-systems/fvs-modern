#!/usr/bin/env Rscript
# elpd_compare_ingrowth.R
# Settle the ingrowth base form (hurdle vs negative binomial) by held-out ELPD via PSIS-LOO.
# Run after the hurdle and negbinom fits complete. Both Stan models emit log_lik, so loo is
# computed from the saved CmdStanMCMC objects. Requires the two fits to cover the same plot
# observations in the same order; the script checks N and warns if they differ.
#
# Usage:
#   Rscript R/elpd_compare_ingrowth.R <hurdle_fit.rds> <negbinom_fit.rds> [out_dir]
suppressPackageStartupMessages({ library(loo); library(posterior) })
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("need: <hurdle_fit.rds> <negbinom_fit.rds> [out_dir]")
hurdle_f <- a[1]; negbin_f <- a[2]
out_dir  <- if (length(a) >= 3) a[3] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

get_loo <- function(path) {
  cat("loading", path, "\n"); fit <- readRDS(path)
  ll <- tryCatch(fit$draws("log_lik"), error = function(e) NULL)
  if (is.null(ll)) stop("no log_lik draws in ", path)
  ll <- posterior::as_draws_matrix(ll)             # iterations x observations
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1:posterior::nchains(fit$draws()),
                                                     each = nrow(ll) / posterior::nchains(fit$draws())))
  loo::loo(ll, r_eff = r_eff)
}

loo_h <- get_loo(hurdle_f)
loo_n <- get_loo(negbin_f)
cat("\n-- hurdle --\n");   print(loo_h$estimates)
cat("\n-- negbinom --\n"); print(loo_n$estimates)

nh <- nrow(loo_h$pointwise); nn <- nrow(loo_n$pointwise)
if (nh != nn) {
  cat(sprintf("\nWARNING: observation counts differ (hurdle %d, negbinom %d). ", nh, nn))
  cat("loo_compare requires identical observations; refit both on one frozen plot set before trusting the comparison.\n")
} else {
  cmp <- loo::loo_compare(list(hurdle = loo_h, negbinom = loo_n))
  cat("\n=== loo_compare (positive elpd_diff favors the first row) ===\n")
  print(cmp)
  saveRDS(list(loo_hurdle = loo_h, loo_negbinom = loo_n, compare = cmp),
          file.path(out_dir, "ingrowth_elpd_compare.rds"))
  utils::write.csv(as.data.frame(cmp), file.path(out_dir, "ingrowth_elpd_compare.csv"))
  best <- rownames(cmp)[1]
  cat(sprintf("\nVerdict: %s has the higher ELPD (better held-out predictive fit).\n", best))
  cat("Decision rule: if elpd_diff exceeds about 2x its SE, the difference is meaningful;\n")
  cat("otherwise prefer the simpler negbinom base.\n")
}
cat("DONE_ELPD\n")
