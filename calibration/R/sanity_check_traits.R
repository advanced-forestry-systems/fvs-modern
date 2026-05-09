#!/usr/bin/env Rscript
## sanity_check_traits.R
##
## Decision-quality summary for the species-free model design choice.
## For each banked base model on Cardinal, extracts:
##   - sigma_sp, sigma_eco, sigma (with rhat)
##   - all gamma_trait* coefficients (loading on each species trait)
##   - sigma_sp / sigma_eco ratio (the "is species RE doing more than ecoregion?" test)
##
## Run on Cardinal via:
##   module load gcc/12.3.0 R/4.4.0
##   Rscript ~/sanity_check_traits.R
##
## Output: prints a tidy summary to stdout, writes a CSV to
## ~/fvs-modern/MEMORY/sanity_check_traits_20260509.csv

suppressPackageStartupMessages({
  library(posterior)
  library(data.table)
})

base_models <- list(
  HG       = "~/fvs-modern/calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds",
  HT_DBH   = "~/fvs-modern/calibration/output/conus/ht_dbh/htdbh_wykoff_lognormal_cspi_traits1_fit.rds",
  CR       = "~/fvs-modern/calibration/output/conus/crown_recession/cr_recession_cspi_traits1_fit.rds",
  Mort     = "~/fvs-modern/calibration/output/conus/mortality/mort_logit_simple_cspi_traits1_fit.rds",
  HCB      = "~/fvs-modern/calibration/output/conus/hcb/hcb_organon_cspi_traits1_fit.rds",
  DG_Kue   = "~/fvs-modern/calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds"
)

cat("FVS-CONUS species-free design sanity check\n")
cat("==========================================\n\n")

results <- list()

for (mod in names(base_models)) {
  path <- base_models[[mod]]
  cat(sprintf("--- %s ---\n  fit: %s\n", mod, path))
  if (!file.exists(path)) {
    cat("  MISSING\n\n")
    next
  }
  size_gb <- round(file.info(path)$size / 1e9, 2)
  cat(sprintf("  size: %.2f GB ... reading\n", size_gb))

  fit <- readRDS(path)
  draws <- tryCatch(fit$draws(), error = function(e) NULL)
  if (is.null(draws)) {
    draws <- tryCatch(as_draws_array(fit), error = function(e) NULL)
  }
  if (is.null(draws)) {
    cat("  could not extract draws; skipping\n\n")
    next
  }

  vars <- variables(draws)
  pick <- function(p) intersect(p, vars)

  ## Pull the variance components
  want <- pick(c("sigma_sp", "sigma_eco", "sigma", "sigma_resid"))
  if (length(want) == 0) {
    cat("  no sigma_sp/sigma_eco/sigma found in posterior\n")
  } else {
    s <- summarize_draws(subset_draws(draws, variable = want),
                         "mean", "median", "sd", ~quantile(.x, c(0.05, 0.95)),
                         "rhat", "ess_bulk")
    print(s, n = 50)
    for (i in seq_len(nrow(s))) {
      results[[length(results) + 1L]] <- data.table(
        model = mod,
        param = s$variable[i],
        mean = s$mean[i],
        rhat = s$rhat[i],
        ess  = s$ess_bulk[i]
      )
    }
  }

  ## Pull all gamma_trait* if present
  trait_vars <- vars[grepl("^gamma_trait|^gamma_alpha|^trait_effect|^gamma_sp|^trait_load", vars)]
  if (length(trait_vars) > 0) {
    cat(sprintf("\n  TRAIT block (%d coefs):\n", length(trait_vars)))
    s2 <- summarize_draws(subset_draws(draws, variable = trait_vars),
                          "mean", "sd", ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk")
    s2 <- s2[order(abs(s2$mean), decreasing = TRUE), ]
    cat("  Top 10 by |mean|:\n")
    print(head(s2, 10), n = 10)
    cat(sprintf("  Coefs with 90%% CI excluding zero: %d / %d\n",
                sum(s2$`5%` * s2$`95%` > 0, na.rm = TRUE), nrow(s2)))
  } else {
    cat("\n  NO gamma_trait* / trait_effect block found\n")
  }
  cat("\n")
  rm(fit, draws); gc(verbose = FALSE)
}

if (length(results) > 0) {
  out <- rbindlist(results)
  cat("\n\n=== Variance components summary across models ===\n")
  print(dcast(out, model ~ param, value.var = "mean"), n = 50)

  ## sigma_sp / sigma_eco ratio
  wide <- dcast(out, model ~ param, value.var = "mean")
  if ("sigma_sp" %in% names(wide) && "sigma_eco" %in% names(wide)) {
    wide[, ratio_sp_eco := sigma_sp / sigma_eco]
    cat("\n=== sigma_sp / sigma_eco ratio ===\n")
    cat("(if <2, species RE adds little beyond ecoregion: B1 fully species-free should converge cleanly)\n")
    cat("(if >5, species RE absorbs a lot: B2 trait-anchored species hybrid is safer)\n\n")
    print(wide[, .(model, sigma_sp, sigma_eco, sigma, ratio_sp_eco)], n = 50)
  }

  out_path <- "~/fvs-modern/MEMORY"
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  fwrite(out, file.path(out_path, "sanity_check_traits_20260509.csv"))
  cat(sprintf("\nWrote %s/sanity_check_traits_20260509.csv\n", out_path))
}
