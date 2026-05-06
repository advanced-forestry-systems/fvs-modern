## =============================================================================
## 20260422_fvs-conus_convergence_diagnostics.R
##
## Post-fit convergence diagnostic for the 7 production Stan fits on Cardinal:
##   HT-DBH Wykoff log-normal (8682996)
##   DG ORGANON             (8682998)
##   DG Kuehne              (8682999)
##   HG ORGANON fixedK      (8694614)
##   HCB ORGANON            (8694618)
##   CR Hann and Hanus      (8694619)
##   Mortality logit        (8694620)
##
## Reads the *_meta.rds and optional *_fit.rds written by each driver,
## summarises rhat / ESS / divergence / wall time, and prints posterior
## medians and 95 percent CIs for the top-level fixed effects. Flags any
## parameter with rhat >= 1.05 or ESS_bulk < 400. Drops RE vectors since
## per-species / per-ecoregion rhat would overwhelm the report.
##
## Run on Cardinal:
##   module load gcc/12.3.0 R/4.4.0
##   cd ~/fvs-modern
##   Rscript calibration/R/diagnostics_20260422.R
##
## Author: A. Weiskittel, 2026-04-22
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
})

## Component configuration -----------------------------------------------------
components <- tribble(
  ~name,        ~dir,                                         ~meta_pattern,                 ~fit_pattern,
  "HT_DBH",     "calibration/output/conus/htdbh",             "htdbh_wykoff_lognormal_.*_meta\\.rds$", "htdbh_wykoff_lognormal_.*_fit\\.rds$",
  "DG_organon", "calibration/output/conus/dg",                "dg_organon_.*_meta\\.rds$",    "dg_organon_.*_fit\\.rds$",
  "DG_kuehne",  "calibration/output/conus/dg",                "dg_kuehne_.*_meta\\.rds$",     "dg_kuehne_.*_fit\\.rds$",
  "HG",         "calibration/output/conus/hg",                "hg_organon_fixedK_.*_meta\\.rds$", "hg_organon_fixedK_.*_fit\\.rds$",
  "HCB",        "calibration/output/conus/hcb",               "hcb_.*_meta\\.rds$",            "hcb_.*_fit\\.rds$",
  "CR",         "calibration/output/conus/crown_recession",   "cr_.*_meta\\.rds$",             "cr_.*_fit\\.rds$",
  "Mortality",  "calibration/output/conus/mortality",         "mort_.*_meta\\.rds$",           "mort_.*_fit\\.rds$"
)

## Helpers --------------------------------------------------------------------
latest_match <- function(dir, pattern) {
  if (!dir.exists(dir)) return(NA_character_)
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0L) return(NA_character_)
  f[order(file.info(f)$mtime, decreasing = TRUE)][1]
}

summarise_one <- function(row) {
  cat("\n==============================================\n")
  cat(sprintf("Component: %s\n", row$name))
  cat("==============================================\n")

  meta_path <- latest_match(row$dir, row$meta_pattern)
  fit_path  <- latest_match(row$dir, row$fit_pattern)

  if (is.na(meta_path)) {
    cat("  [not yet completed: meta rds missing]\n")
    return(invisible(NULL))
  }

  meta <- readRDS(meta_path)
  cat("  Meta:     ", basename(meta_path), "\n")
  cat("  Stan file:", basename(meta$stan_file %||% "?"), "\n")
  cat("  N obs:    ", format(meta$n_obs  %||% NA, big.mark = ","), "\n")
  cat("  Wall (min):", round(meta$wall_min %||% NA, 1), "\n")
  if (!is.null(meta$trait_cols))
    cat("  Traits:   ", paste(meta$trait_cols, collapse = ", "), "\n")

  ## From summary in meta (fixed effects only)
  summ <- meta$summary %>%
    mutate(
      rhat_flag = ifelse(rhat    >= 1.05, "!", ""),
      ess_flag  = ifelse(ess_bulk < 400,  "!", "")
    )

  cat("\nFixed-effect / variance-component summary (rhat, ess):\n")
  summ %>%
    select(variable, mean, median, q5, q95, rhat, ess_bulk, rhat_flag, ess_flag) %>%
    mutate(across(where(is.numeric), ~round(., 3))) %>%
    print(n = Inf)

  bad_rhat <- summ %>% filter(rhat >= 1.05)
  bad_ess  <- summ %>% filter(ess_bulk < 400)
  cat(sprintf("\n  Flagged rhat >= 1.05: %d\n", nrow(bad_rhat)))
  cat(sprintf("  Flagged ess_bulk < 400: %d\n", nrow(bad_ess)))

  ## Full fit object: divergences + per-species RE rhat percentiles
  if (!is.na(fit_path) && file.exists(fit_path)) {
    cat("\n  Fit:     ", basename(fit_path), "\n")
    fit <- readRDS(fit_path)
    diag <- fit$diagnostic_summary()
    cat(sprintf("  Divergences (per chain): %s\n",
                paste(diag$num_divergent, collapse = ", ")))
    cat(sprintf("  Max treedepth hits: %s\n",
                paste(diag$num_max_treedepth, collapse = ", ")))
    cat(sprintf("  E-BFMI per chain: %s\n",
                paste(round(diag$ebfmi, 2), collapse = ", ")))

    ## Species and ecoregion RE rhat distribution
    re_vars <- grep("^z_(sp|L1|L2|L3)(_raw)?\\[", fit$metadata()$variables,
                    value = TRUE)
    if (length(re_vars) > 0) {
      re_summ <- fit$summary(variables = re_vars)
      cat("\n  Random-effect rhat quantiles (sp / L1 / L2 / L3):\n")
      groups <- c("z_sp", "z_L1", "z_L2", "z_L3")
      for (g in groups) {
        idx <- grepl(sprintf("^%s(_raw)?\\[", g), re_summ$variable)
        if (any(idx)) {
          q <- quantile(re_summ$rhat[idx],
                        c(0.5, 0.9, 0.99, 1.0), na.rm = TRUE)
          cat(sprintf("    %-6s: median=%.3f  p90=%.3f  p99=%.3f  max=%.3f  (n=%d)\n",
                      g, q[1], q[2], q[3], q[4], sum(idx)))
        }
      }
    }
  } else {
    cat("\n  [fit rds not present: run with --save_draws to keep full draws]\n")
  }

  invisible(list(meta = meta, bad_rhat = bad_rhat, bad_ess = bad_ess))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("==============================================================\n")
cat("FVS-CONUS Production Convergence Diagnostics\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("==============================================================\n")

results <- vector("list", nrow(components))
for (i in seq_len(nrow(components))) {
  results[[i]] <- summarise_one(components[i, ])
}

## Final roll-up
cat("\n==============================================================\n")
cat("Overall roll-up\n")
cat("==============================================================\n")
roll <- map2_dfr(results, components$name, function(r, nm) {
  if (is.null(r)) return(tibble(component = nm, status = "missing"))
  tibble(
    component       = nm,
    status          = "completed",
    bad_rhat_count  = nrow(r$bad_rhat),
    bad_ess_count   = nrow(r$bad_ess),
    wall_min        = r$meta$wall_min %||% NA_real_,
    n_obs           = r$meta$n_obs    %||% NA_integer_
  )
})
print(roll, n = Inf)

cat("\nDone.\n")
