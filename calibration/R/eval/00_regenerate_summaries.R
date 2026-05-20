# =============================================================================
# Title: Regenerate base-model summary CSVs in a consistent format
# Author: A. Weiskittel
# Date: 2026-05-12
# Description: The banked fits from April have summary CSVs with different
#              column conventions (some have q5/q95, others 2.5%/97.5%, some
#              have median, some don't). This script loads each banked fit
#              and writes a unified _eval_summary.csv at:
#                  output/conus/evaluation_summaries/{model}.csv
#              with columns: variable, mean, median, sd, q5, q95, rhat,
#              ess_bulk, ess_tail
#
#              The convergence dashboard and trait coefficient scripts read
#              from these unified files.
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_00_regen.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/00_regenerate_summaries.R
# =============================================================================

library(tidyverse)
library(posterior)

PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR   <- file.path(PROJ_ROOT, "calibration/output/evaluation/summaries")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# All banked production fits
fit_specs <- tribble(
  ~name,            ~path,
  "hg_b2",          "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds",
  "hg_b1",          "calibration/output/conus/hg/speciesfree_pilot/hg_organon_fixedK_cspi_traits1_fit.rds",
  "dg_kue_b2",      "calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds",
  "dg_kue_b1",      "calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds",
  "htdbh_b2",       "calibration/output/conus/ht_dbh/htdbh_wykoff_lognormal_cspi_traits1_fit.rds",
  "cr_b2",          "calibration/output/conus/crown_recession/cr_recession_cspi_traits1_fit.rds",
  "hcb_b2",         "calibration/output/conus/hcb/hcb_organon_cspi_traits1_fit.rds",
  "mort_b2",        "calibration/output/conus/mortality/mort_logit_simple_cspi_traits1_fit.rds"
)

# For each fit, extract summary on fixed-effect parameters only (a/b coefs,
# gammas, sigma_*, sigma, plus trait_effect for species)
extract_summary <- function(name, path) {
  full <- file.path(PROJ_ROOT, path)
  if (!file.exists(full)) {
    message(sprintf("  SKIP %s: file not found", name))
    return(invisible(NULL))
  }
  message(sprintf("Processing %s (%.1f GB) ...", name,
                  file.info(full)$size / 1e9))
  flush.console()
  fit <- readRDS(full)

  # Variable groups to summarize. Avoid the gigantic random-effect deviations
  # (z_sp, z_L1, z_L2, z_L3) which add billions of cells. Keep their scales.
  vars_to_keep <- function() {
    all_vars <- variables(fit$draws())
    keep <- c(
      grep("^a[0-9]$", all_vars, value = TRUE),
      grep("^b[0-9]$", all_vars, value = TRUE),
      grep("^gamma\\[", all_vars, value = TRUE),
      grep("^sigma$|^sigma_(sp|L1|L2|L3)$|^phi$", all_vars, value = TRUE),
      grep("^trait_effect\\[", all_vars, value = TRUE)
    )
    unique(keep)
  }

  vars <- vars_to_keep()
  message(sprintf("  %d parameters", length(vars)))

  summary_df <- fit$summary(
    variables = vars,
    "mean", "median", "sd",
    ~quantile(.x, c(0.05, 0.95), na.rm = TRUE),
    "rhat", "ess_bulk", "ess_tail"
  ) %>%
    rename(q5 = `5%`, q95 = `95%`)

  out_path <- file.path(OUT_DIR, paste0(name, "_eval_summary.csv"))
  write_csv(summary_df, out_path)
  message(sprintf("  wrote %s", out_path))

  rm(fit, summary_df); gc()
}

for (i in seq_len(nrow(fit_specs))) {
  extract_summary(fit_specs$name[i], fit_specs$path[i])
}

message("\nAll done.")
