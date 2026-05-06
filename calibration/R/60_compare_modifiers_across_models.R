## =============================================================================
## 60_compare_modifiers_across_models.R
##
## Cross-model modifier comparison. Walks the
## calibration/output/conus/<model>/modifier/ directories, loads every
## *_summary.csv it finds, and produces a tidy side-by-side table of
## the modifier alphas with their 90% CIs and convergence diagnostics.
##
## Output: tidy CSV plus a quick text table for deck inclusion.
##
## Usage:
##   Rscript calibration/R/60_compare_modifiers_across_models.R \
##       --base_dir calibration/output/conus \
##       --out      calibration/output/conus/modifier_comparison.csv
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(readr)
})

opts <- OptionParser(option_list = list(
  make_option("--base_dir", type = "character",
              default = "calibration/output/conus"),
  make_option("--out",      type = "character",
              default = "calibration/output/conus/modifier_comparison.csv")
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS modifier comparison across base models\n")
cat("  base_dir:", opts$base_dir, "\n")
cat("  out     :", opts$out, "\n")
cat("==============================================================\n")

# Find all modifier summary CSVs
candidate_dirs <- list.dirs(opts$base_dir, full.names = TRUE, recursive = TRUE)
modifier_dirs  <- candidate_dirs[grepl("modifier($|/|_v2|_traitmed)", candidate_dirs)]
cat("Modifier directories found:\n")
for (d in modifier_dirs) cat("  ", d, "\n")

key_alphas <- c("alpha_0", "alpha_plant",
                "alpha_fire", "alpha_wind", "alpha_harvest",
                "alpha_insect", "alpha_disease",
                "alpha_cutting", "alpha_siteprep",
                "sigma_L1", "sigma_resid")

rows <- list()
for (d in modifier_dirs) {
  csvs <- list.files(d, pattern = "summary\\.csv$", full.names = TRUE)
  for (cs in csvs) {
    df <- tryCatch(read_csv(cs, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    setDT(df)
    if (!"variable" %in% names(df)) next

    # Derive model + variant from path
    model_name <- basename(dirname(d))
    if (model_name == "modifier" || model_name == "modifier_v2" ||
        model_name == "modifier_traitmed") {
      model_name <- basename(dirname(dirname(d)))
    }
    variant <- basename(d)
    fit_tag <- sub("_summary\\.csv$", "", basename(cs))

    df[, model   := model_name]
    df[, variant := variant]
    df[, fit_tag := fit_tag]

    keep <- df[variable %in% key_alphas]
    if (nrow(keep) > 0) rows[[length(rows) + 1L]] <- keep
  }
}

if (length(rows) == 0) {
  stop("No modifier summary CSVs found anywhere under ", opts$base_dir)
}

all <- rbindlist(rows, fill = TRUE)
all[, q90_lo := q5]
all[, q90_hi := q95]

# Pretty wide table: model x variant x fit_tag rows, alphas as columns
out_long <- all[, .(model, variant, fit_tag, variable, mean, q90_lo, q90_hi,
                    rhat, ess_bulk)]
setorder(out_long, model, variant, variable)

write_csv(out_long, opts$out)
cat(sprintf("\nWrote %d rows to %s\n", nrow(out_long), opts$out))

cat("\n=== Headline summary table (mean alpha, 90%% CI) ===\n")
for (m in unique(out_long$model)) {
  cat(sprintf("\n--- %s ---\n", m))
  sub <- out_long[model == m]
  for (v in unique(sub$variant)) {
    sv <- sub[variant == v]
    cat(sprintf("  variant: %s\n", v))
    for (i in seq_len(nrow(sv))) {
      r <- sv[i]
      cat(sprintf("    %-18s mean=%7.4f  CI=[%7.4f, %7.4f]  rhat=%.3f  ess=%4.0f\n",
                  r$variable, r$mean, r$q90_lo, r$q90_hi, r$rhat, r$ess_bulk))
    }
  }
}

cat("\nDone.\n")
