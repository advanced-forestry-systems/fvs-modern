#!/usr/bin/env Rscript
# Compare four ACD/NE benchmark results side-by-side:
#   1. baseline      — pre-refit, no post-pass (current cc00b92 result)
#   2. refit_only    — converged posterior, no post-pass
#   3. refit_postpass_pop  — + population multipliers
#   4. refit_postpass_strat_ny — + stratified post-pass + NY Adirondacks
#
# Reads the tagged CSVs under
# calibration/output/comparisons/manuscript_tables/ and writes a
# comparison artifact under
# calibration/analysis/acd_stand_level_2026-05-16/post_refit_comparison/

suppressMessages({
  library(data.table)
})

project_root <- "/users/PUOM0008/crsfaaron/fvs-modern-acdbridge"
src_root     <- file.path(project_root,
  "calibration/output/comparisons/manuscript_tables")
out_root     <- file.path(project_root,
  "calibration/analysis/acd_stand_level_2026-05-16/post_refit_comparison")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# Source filenames in priority order. We accept whichever exist (so the
# script is useful even before all four passes complete).
file_map <- c(
  baseline                = "fia_benchmark_pctrmse.csv",
  refit_only              = "fia_benchmark_pctrmse_refit_only.csv",
  refit_postpass_pop      = "fia_benchmark_pctrmse_refit_postpass_pop.csv",
  refit_postpass_strat_ny = "fia_benchmark_pctrmse_refit_postpass_strat_ny.csv"
)

read_one <- function(tag, fn) {
  p <- file.path(src_root, fn)
  if (!file.exists(p)) {
    cat(sprintf("[skip] %-26s  not found: %s\n", tag, fn))
    return(NULL)
  }
  d <- fread(p)
  # Keep ACD and NE rows; relabel ACD-flavored variants
  acd_re <- grepl("acd", d$VARIANT, ignore.case = TRUE)
  ne_re  <- d$VARIANT == "NE_acadian" | d$VARIANT == "NE"
  d <- d[acd_re | ne_re]
  if (nrow(d) == 0) {
    cat(sprintf("[skip] %-26s  no ACD/NE rows in %s\n", tag, fn))
    return(NULL)
  }
  d[, source := tag]
  d
}

rows <- rbindlist(Map(read_one, names(file_map), file_map), fill = TRUE)
if (nrow(rows) == 0) {
  cat("No usable input files yet. Re-run after the A/B chain completes.\n")
  quit(status = 0)
}

cols_keep <- intersect(
  c("source","VARIANT","n","mean_obs_BA","BA_pctRMSE_calib","BA_pctRMSE_default",
    "BA_pctBias_calib","BA_pctBias_default","pctRMSE_reduction",
    "VOL_pctRMSE_calib","VOL_pctRMSE_default","VOL_pctRMSE_reduction"),
  names(rows))
tab <- rows[, ..cols_keep]

cat("\n=== ACD vs NE side-by-side across configurations ===\n")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# Markdown summary
md_path <- file.path(out_root, "comparison.md")
md <- c(
  "# Post-refit ACD vs NE comparison",
  "",
  paste0("Source: tagged CSVs under `", src_root, "/`"),
  paste0("Configurations included: ", paste(unique(tab$source), collapse = ", ")),
  "",
  "## Side-by-side table",
  "",
  paste("|", paste(cols_keep, collapse = " | "), "|"),
  paste("|", paste(rep("---", length(cols_keep)), collapse = " | "), "|")
)
for (i in seq_len(nrow(tab))) {
  vals <- as.character(unlist(tab[i]))
  vals <- ifelse(is.na(vals), "-", vals)
  md <- c(md, paste("|", paste(vals, collapse = " | "), "|"))
}
md <- c(md, "",
  "## Headline numbers",
  "")

# Compute the headline ACD %RMSE delta vs baseline
acd_rows <- tab[grepl("acd|ACD", VARIANT, ignore.case = TRUE)]
if (nrow(acd_rows) > 0 && "baseline" %in% acd_rows$source) {
  base_pct <- acd_rows[source == "baseline", BA_pctRMSE_calib][1]
  for (s in setdiff(unique(acd_rows$source), "baseline")) {
    new_pct <- acd_rows[source == s, BA_pctRMSE_calib][1]
    delta <- new_pct - base_pct
    md <- c(md, sprintf("- ACD BA pctRMSE %s: %.2f (baseline %.2f, delta %+.2f pp)",
      s, new_pct, base_pct, delta))
  }
}

ne_rows <- tab[VARIANT %in% c("NE","NE_acadian")]
if (nrow(ne_rows) > 0) {
  md <- c(md, "",
    sprintf("- NE BA pctRMSE (baseline reference): %.2f",
      ne_rows[source == "baseline" | source == "refit_only",
              BA_pctRMSE_calib][1]))
}

writeLines(md, md_path)
cat(sprintf("\nWrote %s\n", md_path))
