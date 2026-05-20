#!/usr/bin/env Rscript
## reconstruct_modifier_summaries.R
##
## Each modifier meta.rds bundle contains a $summary slot with the full
## per-parameter draws_summary (mean, median, sd, q5, q95, rhat, ess_*).
## This script reads each meta.rds in the modifier dirs and writes the
## corresponding *_summary.csv that 60b_pivot_modifier_table.R expects.
## Also copies Kuehne's modifier_kuehne files into the standard modifier
## path so the pivot's regex picks them up.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
})

base <- "~/fvs-modern/calibration/output/conus"
base <- path.expand(base)

# Reconstruct summary CSV for each meta.rds in any modifier* directory
dirs <- list.dirs(base, full.names = TRUE, recursive = TRUE)
dirs <- dirs[grepl("/modifier($|_)", dirs)]

n_written <- 0
for (d in dirs) {
  meta_files <- list.files(d, pattern = "_meta\\.rds$", full.names = TRUE)
  for (mf in meta_files) {
    csv_path <- sub("_meta\\.rds$", "_summary.csv", mf)
    meta <- tryCatch(readRDS(mf), error = function(e) NULL)
    if (is.null(meta) || is.null(meta$summary)) next
    summary_tbl <- as.data.table(meta$summary)
    write_csv(summary_tbl, csv_path)
    n_written <- n_written + 1
    cat(sprintf("  wrote %s\n", csv_path))
  }
}
cat(sprintf("\nReconstructed %d summary CSVs\n\n", n_written))

# For Kuehne: copy from modifier_kuehne to a parallel modifier_kuehne_pivot
# that the regex /modifier($|/)/ will match. Use the path
# .../dg/modifier_kuehne/<files> by also linking into .../dg/modifier/.
# Easier: dump kuehne files into .../dg/modifier/ with a kuehne_ prefix
# so they don't collide with anything (currently empty).
kue_src <- file.path(base, "dg/modifier_kuehne")
kue_dst <- file.path(base, "dg/modifier")
if (dir.exists(kue_src) && dir.exists(kue_dst)) {
  for (f in list.files(kue_src, pattern = "_summary\\.csv$", full.names = TRUE)) {
    dst <- file.path(kue_dst, basename(f))
    file.copy(f, dst, overwrite = TRUE)
    cat(sprintf("  copied %s -> %s\n", f, dst))
  }
}

cat("\nDone.\n")
