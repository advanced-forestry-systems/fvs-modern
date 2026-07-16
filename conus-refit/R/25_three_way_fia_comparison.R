#!/usr/bin/env Rscript
################################################################################
# 25_three_way_fia_comparison.R
# ------------------------------------------------------------------------------
# Build a head-to-head comparison table for FVS-ACD default vs FVS-ACD Bayesian-
# calibrated vs OSM-ACD on the same FIA plot pair universe.
#
# Reads:
#   benchmark_osm_summary.csv         (produced by 24_osm_fia_benchmark.R)
#   fia_benchmark_results.csv         (produced by 19_fia_benchmark_engine.R)
#
# Writes:
#   three_way_acd_comparison.csv      long format: model x variable x stat
#   three_way_acd_summary.csv         wide table: one row per variable
#
# Units: metric throughout. FVS results in fia_benchmark_results.csv are reported
# in mixed metric/imperial; we extract the metric columns where available and
# convert otherwise.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--osm-summary",  type = "character",
              help = "benchmark_osm_summary.csv from 24_osm_fia_benchmark.R"),
  make_option("--fvs-benchmark", type = "character",
              help = "fia_benchmark_results.csv from 19_fia_benchmark_engine.R"),
  make_option("--variant",      type = "character", default = "ACD"),
  make_option("--output-dir",   type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt[["output-dir"]], showWarnings = FALSE, recursive = TRUE)

osm <- fread(opt[["osm-summary"]])
fvs <- fread(opt[["fvs-benchmark"]])

acd <- fvs[toupper(VARIANT) == opt$variant]
if (nrow(acd) == 0) {
  stop(sprintf("No ACD row in %s", opt[["fvs-benchmark"]]))
}

# --- Build long format from FVS benchmark wide row ---
# fia_benchmark_results.csv reports each variable with _calib and _default suffix
# (e.g., BA_RMSE_calib, BA_RMSE_default). It is in imperial native units (ft^2/ac).
FTAC2_TO_M2HA <- 0.2295684
ACRES_PER_HA  <- 2.4710538147
FT_TO_M       <- 0.3048

# Helper to grab a stat for a variable from a specific suffix
get_stat <- function(var, stat, suffix) {
  col <- paste0(var, "_", stat, "_", suffix)
  if (col %in% names(acd)) acd[[col]] else NA_real_
}

# FVS variables we have direct stats for, with unit conversion factor
# (FVS native -> metric reporting unit used in OSM table)
fvs_var_specs <- list(
  list(fvs = "BA",  out_var = "BA_m2ha",  conv = FTAC2_TO_M2HA),
  list(fvs = "TPA", out_var = "TPH",      conv = ACRES_PER_HA),
  list(fvs = "QMD", out_var = "QMD_cm",   conv = 2.54),   # inches -> cm
  list(fvs = "HT_top", out_var = "topht_m", conv = FT_TO_M)
)

fvs_rows <- list()
for (spec in fvs_var_specs) {
  for (suf in c("calib", "default")) {
    n <- get_stat(spec$fvs, "obs_mean", suf)
    if (is.na(n)) next
    obs_mean   <- get_stat(spec$fvs, "obs_mean", suf) * spec$conv
    pred_mean  <- get_stat(spec$fvs, "pred_mean", suf) * spec$conv
    bias       <- get_stat(spec$fvs, "bias", suf) * spec$conv
    bias_pct   <- get_stat(spec$fvs, "bias_pct", suf)
    rmse       <- get_stat(spec$fvs, "RMSE", suf) * spec$conv
    rmse_pct   <- get_stat(spec$fvs, "RMSE_pct", suf)
    mae        <- get_stat(spec$fvs, "MAE", suf) * spec$conv
    r2         <- get_stat(spec$fvs, "r2", suf)
    equiv      <- get_stat(spec$fvs, "equiv", suf)
    model_lab  <- if (suf == "calib") "FVS-ACD calibrated" else "FVS-ACD default"
    fvs_rows[[length(fvs_rows) + 1]] <- data.table(
      model = model_lab, variant = "ACD",
      variable = spec$out_var,
      n = acd$n_conditions,
      obs_mean = obs_mean, pred_mean = pred_mean,
      bias = bias, bias_pct = bias_pct,
      rmse = rmse, rmse_pct = rmse_pct,
      mae = mae, r2 = r2, equiv_pct = equiv
    )
  }
}
fvs_long <- rbindlist(fvs_rows, use.names = TRUE)

# Combine with OSM
all_long <- rbind(fvs_long, osm, use.names = TRUE, fill = TRUE)
all_long[, variant := "ACD"]
all_long[, model := factor(model, levels = c(
  "FVS-ACD default", "FVS-ACD calibrated", "OSM-ACD"))]
setorder(all_long, variable, model)

fwrite(all_long, file.path(opt[["output-dir"]], "three_way_acd_comparison.csv"))

# Wide summary: rows = variable, cols = model x stat
wide <- dcast(all_long,
              variable + n + obs_mean ~ model,
              value.var = c("pred_mean", "bias_pct", "rmse_pct", "r2", "equiv_pct"))
fwrite(wide, file.path(opt[["output-dir"]], "three_way_acd_summary.csv"))

cat("=== Three-way ACD comparison ===\n")
cat("(FVS-modern obs_mean from imperial -> metric conversion; OSM obs_mean from direct metric aggregation. Numbers may not match exactly across models because they were computed on slightly different plot subsets and slightly different observed-stand aggregation rules. Use the percent statistics for cross-model comparison.)\n\n")
print(all_long[, .(model, variable, n, obs_mean = round(obs_mean, 2),
                    pred_mean = round(pred_mean, 2),
                    bias_pct = round(bias_pct, 2),
                    rmse_pct = round(rmse_pct, 2),
                    r2 = round(r2, 3),
                    equiv_pct = round(equiv_pct, 1))])
cat(sprintf("\nWrote %s\n", file.path(opt[["output-dir"]], "three_way_acd_comparison.csv")))
cat(sprintf("Wrote %s\n", file.path(opt[["output-dir"]], "three_way_acd_summary.csv")))
