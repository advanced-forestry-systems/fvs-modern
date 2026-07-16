#!/usr/bin/env Rscript
################################################################################
# 24_osm_fia_benchmark.R
# ------------------------------------------------------------------------------
# Parse OSM-ACD output for a batch of FIA plot pairs and compute the standard
# benchmark metric block (RMSE, bias, %bias, R^2, Willmott d, equivalence).
#
# Inputs:
#   --osm-stand      OSM StandListProjections.csv (one row per SurveyID x Period)
#   --osm-tree       OSM TreeListProjections.csv (one row per SurveyID x Period x tree)
#   --lookup         acd_plt_cn_lookup.csv produced by osm_input_builder.py
#   --observed-dg    fvs-modern diameter_growth.csv (the t1/t2 paired data)
#   --output-dir     Where to write benchmark_long.csv and figures
#
# Output: benchmark_osm_long.csv with columns
#   (variant, model, statistic, variable, value, n)
#
# Metric set (matched to 19_fia_benchmark_engine.R conventions):
#   variables: BA (m^2/ha), TPH, QMD (cm), top_ht (m), GTV (m^3/ha)
#   stats:     RMSE, MAE, bias, bias_pct, RMSE_pct, r2, equiv (10% threshold)
#
# Units: all metric. Output BA in m^2/ha, GTV in m^3/ha, etc.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--osm-stand",   type = "character", help = "OSM StandListProjections.csv"),
  make_option("--osm-tree",    type = "character", help = "OSM TreeListProjections.csv"),
  make_option("--lookup",      type = "character", help = "acd_plt_cn_lookup.csv"),
  make_option("--observed-dg", type = "character", help = "diameter_growth.csv with t1/t2"),
  make_option("--output-dir",  type = "character", help = "Output directory"),
  make_option("--variant",     type = "character", default = "acd"),
  make_option("--ypc",         type = "integer",   default = 5L)
)
opt <- parse_args(OptionParser(option_list = option_list))

INCH_TO_CM <- 2.54
FT_TO_M    <- 0.3048
FTAC2_TO_M2HA <- 0.2295684
ACRES_PER_HA <- 2.4710538147

cat("[24_osm_fia_benchmark] reading inputs (with pre-clean to drop OSM-beta corrupt lines)\n")

# Pre-clean OSM output: drop lines whose comma count doesn't match the header.
# OSM Linux beta v2.26.0 occasionally concatenates two records onto one line under heavy I/O.
clean_osm_csv <- function(path) {
  tmp <- tempfile(fileext = ".csv")
  cmd <- sprintf(
    "awk -F',' 'NR==1 {nc=NF; print; next} NF==nc {print}' %s > %s",
    shQuote(path), shQuote(tmp)
  )
  system(cmd)
  tmp
}
stand_clean <- clean_osm_csv(opt[["osm-stand"]])
tree_clean  <- clean_osm_csv(opt[["osm-tree"]])
osm_stand <- fread(stand_clean)
osm_tree  <- fread(tree_clean)
lookup    <- fread(opt$lookup)
observed  <- fread(opt[["observed-dg"]])
cat(sprintf("[24_osm_fia_benchmark] cleaned: %d stand rows, %d tree rows\n",
            nrow(osm_stand), nrow(osm_tree)))

# ---- Restrict observed to plots actually in this OSM batch ----
observed <- observed[PLT_CN %in% lookup$PLT_CN]
cat(sprintf("[24_osm_fia_benchmark] %d OSM plots, %d observed-tree rows\n",
            nrow(lookup), nrow(observed)))

# ---- Select OSM projection period matching each plot's years_interval ----
# Each plot's `years_interval` is in years; OSM YPC is 5; so Period = round(years_interval / YPC).
# (Use nearest cycle. For intervals not on a YPC multiple, this picks the nearer of the two cycles
#  available; for a tighter validation, interpolation could be added.)
lookup[, target_period := round(years_interval / opt$ypc)]
lookup[target_period < 1, target_period := 1L]

# Join lookup -> osm_stand to pick projected t2 stand metrics
setnames(osm_stand, c("SurveyID", "Period"), c("SurveyID", "Period"))
proj <- merge(
  lookup[, .(PLT_CN, SurveyID, years_interval, target_period, ba_t1_m2ha, tph_t1)],
  osm_stand[, .(SurveyID, Period, BA_proj = BA, Trees_proj = Trees,
                 QMD_proj = QMD, LryH_proj = LryH, T50HT_proj = T50HT,
                 GMV_proj = GMV)],
  by.x = c("SurveyID", "target_period"),
  by.y = c("SurveyID", "Period"),
  all.x = TRUE
)
cat(sprintf("[24_osm_fia_benchmark] joined %d plots with OSM projection rows\n", nrow(proj)))

# ---- Aggregate OSM TreeList for GTV (total volume) at the target period ----
# osm_tree has columns including GTV (m^3 per tree) and Stems (TPH).
osm_tree_target <- merge(
  osm_tree[, .(SurveyID, Cycle, GTV, Stems)],
  lookup[, .(SurveyID, target_period)],
  by = "SurveyID"
)[Cycle == target_period]
gtv_by_plot <- osm_tree_target[, .(GTV_proj_m3ha = sum(GTV * Stems, na.rm = TRUE)),
                                by = SurveyID]
proj <- merge(proj, gtv_by_plot, by = "SurveyID", all.x = TRUE)

# ---- Compute observed t2 stand metrics from the FVS-style diameter_growth.csv ----
# Each tree row contributes (DIA_t2 inches, HT_t2 ft) to the t2 stand if STATUSCD_t2 == 1.
# We use the same FIA TPA design used in the OSM input (per-tree assignment by size class).
obs_alive <- observed[STATUSCD_t2 == 1 & !is.na(DIA_t2) & DIA_t2 > 0]
obs_alive[, DBH_cm := DIA_t2 * INCH_TO_CM]
obs_alive[, HT_m := ifelse(!is.na(HT_t2) & HT_t2 > 0, HT_t2 * FT_TO_M, NA_real_)]
obs_alive[, tpa_design := fifelse(DBH_cm < 12.7, 74.965,
                          fifelse(DBH_cm < 60.96, 6.018, 0.25))]
obs_alive[, stems_tph := tpa_design * ACRES_PER_HA]
obs_summary <- obs_alive[, .(
  BA_obs    = sum(0.00007854 * DBH_cm^2 * stems_tph),
  Trees_obs = sum(stems_tph),
  GTV_obs_partial = NA_real_,  # we don't have observed volume in dg.csv
  topht_obs = mean(HT_m[!is.na(HT_m) & DBH_cm >= quantile(DBH_cm, 0.85, na.rm = TRUE)],
                    na.rm = TRUE),
  n_trees_obs = .N
), by = PLT_CN]
obs_summary[, QMD_obs := sqrt((BA_obs / Trees_obs) * 40000 / pi)]
proj <- merge(proj, obs_summary, by = "PLT_CN", all.x = TRUE)
cat(sprintf("[24_osm_fia_benchmark] %d plots with both projected and observed t2 data\n",
            sum(!is.na(proj$BA_obs) & !is.na(proj$BA_proj))))

# ---- Compute metric block ----
metric_block <- function(obs, pred, name) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]; pred <- pred[ok]
  if (length(obs) < 2) return(NULL)
  resid <- pred - obs
  bias <- mean(resid)
  rmse <- sqrt(mean(resid^2))
  mae  <- mean(abs(resid))
  ss_tot <- sum((obs - mean(obs))^2)
  ss_res <- sum(resid^2)
  r2  <- 1 - ss_res / ss_tot
  obs_mean <- mean(obs)
  # 10% equivalence test (proportion of |resid| within 10% of obs_mean)
  equiv <- mean(abs(resid) <= 0.10 * obs_mean) * 100
  data.table(
    model = "OSM-ACD",
    variant = "ACD",
    variable = name,
    n        = length(obs),
    obs_mean = obs_mean,
    pred_mean = mean(pred),
    bias     = bias,
    bias_pct = 100 * bias / obs_mean,
    rmse     = rmse,
    rmse_pct = 100 * rmse / obs_mean,
    mae      = mae,
    r2       = r2,
    equiv_pct = equiv
  )
}

results <- rbindlist(list(
  metric_block(proj$BA_obs,    proj$BA_proj,    "BA_m2ha"),
  metric_block(proj$Trees_obs, proj$Trees_proj, "TPH"),
  metric_block(proj$QMD_obs,   proj$QMD_proj,   "QMD_cm"),
  metric_block(proj$topht_obs, proj$T50HT_proj, "topht_m")
))

# ---- Write outputs ----
dir.create(opt[["output-dir"]], showWarnings = FALSE, recursive = TRUE)
fwrite(results, file.path(opt[["output-dir"]], "benchmark_osm_summary.csv"))
fwrite(proj,    file.path(opt[["output-dir"]], "benchmark_osm_per_plot.csv"))
print(results)
cat(sprintf("\n[24_osm_fia_benchmark] outputs in %s\n", opt[["output-dir"]]))
