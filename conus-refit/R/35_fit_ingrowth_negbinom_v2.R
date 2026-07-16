##=============================================================================
## 35_fit_ingrowth_negbinom_v2.R
##
## Fixed ingrowth driver: uses FIA TREE_GRM_COMPONENT.csv files (which have
## explicit COMPONENT = 'INGROWTH' rows) instead of the matched-pairs file
## (which has zero ingrowth by construction).
##
## Currently available state files (May 15 2026):
##   AR_TREE_GRM_COMPONENT.csv (~654k rows)
##   CA_TREE_GRM_COMPONENT.csv (~213k rows)
## ~13k INGROWTH events in CA alone. Pilot on these two states.
## Full CONUS coverage requires downloading the remaining 46 state files.
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-15
##=============================================================================

library(data.table)
library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/ingrowth_negbinom.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/ingrowth")
OUT_NAME  <- get_arg("outname",   "ingrowth_negbinom_v2")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 35_fit_ingrowth_negbinom_v2.R ==\n\n")

## 1. Load FIA GRM component files (ingrowth source) -------------------------
GRM_DIR <- "calibration/data/raw_fia"
grm_files <- list.files(GRM_DIR, pattern = "_TREE_GRM_COMPONENT\\.csv$",
                        full.names = TRUE)
cat("GRM component files found:", length(grm_files), "\n")
print(basename(grm_files))

# Select only needed columns to save memory
gr_cols <- c("TRE_CN", "PREV_TRE_CN", "PLT_CN", "STATECD",
             "DIA_BEGIN", "DIA_MIDPT", "DIA_END",
             "SUBP_COMPONENT_AL_FOREST", "SUBP_TPAGROW_UNADJ_AL_FOREST")

cat("\nLoading GRM data ..."); flush.console()
grm_list <- lapply(grm_files, function(f) {
  fread(f, select = gr_cols, showProgress = FALSE)
})
grm <- rbindlist(grm_list)
cat(" done. Total rows:", nrow(grm), "\n")

# Filter to INGROWTH plus get a plot-level base set
ingrowth_events <- grm[SUBP_COMPONENT_AL_FOREST == "INGROWTH"]
cat("INGROWTH events:", nrow(ingrowth_events), "\n")

## 2. Aggregate to plot-level ingrowth counts --------------------------------
# Sum TPA-weighted ingrowth per plot
plot_ingrowth <- ingrowth_events[, .(
  n_recruits_tpa = sum(SUBP_TPAGROW_UNADJ_AL_FOREST, na.rm = TRUE),
  n_recruits     = .N
), by = .(PLT_CN, STATECD)]

# Also build the full plot set (includes plots with zero ingrowth)
all_plots <- unique(grm[, .(PLT_CN, STATECD)])
all_plots <- merge(all_plots, plot_ingrowth, by = c("PLT_CN", "STATECD"),
                   all.x = TRUE)
all_plots[is.na(n_recruits), n_recruits := 0L]
all_plots[is.na(n_recruits_tpa), n_recruits_tpa := 0]

cat("Total plots:", nrow(all_plots), "\n")
cat("Plots with >= 1 ingrowth:", sum(all_plots$n_recruits > 0), "\n")
cat("Mean ingrowth per plot:", round(mean(all_plots$n_recruits), 2), "\n")
cat("Max ingrowth per plot:", max(all_plots$n_recruits), "\n\n")

## 3. Join with the matched-pairs file for plot covariates -------------------
pairs_path <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
if (file.exists(pairs_path)) {
  cat("Loading pairs data for covariates ..."); flush.console()
  trees <- as.data.table(readRDS(pairs_path))
  cat(" done\n")

  # The pairs file has CN_cond1 and plot_key. We need to match PLT_CN to one of these
  # Pairs are at the TREE level; aggregate to plot level for covariates.
  plot_covs <- trees[, .(
    ba_t1       = first(BA1),
    bal_mean_t1 = mean(BAL1, na.rm = TRUE),
    qmd_t1      = first(QMD1),
    sdi_t1      = first(SDI1),
    rd_t1       = first(RD1),
    stand_age   = first(STDAGE),
    years       = first(YEARS),
    cspi        = first(cspi),
    bgi         = first(bgi),
    climate_si  = first(climate_si),
    clim_pca1   = first(clim_pca1),
    EPA_L1_CODE = first(EPA_L1_CODE),
    EPA_L2_CODE = first(EPA_L2_CODE),
    EPA_L3_CODE = first(EPA_L3_CODE),
    dom_spcd    = SPCD[which.max(DBH1 * (TREESTATUS1 == 1))]
  ), by = plot_key]

  # Note: plot_key in pairs file != PLT_CN in GRM. We can't merge cleanly without
  # a plot ID translation table. For the pilot, work with whichever plots happen
  # to be in both via state code matching as a rough approximation.
  cat("WARNING: plot ID matching between GRM (PLT_CN) and pairs (plot_key)\n")
  cat("not implemented. The pilot fits without covariates from the pairs file.\n\n")
}

## 4. Quick analysis: distribution of ingrowth counts ------------------------
cat("=== Ingrowth count distribution ===\n")
print(table(cut(all_plots$n_recruits,
                breaks = c(-0.5, 0.5, 1.5, 3.5, 7.5, Inf),
                labels = c("0", "1", "2-3", "4-7", "8+"))))
cat("\n")
cat("Total ingrowth events (tree-level):", nrow(ingrowth_events), "\n")
cat("Total plots with ingrowth:", sum(all_plots$n_recruits > 0), "\n")
cat("Aggregate rate (ingrowth / plot):",
    round(sum(all_plots$n_recruits) / nrow(all_plots), 3), "\n")

## 5. Save pilot data --------------------------------------------------------
data_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_plot_data.csv"))
fwrite(all_plots, data_path)
cat("\nSaved plot data to:", data_path, "\n")
cat("\nNote: This is a *pilot* on AR + CA only. Full CONUS ingrowth model\n")
cat("requires downloading TREE_GRM_COMPONENT.csv for the remaining 46 states\n")
cat("from the FIA DataMart (https://apps.fs.usda.gov/fia/datamart/datamart.html).\n")
cat("\nDone.\n")
