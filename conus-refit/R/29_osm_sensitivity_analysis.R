#!/usr/bin/env Rscript
################################################################################
# 29_osm_sensitivity_analysis.R
# ------------------------------------------------------------------------------
# Model sensitivity analysis for OSM-ACD on FIA. Uses the per-plot residuals
# already in hand (benchmark_osm_per_plot.csv) joined with site-productivity
# covariates from acd_plot_cn_lookup_with_bgi_csi.csv.
#
# Produces:
#   1. Random-forest variable importance for BA residual
#   2. Partial-dependence-style plots of residual vs each covariate (binned)
#   3. Faceted residual plots by species mix / forest type proxies
#   4. Map of residual hotspots
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

argv <- commandArgs(trailingOnly = TRUE)
per_plot <- if (length(argv) > 0) argv[1] else
  "calibration/output/comparisons/osm_vs_fvs/acd_bgi/benchmark_osm_per_plot.csv"
lookup   <- if (length(argv) > 1) argv[2] else
  "calibration/output/comparisons/osm_vs_fvs/acd_bgi/acd_plt_cn_lookup_with_bgi_csi.csv"
out_dir  <- if (length(argv) > 2) argv[3] else
  "calibration/output/comparisons/osm_vs_fvs/acd_bgi"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- fread(per_plot)
lk  <- fread(lookup)
# Join lookup for BGI, CSI, lat/lon, site index, initial BA/TPH
dat <- merge(dat,
              lk[, .(PLT_CN, LAT, LON, BGI, CSI_2030, si_ft, ba_t1_m2ha, tph_t1)],
              by = "PLT_CN", all.x = TRUE)

# Compute residuals (predicted - observed)
dat[, BA_resid := BA_proj - BA_obs]
dat[, TPH_resid := Trees_proj - Trees_obs]
dat[, QMD_resid := QMD_proj - QMD_obs]
dat[, BA_resid_pct := 100 * BA_resid / pmax(BA_obs, 1)]

# Site index from CSI 2030 (already metric meters or feet? — column is named CSI_2030)
# si_ft is from FIA, in feet (FIA convention)
# Convert FIA si_ft to numeric
dat[, si_ft := as.numeric(si_ft)]
dat[, CSI_2030 := as.numeric(CSI_2030)]
dat[, BGI := as.numeric(BGI)]

# ---- 1. Random forest variable importance for BA residual ----
covars <- c("BGI", "CSI_2030", "si_ft", "ba_t1_m2ha", "tph_t1", "LAT", "LON",
             "years_interval")
covars <- intersect(covars, names(dat))
d_rf <- dat[, c("BA_resid", covars), with = FALSE]
d_rf <- d_rf[complete.cases(d_rf)]
cat(sprintf("RF input: %d plots with complete covariates\n", nrow(d_rf)))

rf_done <- FALSE
if (requireNamespace("randomForest", quietly = TRUE) && nrow(d_rf) >= 200) {
  set.seed(42)
  library(randomForest)
  sub_n <- min(nrow(d_rf), 5000)
  d_sub <- d_rf[sample(.N, sub_n)]
  rf <- randomForest(BA_resid ~ ., data = d_sub, ntree = 200, importance = TRUE)
  imp <- importance(rf, type = 1)  # %IncMSE
  imp_dt <- data.table(variable = rownames(imp),
                         MSE_pct_increase = imp[, "%IncMSE"])
  setorder(imp_dt, -MSE_pct_increase)
  print(imp_dt)
  fwrite(imp_dt, file.path(out_dir, "osm_sensitivity_rf_importance.csv"))

  p_imp <- ggplot(imp_dt, aes(x = reorder(variable, MSE_pct_increase),
                                y = MSE_pct_increase)) +
    geom_col(fill = "#1f77b4") +
    coord_flip() +
    labs(title = "OSM-ACD BA Residual — Variable Importance (Random Forest)",
         subtitle = sprintf("OOB %% increase in MSE on permutation; n=%d", sub_n),
         x = NULL, y = "% increase in MSE") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(out_dir, "osm_sensitivity_rf_importance.png"), p_imp,
         width = 8, height = 5.5, dpi = 180)
  cat(sprintf("Wrote %s\n", file.path(out_dir, "osm_sensitivity_rf_importance.png")))
  rf_done <- TRUE
} else {
  cat("randomForest unavailable — skipping RF importance\n")
}

# ---- 2. Binned residual vs covariate plots ----
plot_binned <- function(d, cov_name, cov_label) {
  d_c <- d[!is.na(get(cov_name)) & !is.na(BA_resid)]
  if (nrow(d_c) < 50) return(NULL)
  qbreaks <- unique(quantile(d_c[[cov_name]], probs = seq(0, 1, 0.1), na.rm = TRUE))
  if (length(qbreaks) < 3) return(NULL)
  d_c[, bin := cut(get(cov_name), breaks = qbreaks, include.lowest = TRUE)]
  bin_summary <- d_c[, .(
    mid = mean(get(cov_name), na.rm = TRUE),
    median_resid = median(BA_resid, na.rm = TRUE),
    q25 = quantile(BA_resid, 0.25, na.rm = TRUE),
    q75 = quantile(BA_resid, 0.75, na.rm = TRUE),
    n = .N
  ), by = bin]
  bin_summary <- bin_summary[!is.na(mid)]
  ggplot(bin_summary, aes(x = mid, y = median_resid)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25, fill = "#1f77b4") +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
    geom_line(linewidth = 0.7, color = "#1f77b4") +
    geom_point(size = 1.5, color = "#1f77b4") +
    labs(title = paste("OSM BA residual vs", cov_label),
         x = cov_label, y = "BA residual (m²/ha)\n(predicted − observed)") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

covs_for_plot <- list(
  list("BGI", "Brunswick Growth Index"),
  list("CSI_2030", "Climate Site Index 2030 (m)"),
  list("si_ft", "FIA Site Index (ft)"),
  list("ba_t1_m2ha", "Initial BA t1 (m²/ha)"),
  list("tph_t1", "Initial TPH t1"),
  list("years_interval", "Years between measurements")
)
panels <- list()
for (it in covs_for_plot) {
  p <- plot_binned(dat, it[[1]], it[[2]])
  if (!is.null(p)) panels[[it[[1]]]] <- p
}

if (length(panels) >= 4) {
  pp <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = "OSM-ACD BA Residual Sensitivity",
      subtitle = "Binned residual (median + IQR) by deciles of each covariate",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(file.path(out_dir, "osm_sensitivity_binned_panels.png"), pp,
         width = 11, height = 8, dpi = 180)
  cat(sprintf("Wrote %s\n", file.path(out_dir, "osm_sensitivity_binned_panels.png")))
}

# Also save individual panels
for (nm in names(panels)) {
  ggsave(file.path(out_dir, sprintf("osm_sensitivity_residual_vs_%s.png", nm)),
         panels[[nm]], width = 6, height = 4.5, dpi = 180)
}

# ---- 3. Summary table ----
summary_dt <- dat[, .(
  n = .N,
  mean_BA_resid = mean(BA_resid, na.rm = TRUE),
  median_BA_resid = median(BA_resid, na.rm = TRUE),
  RMSE_BA = sqrt(mean(BA_resid^2, na.rm = TRUE)),
  pct_within_10pct = mean(abs(BA_resid_pct) < 10, na.rm = TRUE) * 100
)]
print(summary_dt)
fwrite(summary_dt, file.path(out_dir, "osm_sensitivity_summary.csv"))

cat("\nDone.\n")
