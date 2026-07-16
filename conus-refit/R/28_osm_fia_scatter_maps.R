#!/usr/bin/env Rscript
################################################################################
# 28_osm_fia_scatter_maps.R
# ------------------------------------------------------------------------------
# Build publication-quality FIA benchmark plots for OSM-ACD:
#   1. Scatter: predicted vs observed for BA, TPH, QMD, top height (4-panel)
#   2. Residual maps: lat/lon hexbin of (pred - obs) for BA and biomass
#   3. Composite dashboard combining trajectory + scatter + residual
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
# Join to get LAT/LON for residual maps
dat <- merge(dat, lk[, .(PLT_CN, LAT, LON)],
              by.x = "PLT_CN", by.y = "PLT_CN", all.x = TRUE)

# ---- 1. Predicted vs observed scatter (4 panels) ----
plot_one <- function(obs_col, pred_col, label, units) {
  d <- dat[!is.na(get(obs_col)) & !is.na(get(pred_col))]
  if (nrow(d) < 10) {
    return(ggplot() + labs(title = paste("OSM-ACD —", label, "(insufficient data)")) +
             theme_minimal())
  }
  obs <- d[[obs_col]]
  pred <- d[[pred_col]]
  resid <- pred - obs
  rmse <- sqrt(mean(resid^2))
  bias <- mean(resid)
  r2 <- 1 - sum(resid^2) / sum((obs - mean(obs))^2)
  rmse_pct <- 100 * rmse / mean(obs)
  bias_pct <- 100 * bias / mean(obs)
  ann <- sprintf("n=%d | R²=%.3f\nbias=%+.2f (%+.1f%%)\nRMSE=%.2f (%.1f%%)",
                 nrow(d), r2, bias, bias_pct, rmse, rmse_pct)
  d_plot <- data.table(obs = obs, pred = pred)
  ggplot(d_plot, aes(x = obs, y = pred)) +
    geom_hex(bins = 60) +
    geom_abline(intercept = 0, slope = 1, color = "red",
                  linetype = "dashed", linewidth = 0.4) +
    scale_fill_gradient(low = "lightblue", high = "darkblue",
                         trans = "log10") +
    annotate("label", x = -Inf, y = Inf, label = ann,
              hjust = -0.05, vjust = 1.1, size = 3, alpha = 0.85) +
    labs(title = paste("OSM-ACD —", label),
         x = paste("Observed", label, units),
         y = paste("Predicted", label, units)) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 11))
}

p_ba   <- plot_one("BA_obs",   "BA_proj",   "Basal Area", "(m²/ha)")
p_tph  <- plot_one("Trees_obs","Trees_proj","TPH",        "(trees/ha)")
p_qmd  <- plot_one("QMD_obs",  "QMD_proj",  "QMD",        "(cm)")
p_topht<- plot_one("topht_obs","T50HT_proj","Top Height", "(m)")

# Save each panel individually plus a combined 2x2 grid using cowplot
for (item in list(list("BA", p_ba), list("TPH", p_tph),
                   list("QMD", p_qmd), list("topHT", p_topht))) {
  out_f <- file.path(out_dir, sprintf("osm_fia_scatter_%s.png", item[[1]]))
  tryCatch({
    ggsave(out_f, item[[2]], width = 6, height = 5, dpi = 180)
    cat(sprintf("Wrote %s\n", out_f))
  }, error = function(e) cat(sprintf("Skipped %s: %s\n", out_f, e$message)))
}
# Combined grid (skip on patchwork error)
scatter_grid <- tryCatch({
  (p_ba + p_tph) / (p_qmd + p_topht) +
    plot_annotation(title = "OSM-ACD FIA Benchmark: Predicted vs Observed",
                     subtitle = sprintf("%d ACD plot pairs", nrow(dat)))
}, error = function(e) NULL)
if (!is.null(scatter_grid)) {
  tryCatch({
    ggsave(file.path(out_dir, "osm_fia_scatter_4panel.png"),
           scatter_grid, width = 11, height = 9, dpi = 180)
    cat(sprintf("Wrote %s\n", file.path(out_dir, "osm_fia_scatter_4panel.png")))
  }, error = function(e) cat(sprintf("Composite skipped: %s\n", e$message)))
}

# ---- 2. Residual maps (BA) ----
dat[, BA_resid := BA_proj - BA_obs]
dat[, BA_resid_pct := 100 * BA_resid / pmax(BA_obs, 1)]
map_dat <- dat[!is.na(LAT) & !is.na(LON) & !is.na(BA_resid)]

p_map_abs <- ggplot(map_dat, aes(x = LON, y = LAT, z = BA_resid)) +
  stat_summary_hex(fun = "median", bins = 60) +
  scale_fill_gradient2(low = "#1f77b4", mid = "#ffffff", high = "#d62728",
                       midpoint = 0, limits = c(-10, 10),
                       oob = scales::squish,
                       name = "Median residual\n(m²/ha)") +
  coord_quickmap() +
  labs(title = "BA Residual (OSM predicted − observed) by plot location",
        x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "osm_fia_residual_map_BA.png"), p_map_abs,
       width = 10, height = 7, dpi = 180)
cat(sprintf("Wrote %s\n", file.path(out_dir, "osm_fia_residual_map_BA.png")))

# ---- 3. Composite dashboard ----
if (!is.null(scatter_grid)) {
  tryCatch({
    dashboard <- scatter_grid / p_map_abs +
      plot_layout(heights = c(2, 1.4)) +
      plot_annotation(title = "OSM-ACD vs FIA Observed — Benchmark Dashboard")
    ggsave(file.path(out_dir, "osm_fia_dashboard.png"), dashboard,
           width = 12, height = 13, dpi = 180)
    cat(sprintf("Wrote %s\n", file.path(out_dir, "osm_fia_dashboard.png")))
  }, error = function(e) cat(sprintf("Dashboard skipped: %s\n", e$message)))
}

cat("\nDone.\n")
