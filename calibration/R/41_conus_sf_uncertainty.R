#!/usr/bin/env Rscript
# 41_conus_sf_uncertainty.R
# Registers the unified species-free posterior draws + CSPI site-surface draws with the
# existing 21_uncertainty_propagation.R / Bakuzis credible-band machinery.
#
# Two uncertainty sources, combined by Monte Carlo over draw index j:
#   (1) parameter draws : config$components[[c]]$posterior$path  (thinned parquet)
#   (2) site draws      : CSPI v7 QRF predictive interval per plot
#
# STATUS: scaffold. Reuses 21_uncertainty_propagation.R; only the source loaders are new.

suppressMessages({ library(arrow); library(data.table) })

load_param_draws <- function(cfg_path = "config/calibrated/conus_sf.json") {
  stop("TODO: read each component posterior parquet referenced in the config")
}
draw_site_index <- function(plot_xy, j, qrf_path) {
  # TODO: sample SI from the v7 QRF predictive interval for these plots, draw j
  stop("TODO: implement QRF predictive draw")
}
# hand the combined draws to 21_uncertainty_propagation.R's projection loop
