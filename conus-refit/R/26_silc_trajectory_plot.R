#!/usr/bin/env Rscript
################################################################################
# 26_silc_trajectory_plot.R
# Plot OSM-ACD vs FVS-NE 100-year trajectories on Seven Islands managed-forest
# stands. Reads trajectory_long.csv from compare_silc_osm_vs_fvsne.py.
################################################################################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

argv <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(argv) > 0) argv[1] else
  "calibration/output/comparisons/osm_vs_fvs/seven_islands_100yr/trajectory_long.csv"
out_dir <- if (length(argv) > 1) argv[2] else
  "calibration/output/comparisons/osm_vs_fvs/seven_islands_100yr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- fread(in_csv)

# Filter to only stands present in BOTH models (some stands are OSM-only or FVS-only)
stands_both <- dat[, .N, by = .(STAND_ID, model)][, .(n_models = .N), by = STAND_ID][n_models == 2, STAND_ID]
dat <- dat[STAND_ID %in% stands_both]

# Make a long form for the 4 metrics
dat_long <- melt(dat,
                 id.vars = c("model", "STAND_ID", "Cycle", "Year"),
                 measure.vars = intersect(
                   c("BA_m2ha", "TPH", "QMD_cm", "GMV_m3ha",
                     "MerchVol_m3ha", "RD",
                     "Biomass_Mgha", "Carbon_Mgha"),
                   names(dat)
                 ),
                 variable.name = "metric", value.name = "value")
dat_long[, metric := factor(metric,
  levels = c("BA_m2ha", "TPH", "QMD_cm", "GMV_m3ha", "MerchVol_m3ha", "RD",
              "Biomass_Mgha", "Carbon_Mgha"),
  labels = c("Basal Area (m²/ha)", "TPH (trees/ha)",
              "QMD (cm)", "Gross Volume (m³/ha)",
              "Merch Volume 1ft-4in top (m³/ha)",
              "Relative Density (RD)",
              "Biomass (Mg/ha)", "Carbon (Mg/ha)"))]
dat_long[, model := factor(model, levels = c("OSM-ACD", "FVS-NE"))]

# Plot 1: Mean trajectories with stand-level ribbon
mean_traj <- dat_long[, .(value_mean = mean(value, na.rm = TRUE),
                          q25 = quantile(value, 0.25, na.rm = TRUE),
                          q75 = quantile(value, 0.75, na.rm = TRUE),
                          n_stands = .N),
                       by = .(model, Year, metric)]

p1 <- ggplot(mean_traj, aes(x = Year, y = value_mean, color = model, fill = model)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("OSM-ACD" = "#1f77b4", "FVS-NE" = "#d62728")) +
  scale_fill_manual(values = c("OSM-ACD" = "#1f77b4", "FVS-NE" = "#d62728")) +
  labs(title = "OSM-ACD vs FVS-NE 100-year projection on Seven Islands managed forest",
       subtitle = sprintf("%d stands; mean line, IQR ribbon across stands",
                           length(stands_both)),
       x = "Year", y = NULL, color = "Model", fill = "Model") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "silc_100yr_trajectories.png"), p1,
       width = 9, height = 7, dpi = 200)
ggsave(file.path(out_dir, "silc_100yr_trajectories.pdf"), p1,
       width = 9, height = 7)
cat(sprintf("Wrote %s and .pdf\n",
            file.path(out_dir, "silc_100yr_trajectories.png")))

# Plot 2: Per-stand trajectories of BA, color by model
ba_by_stand <- dat_long[metric == "Basal Area (m²/ha)"]
p2 <- ggplot(ba_by_stand, aes(x = Year, y = value,
                               color = model, group = interaction(model, STAND_ID))) +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  facet_wrap(~ STAND_ID, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("OSM-ACD" = "#1f77b4", "FVS-NE" = "#d62728")) +
  labs(title = "Per-stand BA trajectories — OSM-ACD vs FVS-NE",
       x = "Year", y = "Basal Area (m²/ha)", color = "Model") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "silc_100yr_BA_by_stand.png"), p2,
       width = 10, height = 7, dpi = 200)
cat(sprintf("Wrote %s\n", file.path(out_dir, "silc_100yr_BA_by_stand.png")))

# Plot 3: OSM-FVS-NE divergence over time per metric
diverge <- dcast(dat_long, STAND_ID + Cycle + Year + metric ~ model,
                 value.var = "value")
diverge[, diff := `OSM-ACD` - `FVS-NE`]
diverge[, pct_diff := 100 * diff / `FVS-NE`]
divg_summary <- diverge[, .(mean_pct = mean(pct_diff, na.rm = TRUE)),
                         by = .(metric, Year)]
p3 <- ggplot(divg_summary, aes(x = Year, y = mean_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
  geom_line(linewidth = 0.9, color = "#2ca02c") +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(title = "OSM-ACD minus FVS-NE divergence over time",
       subtitle = "Mean % difference across stands. Negative = OSM-ACD lower than FVS-NE.",
       x = "Year", y = "Mean % difference (OSM-ACD − FVS-NE)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "silc_100yr_divergence.png"), p3,
       width = 9, height = 7, dpi = 200)
cat(sprintf("Wrote %s\n", file.path(out_dir, "silc_100yr_divergence.png")))

# Plot 4: SDI vs TPH (stand density framework) per stand, colored by model
if (all(c("SDI", "TPH") %in% names(dat))) {
  sdmd <- dat[!is.na(SDI) & SDI > 0 & TPH > 0]
  # Reineke maximum density reference line (assume SDImax ~ 1200 metric for Acadian)
  sdimax_ref <- 1200
  ref_line <- data.table(TPH = seq(100, 20000, length.out = 200))
  ref_line[, SDI_max := sdimax_ref]
  ref_line[, SDI_60 := 0.6 * sdimax_ref]
  ref_line[, SDI_35 := 0.35 * sdimax_ref]
  p4 <- ggplot(sdmd, aes(x = TPH, y = SDI, color = model, group = interaction(model, STAND_ID))) +
    geom_path(alpha = 0.6, linewidth = 0.5,
              arrow = arrow(angle = 18, length = unit(0.10, "inches"), type = "closed")) +
    geom_point(data = sdmd[Cycle == 0], size = 1.4, shape = 16) +
    geom_hline(yintercept = sdimax_ref, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0.6 * sdimax_ref, linetype = "dotted", color = "grey60") +
    geom_hline(yintercept = 0.35 * sdimax_ref, linetype = "dotted", color = "grey60") +
    annotate("text", x = 100, y = sdimax_ref, label = "SDImax (1200)",
              hjust = 0, vjust = -0.3, size = 3, color = "grey40") +
    annotate("text", x = 100, y = 0.6 * sdimax_ref, label = "60% (full stocking)",
              hjust = 0, vjust = -0.3, size = 3, color = "grey50") +
    annotate("text", x = 100, y = 0.35 * sdimax_ref, label = "35% (management lower)",
              hjust = 0, vjust = -0.3, size = 3, color = "grey50") +
    scale_x_log10() +
    scale_color_manual(values = c("OSM-ACD" = "#1f77b4", "FVS-NE" = "#d62728")) +
    labs(title = "SDI vs TPH — Stand Density Framework over 100-yr Projection",
         subtitle = "Per-stand trajectory paths; dots = year 0; arrows show direction over time",
         x = "Trees per hectare (log scale)", y = "Stand Density Index (Reineke metric)",
         color = "Model") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "silc_100yr_SDI_vs_TPH.png"), p4,
         width = 9, height = 6.5, dpi = 200)
  cat(sprintf("Wrote %s\n", file.path(out_dir, "silc_100yr_SDI_vs_TPH.png")))
}

# Plot 5: Merch volume detail (1ft stump, 4in top) — comparable industry metric
if ("MerchVol_m3ha" %in% names(dat)) {
  mv <- dat_long[metric == "Merch Volume 1ft-4in top (m³/ha)"]
  mv_summary <- mv[, .(value_mean = mean(value, na.rm = TRUE),
                       q25 = quantile(value, 0.25, na.rm = TRUE),
                       q75 = quantile(value, 0.75, na.rm = TRUE)),
                    by = .(model, Year)]
  p5 <- ggplot(mv_summary, aes(x = Year, y = value_mean, color = model, fill = model)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c("OSM-ACD" = "#1f77b4", "FVS-NE" = "#d62728")) +
    scale_fill_manual(values = c("OSM-ACD" = "#1f77b4", "FVS-NE" = "#d62728")) +
    labs(title = "Merchantable Volume Trajectory — Seven Islands",
         subtitle = "1-ft stump to 4-inch top diameter inside bark, NSVB-derived",
         x = "Year", y = "Merch Volume (m³/ha)",
         color = "Model", fill = "Model") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "silc_100yr_merch_volume.png"), p5,
         width = 8, height = 5.5, dpi = 200)
  cat(sprintf("Wrote %s\n", file.path(out_dir, "silc_100yr_merch_volume.png")))
}
