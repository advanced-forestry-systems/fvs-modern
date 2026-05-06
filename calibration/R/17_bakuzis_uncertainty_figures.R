# =============================================================================
# Title: Bakuzis matrix uncertainty figures
# Author: A. Weiskittel
# Date: 2026-04-23
# Description: Publication quality figures from the Bakuzis uncertainty
#   pipeline. Consumes the long summary frame produced by
#   calibration/python/bakuzis_uncertainty_aggregate.py and renders four
#   figure panels:
#     1. 100 year trajectories with default, calibrated MAP, and posterior
#        95 percent credible band per scenario grid
#     2. Year 100 divergence scatter with band width as a size aesthetic
#     3. Bakuzis four laws compliance panel for default vs calibrated
#     4. Posterior band width over horizon, one line per scenario grouping
#
#   Inputs expected under calibration/output/bakuzis/
#     - bakuzis_uncertainty_summary_long.csv
#     - bakuzis_benchmark_wide.csv
#     - bakuzis_laws_compliance.csv
#
#   Output directory: calibration/output/comparisons/manuscript_figures/
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

# --- Paths -------------------------------------------------------------------

project_root <- Sys.getenv(
  "FVS_PROJECT_ROOT",
  normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."),
                mustWork = FALSE)
)
bakuzis_dir <- file.path(project_root, "calibration", "output", "bakuzis")
fig_dir <- file.path(project_root, "calibration", "output",
                     "comparisons", "manuscript_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

summary_path <- file.path(bakuzis_dir, "bakuzis_uncertainty_summary_long.csv")
bench_path <- file.path(bakuzis_dir, "bakuzis_benchmark_wide.csv")
laws_path <- file.path(bakuzis_dir, "bakuzis_laws_compliance.csv")

stopifnot(file.exists(summary_path))

# --- Theme -------------------------------------------------------------------

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 13, face = "bold")
  )

config_colors <- c(
  "default" = "#999999",
  "calibrated_map" = "#d95f02",
  "posterior_band" = "#1b9e77"
)
config_labels <- c(
  "default" = "Default parameters",
  "calibrated_map" = "Calibrated MAP",
  "posterior_band" = "Calibrated posterior (95 percent band)"
)

# --- Load --------------------------------------------------------------------

summary_long <- read_csv(summary_path, show_col_types = FALSE) %>%
  mutate(
    config = factor(config, levels = names(config_colors)),
    species_group = factor(species_group),
    site_class = factor(site_class, levels = c("Low", "Medium", "High")),
    density_class = factor(density_class, levels = c("Low", "Medium", "High"))
  ) %>%
  group_by(variant, scenario) %>%
  mutate(horizon_yr = year - min(year)) %>%
  ungroup()

# --- Figure 1: trajectory grid with bands ------------------------------------

ba <- summary_long %>% filter(variable == "atba")

p1 <- ggplot(ba, aes(x = horizon_yr)) +
  geom_ribbon(
    data = ba %>% filter(config == "posterior_band"),
    aes(ymin = q025, ymax = q975, fill = config),
    alpha = 0.30
  ) +
  geom_line(
    aes(y = median, color = config, linetype = config),
    linewidth = 0.8
  ) +
  scale_color_manual(values = config_colors, labels = config_labels) +
  scale_fill_manual(values = config_colors, labels = config_labels) +
  scale_linetype_manual(
    values = c("default" = "dashed", "calibrated_map" = "solid",
               "posterior_band" = "solid"),
    labels = config_labels
  ) +
  facet_grid(site_class + density_class ~ species_group,
             labeller = labeller(.multi_line = FALSE)) +
  labs(
    title = "Bakuzis matrix: 100 year BA trajectories with parametric uncertainty",
    subtitle = "Ribbon is 95 percent posterior credible band; dashed line is default",
    x = "Projection year",
    y = expression(paste("Basal area (ft"^2, "/ac)")),
    color = "Configuration", fill = "Configuration", linetype = "Configuration"
  ) +
  theme_pub +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_dir, "fig_bakuzis_trajectories.pdf"),
  p1, width = 12, height = 10, device = cairo_pdf
)
ggsave(
  file.path(fig_dir, "fig_bakuzis_trajectories.png"),
  p1, width = 12, height = 10, dpi = 300
)

# --- Figure 2: year 100 divergence scatter -----------------------------------

if (file.exists(bench_path)) {
  bench <- read_csv(bench_path, show_col_types = FALSE)
  ba_wide <- bench %>%
    filter(variable == "atba") %>%
    mutate(
      band_pct = 100 * post_band_y100 / post_median_y100,
      diff_pct = diff_pct_y100
    )

  p2 <- ggplot(ba_wide,
               aes(x = def_median_y100, y = diff_pct,
                   color = species_group, shape = site_class,
                   size = band_pct)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 8),
                          name = "95 percent band width\n(% of posterior median)") +
    labs(
      title = "Year 100 BA divergence: calibrated posterior median vs default",
      x = expression(paste("Default year 100 BA (ft"^2, "/ac)")),
      y = "Calibrated minus default (percent of default)",
      color = "Species group", shape = "Site class"
    ) +
    theme_pub

  ggsave(file.path(fig_dir, "fig_bakuzis_divergence.pdf"),
         p2, width = 9, height = 7, device = cairo_pdf)
  ggsave(file.path(fig_dir, "fig_bakuzis_divergence.png"),
         p2, width = 9, height = 7, dpi = 300)
}

# --- Figure 3: Bakuzis four laws ---------------------------------------------

if (file.exists(laws_path)) {
  laws <- read_csv(laws_path, show_col_types = FALSE)
  laws_long <- laws %>%
    pivot_longer(starts_with("law_"),
                 names_to = "law", values_to = "fraction") %>%
    mutate(
      law = recode(
        law,
        law_sukachev_fraction = "Sukachev effect",
        law_eichhorn_fraction = "Eichhorn's rule",
        law_density_mortality_fraction = "Density drives mortality",
        law_mortality_size_fraction = "Mortality size pattern"
      )
    )

  p3 <- ggplot(laws_long, aes(x = law, y = fraction, fill = config)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = config_colors, labels = config_labels) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    facet_wrap(~ variant, ncol = 2) +
    labs(
      title = "Bakuzis biological law compliance at year 100",
      x = NULL, y = "Fraction of scenario pairs satisfying law",
      fill = "Configuration"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(fig_dir, "fig_bakuzis_laws.pdf"),
         p3, width = 10, height = 6, device = cairo_pdf)
  ggsave(file.path(fig_dir, "fig_bakuzis_laws.png"),
         p3, width = 10, height = 6, dpi = 300)
}

# --- Figure 4: band width over horizon ---------------------------------------

band_over_time <- summary_long %>%
  filter(variable == "atba", config == "posterior_band") %>%
  mutate(band_width = q975 - q025,
         band_pct = 100 * band_width / pmax(median, 1e-6))

p4 <- ggplot(band_over_time,
             aes(x = horizon_yr, y = band_pct,
                 color = site_class, group = scenario)) +
  geom_line(alpha = 0.5, linewidth = 0.6) +
  facet_wrap(~ species_group, ncol = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Posterior uncertainty growth with projection horizon",
    x = "Projection year",
    y = "BA 95 percent band width (% of posterior median)",
    color = "Site class"
  ) +
  theme_pub

ggsave(file.path(fig_dir, "fig_bakuzis_band_growth.pdf"),
       p4, width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig_bakuzis_band_growth.png"),
       p4, width = 10, height = 7, dpi = 300)

cat("\nFigures written to:", fig_dir, "\n")
