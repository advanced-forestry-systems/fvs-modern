# =============================================================================
# Title: Compare site productivity stress test fits (cspi vs bgi vs climate_si)
# Author: A. Weiskittel
# Date: 2026-05-13
# Description: After three fits of the same base model with different site
#              productivity variables, read the three banked summaries and
#              produce a side-by-side comparison: sigma, gamma vector, site
#              coefficient, and convergence.
#
# Run on Cardinal (login node is fine, reads CSVs only):
#   module load gcc/12.3.0 R/4.4.0
#   cd ~/fvs-modern
#   Rscript --vanilla calibration/R/eval/90b_site_stress_compare.R --model=dg_kue
# =============================================================================

library(tidyverse)
library(ggsci)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
MODEL <- get_arg("model", "dg_kue")

PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
STRESS_DIR <- file.path(PROJ_ROOT, "calibration/output/conus", MODEL, "site_stress")
OUT_DIR    <- file.path(PROJ_ROOT, "calibration/output/evaluation/site_stress")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

site_vars <- c("cspi", "bgi", "climate_si")

# --- Read summaries ---------------------------------------------------------
read_one <- function(site_var) {
  path <- file.path(STRESS_DIR, sprintf("%s_%s_b1_summary.csv", MODEL, site_var))
  if (!file.exists(path)) {
    message(sprintf("  skip: %s | file not present (%s)", site_var, path))
    return(NULL)
  }
  message(sprintf("  read: %s", site_var))
  read_csv(path, show_col_types = FALSE) |>
    mutate(site_var = site_var)
}

all_summaries <- map_dfr(site_vars, read_one)
if (nrow(all_summaries) == 0) {
  stop("No site stress summaries found. Run 90_site_productivity_stress_test.R first.")
}

# --- Sigma comparison -------------------------------------------------------
sigma_cmp <- all_summaries |>
  filter(variable == "sigma") |>
  select(site_var, mean, sd, q5, q95, rhat)
write_csv(sigma_cmp, file.path(OUT_DIR, sprintf("%s_sigma_comparison.csv", MODEL)))

cat("=== Sigma comparison ===\n")
print(sigma_cmp)

# --- Gamma comparison -------------------------------------------------------
gamma_cmp <- all_summaries |>
  filter(grepl("^gamma\\[", variable)) |>
  select(site_var, variable, mean, sd, q5, q95, rhat)
write_csv(gamma_cmp, file.path(OUT_DIR, sprintf("%s_gamma_comparison.csv", MODEL)))

# --- b6 (site coefficient) comparison ---------------------------------------
# For dg_kue, b6 is ln_csi (site coefficient). The interpretation depends on
# how each site variable was transformed; comparing b6 magnitudes across
# site_var is informative as direction and significance.
site_coef <- all_summaries |>
  filter(variable == "b6") |>
  select(site_var, mean, sd, q5, q95, rhat)
write_csv(site_coef, file.path(OUT_DIR, sprintf("%s_site_coef_comparison.csv", MODEL)))

cat("\n=== Site coefficient (b6) comparison ===\n")
print(site_coef)

# --- Figure: sigma + site_coef + gamma side-by-side ------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold")
  )

p_sigma <- ggplot(sigma_cmp, aes(site_var, mean)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2, color = "#3B4992") +
  geom_point(size = 3, color = "#3B4992") +
  labs(x = "Site productivity variable", y = "Residual sigma posterior mean",
       title = "A) Residual sigma by site variable") +
  theme_pub

p_site <- ggplot(site_coef, aes(site_var, mean)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2, color = "#EE0000") +
  geom_point(size = 3, color = "#EE0000") +
  labs(x = "Site productivity variable", y = "b6 (site coef) posterior mean",
       title = "B) Site coefficient (b6) by site variable") +
  theme_pub

gamma_labels <- c(
  "gamma[1]" = "Wood specific gravity",
  "gamma[2]" = "Shade tolerance",
  "gamma[3]" = "Softwood",
  "gamma[4]" = "Leaf longevity",
  "gamma[5]" = "Max height",
  "gamma[6]" = "Max DBH",
  "gamma[7]" = "Vulnerability score",
  "gamma[8]" = "Climate sensitivity"
)
gamma_plot_df <- gamma_cmp |>
  mutate(trait = gamma_labels[variable])

p_gamma <- ggplot(gamma_plot_df, aes(mean, trait, color = site_var)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0, linewidth = 0.6,
                 position = position_dodge(0.4)) +
  geom_point(size = 2.2, position = position_dodge(0.4)) +
  scale_color_aaas() +
  labs(x = "Gamma posterior mean (90% CI)", y = NULL,
       color = "Site variable",
       title = "C) Trait coefficients (gamma) by site variable") +
  theme_pub

combined <- (p_sigma + p_site) / p_gamma +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    title = sprintf("Site productivity stress test: %s B1 species-free",
                    toupper(MODEL)),
    subtitle = "Same model refit with three different site variables. Comparable sigma + consistent gamma direction = robust."
  )

fig_png <- file.path(OUT_DIR, sprintf("%s_site_stress_panels.png", MODEL))
fig_pdf <- file.path(OUT_DIR, sprintf("%s_site_stress_panels.pdf", MODEL))
ggsave(fig_png, combined, width = 22, height = 22, units = "cm",
       dpi = 300, bg = "white")
ggsave(fig_pdf, combined, width = 22, height = 22, units = "cm")
cat("Wrote", fig_png, "\n")

# --- Verdict ---------------------------------------------------------------
sigma_max <- max(sigma_cmp$mean, na.rm = TRUE)
sigma_min <- min(sigma_cmp$mean, na.rm = TRUE)
sigma_range_pct <- 100 * (sigma_max - sigma_min) / sigma_min

cat("\n=== Verdict ===\n")
cat(sprintf("Sigma range across site variables: %.2f to %.2f (%.1f%% spread)\n",
            sigma_min, sigma_max, sigma_range_pct))
if (sigma_range_pct < 5) {
  cat("VERDICT: Sigma is stable across site variables (within 5%%).\n")
  cat("The model is robust to site productivity metric choice.\n")
} else if (sigma_range_pct < 15) {
  cat("VERDICT: Modest sigma sensitivity to site variable choice.\n")
  cat("CSPI remains a defensible default but discuss in manuscript.\n")
} else {
  cat("VERDICT: Sigma varies meaningfully across site variables.\n")
  cat("The choice of CSPI is material to the manuscript's conclusions.\n")
}

# Direction consistency of gamma across site variables
gamma_signs <- gamma_cmp |>
  mutate(sign = sign(mean)) |>
  group_by(variable) |>
  summarise(n_signs = n_distinct(sign), .groups = "drop")
n_flips <- sum(gamma_signs$n_signs > 1)
cat(sprintf("\nGamma sign flips across site variables: %d of %d trait coefficients\n",
            n_flips, nrow(gamma_signs)))

cat("\nDone.\n")
