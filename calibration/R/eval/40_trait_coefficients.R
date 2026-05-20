# =============================================================================
# Title: Trait coefficient forest plot across base models
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: Visualizes the eight-trait fixed-effect coefficient vector
#              (gamma) for each base model where the trait architecture is
#              implemented. Outputs a forest plot showing which traits matter
#              for which response, plus a wide-format table for the manuscript.
#
# Run on Cardinal:
#   module load gcc/12.3.0 R/4.4.0
#   cd ~/fvs-modern
#   Rscript --vanilla calibration/R/eval/40_trait_coefficients.R
#
# Reads each model's summary CSV (low memory; does not load the fit objects).
# The trait architecture is currently in HG, HT-DBH, DG Kuehne. CR, HCB, and
# Mortality use random intercepts only; they are excluded from this figure.
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggsci)

# --- Paths and constants ----------------------------------------------------
PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR   <- file.path(PROJ_ROOT, "calibration/output/evaluation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

trait_labels <- c(
  "gamma[1]" = "Wood specific gravity",
  "gamma[2]" = "Shade tolerance",
  "gamma[3]" = "Softwood (indicator)",
  "gamma[4]" = "Leaf longevity (months)",
  "gamma[5]" = "Max height (m)",
  "gamma[6]" = "Max DBH (cm)",
  "gamma[7]" = "Vulnerability score",
  "gamma[8]" = "Climate sensitivity"
)

models_with_traits <- tribble(
  ~model,        ~variant,            ~path,
  "HG",          "B2 species-aware",  "calibration/output/evaluation/summaries/hg_b2_eval_summary.csv",
  "HG",          "B1 species-free",   "calibration/output/evaluation/summaries/hg_b1_eval_summary.csv",
  "HT-DBH",      "B2 species-aware",  "calibration/output/evaluation/summaries/htdbh_b2_eval_summary.csv",
  "DG Kuehne",   "B2 species-aware",  "calibration/output/evaluation/summaries/dg_kue_b2_eval_summary.csv",
  "DG Kuehne",   "B1 species-free",   "calibration/output/evaluation/summaries/dg_kue_b1_eval_summary.csv"
) %>%
  mutate(full_path = file.path(PROJ_ROOT, path))

# --- Read and reshape -------------------------------------------------------

read_gammas <- function(model, variant, full_path) {
  if (!file.exists(full_path)) {
    message(sprintf("  skip: %s | %s (file not present)", model, variant))
    return(NULL)
  }
  message(sprintf("  read: %s | %s", model, variant))
  df <- read_csv(full_path, show_col_types = FALSE)
  df %>%
    filter(grepl("^gamma\\[", variable)) %>%
    mutate(
      model        = model,
      variant      = variant,
      trait_label  = trait_labels[variable],
      excludes_zero = (q5 * q95) > 0
    )
}

gamma_long <- models_with_traits %>%
  pmap_dfr(~read_gammas(..1, ..2, ..4))

if (nrow(gamma_long) == 0) {
  stop("No trait coefficient summaries found. Run base model fits first.")
}

# --- Forest plot ------------------------------------------------------------

theme_forest <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.title         = element_text(face = "bold")
  )

p_forest <- gamma_long %>%
  mutate(model_variant = factor(paste(model, variant, sep = " | "))) %>%
  ggplot(aes(x = mean, y = trait_label, color = variant, shape = variant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  facet_wrap(~model, ncol = 1, scales = "free_x") +
  scale_color_aaas() +
  scale_shape_manual(values = c("B1 species-free" = 17, "B2 species-aware" = 16)) +
  labs(
    x = "Trait coefficient (gamma, log scale)",
    y = NULL,
    color = "Variant", shape = "Variant",
    title = "Trait fixed-effect coefficients across base models",
    subtitle = "Points are posterior means, bars are 90% credible intervals"
  ) +
  theme_forest

# --- Save plots -------------------------------------------------------------

n_models <- length(unique(gamma_long$model))
fig_height <- 6 + 6 * n_models
fig_png <- file.path(OUT_DIR, "trait_coefficients_forest.png")
fig_pdf <- file.path(OUT_DIR, "trait_coefficients_forest.pdf")

ggsave(fig_png, p_forest,
       width = 18, height = fig_height, units = "cm", dpi = 300, bg = "white")
ggsave(fig_pdf, p_forest,
       width = 18, height = fig_height, units = "cm")
cat("Wrote", fig_png, "\n")
cat("Wrote", fig_pdf, "\n")

# --- Wide-format table for manuscript ---------------------------------------

gamma_wide <- gamma_long %>%
  mutate(
    model_variant = paste(model, str_extract(variant, "B[12]"), sep = " "),
    ci_str = sprintf("%+.3f [%+.3f, %+.3f]", mean, q5, q95),
    ci_str = ifelse(excludes_zero, paste0("**", ci_str, "**"), ci_str)
  ) %>%
  select(trait_label, model_variant, ci_str) %>%
  pivot_wider(names_from = model_variant, values_from = ci_str)

gamma_table_path <- file.path(OUT_DIR, "trait_coefficients_table.csv")
write_csv(gamma_wide, gamma_table_path)
cat("Wrote", gamma_table_path, "\n")

# --- Side analysis: B1 vs B2 gamma stability --------------------------------
# Where both B1 and B2 exist (HG, DG Kuehne when ready), compute the mean
# difference and SD ratio per trait.

stability <- gamma_long %>%
  select(model, variant, variable, trait_label, mean, sd, q5, q95) %>%
  pivot_wider(names_from = variant, values_from = c(mean, sd, q5, q95),
              names_glue = "{.value}_{ifelse(variant == 'B1 species-free', 'B1', 'B2')}") %>%
  filter(!is.na(mean_B1), !is.na(mean_B2)) %>%
  mutate(
    mean_diff = round(mean_B1 - mean_B2, 4),
    sd_ratio  = round(sd_B1 / sd_B2, 3)
  ) %>%
  select(model, trait_label,
         mean_B1, mean_B2, mean_diff,
         sd_B1, sd_B2, sd_ratio,
         q5_B1, q95_B1, q5_B2, q95_B2)

if (nrow(stability) > 0) {
  stability_path <- file.path(OUT_DIR, "trait_coefficients_b1_b2_stability.csv")
  write_csv(stability, stability_path)
  cat("Wrote", stability_path, "\n")
  cat("\n=== B1 vs B2 gamma stability ===\n")
  print(stability, n = nrow(stability))
}

cat("\nDone.\n")
