# =============================================================================
# Title: Figure 1. Trait coefficient comparison, B1 species-free vs B2 species-aware
# Author: A. Weiskittel
# Date: 2026-05-13
# Description: Two-panel ggplot showing posterior mean and 90% CI for the eight
#   trait coefficients in the species-free (B1) and species-aware (B2) base
#   models for HG and DG Kuehne. Validates the trait-driven hypothesis across
#   two base equations. The DG Kue B1 has tighter CIs than HG B1 because the
#   100K subsample covers 62 percent of the DG dataset (n = 160,089) while
#   the 500K HG subsample covers only 25 percent of the HG dataset (n = 1.98M).
# Dependencies:
#   b1_b2_hg_gamma.csv     (from job 9333097 on Cardinal)
#   b1_b2_dg_kue_gamma.csv (from job 9470657 on Cardinal)
# Output:
#   figure1_trait_coefficients.pdf  (publication, vector, 17.5 x 12 cm)
#   figure1_trait_coefficients.png  (presentations, 300 dpi)
# =============================================================================

# --- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# --- Paths -------------------------------------------------------------------
hg_path     <- "b1_b2_hg_gamma.csv"
dg_kue_path <- "b1_b2_dg_kue_gamma.csv"
fig_pdf     <- "figure1_trait_coefficients.pdf"
fig_png     <- "figure1_trait_coefficients.png"

# --- Trait labels ------------------------------------------------------------
# Order in the W matrix in the cspi_traits1 Stan models (see v1 Methods 2.3)
trait_names <- c(
  "gamma[1]" = "Wood specific gravity",
  "gamma[2]" = "Shade tolerance",
  "gamma[3]" = "Softwood indicator",
  "gamma[4]" = "Leaf longevity",
  "gamma[5]" = "Maximum height",
  "gamma[6]" = "Maximum DBH",
  "gamma[7]" = "Vulnerability",
  "gamma[8]" = "Sensitivity"
)

# --- Data --------------------------------------------------------------------
# Parse the bracketed CI strings into numeric lower/upper bounds.
parse_ci <- function(ci_string) {
  m <- str_match(ci_string, "\\[\\s*([+-]?[0-9.]+),\\s*([+-]?[0-9.]+)\\s*\\]")
  tibble(lower = as.numeric(m[, 2]),
         upper = as.numeric(m[, 3]))
}

read_gamma <- function(path, base_label) {
  raw <- read_csv(path, show_col_types = FALSE)
  bind_cols(
    raw |> select(trait_idx, B1_mean, B2_mean, B1_excludes_zero, B2_excludes_zero),
    raw |> pull(B1_ci) |> parse_ci() |> rename(B1_lower = lower, B1_upper = upper),
    raw |> pull(B2_ci) |> parse_ci() |> rename(B2_lower = lower, B2_upper = upper)
  ) |>
    mutate(base_model = base_label,
           trait_label = factor(trait_names[trait_idx],
                                levels = unname(trait_names))) |>
    select(base_model, trait_idx, trait_label,
           B1_mean, B1_lower, B1_upper, B1_excludes_zero,
           B2_mean, B2_lower, B2_upper, B2_excludes_zero)
}

hg     <- read_gamma(hg_path,     "HG (height growth, ORGANON form)")
dg_kue <- read_gamma(dg_kue_path, "DG Kuehne (diameter growth)")

# Fix panel order: HG on the left to match manuscript Methods narrative flow
panel_levels <- c(
  "HG (height growth, ORGANON form)",
  "DG Kuehne (diameter growth)"
)

# Long form for plotting: one row per (trait, model)
long_df <- bind_rows(hg, dg_kue) |>
  mutate(base_model = factor(base_model, levels = panel_levels)) |>
  pivot_longer(cols = matches("^(B1|B2)_(mean|lower|upper|excludes_zero)$"),
               names_to = c("model", ".value"),
               names_pattern = "(B1|B2)_(.+)") |>
  rename(posterior_mean = mean,
         ci_lower = lower,
         ci_upper = upper,
         excludes_zero = excludes_zero) |>
  mutate(model_label = recode(model,
                              B1 = "B1 (species free)",
                              B2 = "B2 (species aware)"),
         model_label = factor(model_label,
                              levels = c("B2 (species aware)",
                                         "B1 (species free)")))

# Print the data we will plot, for the manuscript Section 3.7 numbers
cat("\n=== Long-form data going into Figure 1 ===\n")
print(long_df, n = Inf)

# Summary counts for the figure caption
caption_stats <- long_df |>
  group_by(base_model, model_label) |>
  summarise(excludes = sum(excludes_zero),
            total    = n(),
            .groups = "drop") |>
  mutate(label = sprintf("%d of %d trait coefficients exclude zero at 90%% CI",
                        excludes, total))

cat("\n=== Caption stats ===\n")
print(caption_stats)

# --- Theme -------------------------------------------------------------------
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_blank(),
    axis.title          = element_text(size = 11),
    axis.text           = element_text(size = 10),
    legend.position     = "bottom",
    legend.title        = element_text(size = 10),
    legend.text         = element_text(size = 9),
    strip.text          = element_text(size = 11, face = "bold"),
    plot.title          = element_text(size = 13, face = "bold"),
    plot.margin         = margin(8, 8, 8, 8)
  )

# Colors: B2 a muted grey (legacy), B1 a saturated orange (the new pivot)
model_colors <- c(
  "B2 (species aware)" = "#737373",
  "B1 (species free)"  = "#D55E00"
)

# --- Figure ------------------------------------------------------------------
# Dodge so B1 and B2 sit side by side per trait
pd <- position_dodge(width = 0.55)

p <- ggplot(long_df,
            aes(x = posterior_mean,
                y = trait_label,
                color = model_label,
                shape = excludes_zero)) +
  facet_wrap(~ base_model, scales = "free_x", ncol = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                width = 0.25, position = pd, linewidth = 0.55) +
  geom_point(size = 2.6, position = pd, fill = "white", stroke = 0.7) +
  scale_color_manual(values = model_colors, name = "Model") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "Excludes 0 at 90%",
                                `FALSE` = "Includes 0"),
                     name = "Posterior") +
  scale_y_discrete(limits = rev) +
  labs(x = "Trait coefficient γ (link scale, posterior mean and 90% CI)",
       y = NULL) +
  theme_pub +
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0.05, "cm"))

# --- Save --------------------------------------------------------------------
ggsave(fig_pdf, plot = p,
       width = 17.5, height = 11, units = "cm")
ggsave(fig_png, plot = p,
       width = 17.5, height = 11, units = "cm",
       dpi = 300, bg = "white")

cat(sprintf("\nWrote %s\nWrote %s\n", fig_pdf, fig_png))
