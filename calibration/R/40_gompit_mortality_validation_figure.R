#!/usr/bin/env Rscript
# Gompit mortality validation figure (manuscript).
#
# Visualises the model-level comparison of Greg Johnson's gompit survival
# (hazard on crown ratio + crown closure at tree tip) against a per-species base
# survival rate, computed by compare_new_mortality.R on the held national FIA
# remeasurement panel (7.6M tree records, 133 species). No FVS engine is
# involved, so there is no growth / density-feedback confound: this is a pure
# prediction check of the mortality model.
#
# Inputs (from conus_mort/full_out/):
#   mortality_compare_overall.csv  n, obs/pred survival, AUC, logloss (new vs base)
#   mortality_compare_by_cch.csv   survival by crown-closure-at-tip quartile
#
# Output: gompit_mortality_validation.png (two panels)
#   A  survival by cch quartile: the base rate misses crowding mortality; gompit
#      tracks observed survival down into the most crowded quartile.
#   B  discrimination (AUC) and calibration (log loss), gompit vs base.
#
# Usage:
#   Rscript 40_gompit_mortality_validation_figure.R \
#     --indir conus_mort/full_out --out gompit_mortality_validation.png

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(patchwork)
}))

args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d = NULL) {
  i <- grep(paste0("^--", f, "="), args, value = TRUE)
  if (!length(i)) return(d)
  sub(paste0("^--", f, "="), "", i[1])
}
INDIR <- ga("indir", ".")
OUT   <- ga("out", "gompit_mortality_validation.png")

overall <- read_csv(file.path(INDIR, "mortality_compare_overall.csv"),
                    show_col_types = FALSE)
by_cch  <- read_csv(file.path(INDIR, "mortality_compare_by_cch.csv"),
                    show_col_types = FALSE)

# colour-blind-safe palette
col_obs  <- "#222222"
col_new  <- "#0072B2"   # gompit
col_base <- "#D55E00"   # base rate

# ---- Panel A: survival by crown-closure-at-tip quartile ----
qlev <- by_cch$cch_q
panelA_long <- by_cch %>%
  mutate(cch_q = factor(cch_q, levels = qlev)) %>%
  select(cch_q, Observed = obs_surv, Gompit = pred_new, `Base rate` = pred_base) %>%
  pivot_longer(-cch_q, names_to = "series", values_to = "surv") %>%
  mutate(series = factor(series, levels = c("Observed", "Gompit", "Base rate")))

pA <- ggplot(panelA_long, aes(cch_q, surv, colour = series, group = series)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(Observed = col_obs, Gompit = col_new,
                                 `Base rate` = col_base), name = NULL) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Crown closure at tree tip (quartile)",
       y = "Period survival probability",
       title = "A  Survival tracks crown closure",
       subtitle = "Base rate misses crowding mortality; gompit follows observed") +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.02, 0.05), legend.justification = c(0, 0),
        legend.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

# ---- Panel B: discrimination + calibration ----
panelB <- tibble(
  metric = c("AUC", "AUC", "Log loss", "Log loss"),
  model  = c("Gompit", "Base rate", "Gompit", "Base rate"),
  value  = c(overall$auc_new, overall$auc_base,
             overall$logloss_new, overall$logloss_base)
) %>% mutate(model = factor(model, levels = c("Gompit", "Base rate")))

pB <- ggplot(panelB, aes(model, value, fill = model)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.4, size = 3.6) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(Gompit = col_new, `Base rate` = col_base),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL,
       title = "B  Discrimination and calibration",
       subtitle = sprintf("National FIA remeasurement panel (n = %s, %d species)",
                           format(overall$n, big.mark = ","), 133L)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

fig <- pA + pB + plot_layout(widths = c(1.15, 1))
ggsave(OUT, fig, width = 11, height = 4.6, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
