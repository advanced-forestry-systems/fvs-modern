#!/usr/bin/env Rscript
# fig_strata_trends.R -- landowner carbon trajectories, reserve vs managed,
# with bootstrap uncertainty ribbons (default engine).
suppressMessages({library(data.table); library(ggplot2)})
SD <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/strata_trends"
d  <- fread(file.path(SD, "strata_trends_default.csv"))
o  <- d[scale == "owner" & key != "Unknown"]
o[, key := factor(key, levels = c("Industrial","NIPF","State","Public-Other"))]

p <- ggplot(o, aes(year, total_TgC, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = total_lo, ymax = total_hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ key, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("reserve (no harvest)" = "#0072B2",
                                 "managed (harvest)" = "#D55E00"),
                      aesthetics = c("colour","fill")) +
  labs(x = NULL, y = "live aboveground carbon (Tg C)",
       title = "CONUS forest carbon by landowner: reserve vs managed",
       subtitle = "FVS default engine; ribbons = bootstrap 95% CI (plot resampling)",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold"))
ggsave(file.path(SD, "fvs_owner_trends.png"), p, width = 8, height = 6, dpi = 200)
cat("wrote fvs_owner_trends.png\n")
