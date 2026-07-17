#!/usr/bin/env Rscript
# fig_treemap_conus.R -- CONUS FVS x (TreeMap spatial vs FIADB uniform):
# the area-expansion choice washes out at CONUS, diverges at finer scales.
suppressMessages({library(data.table); library(ggplot2)})

SD  <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/treemap_conus"
d   <- fread(file.path(SD, "fvs_treemap_vs_fiadb.csv"))
d30 <- d[year == 2030 & is.finite(tm_over_fia)]

lab <- c(state = "State\n(48)", fortyp = "Forest type\n(143)")
box <- d30[scale %in% c("state", "fortyp")]
box[, sc := factor(lab[scale], levels = lab)]
conus_ratio <- d30[scale == "CONUS", tm_over_fia]

p <- ggplot(box, aes(sc, tm_over_fia)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_hline(yintercept = conus_ratio, colour = "#D55E00", linewidth = 0.7) +
  geom_boxplot(width = 0.5, outlier.size = 0.8, fill = "#0072B2", alpha = 0.25) +
  coord_cartesian(ylim = c(0.2, 2.0)) +
  labs(x = "aggregation scale",
       y = "TreeMap-spatial / FIADB-uniform carbon (2030)",
       title = "Area-expansion choice: negligible at CONUS, grows with resolution",
       subtitle = sprintf("FVS reserve live carbon; ratio 1 = agreement. CONUS = %.3f (orange line)",
                          conus_ratio)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(SD, "fvs_treemap_vs_fiadb_scales.png"), p,
       width = 7.2, height = 5, dpi = 200)
cat("wrote fvs_treemap_vs_fiadb_scales.png\n")
