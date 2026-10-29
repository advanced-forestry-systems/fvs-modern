#!/usr/bin/env Rscript
## paper2_fig2_bakuzis.R (2026-10-29)
## Paper 2 Figure 2: long-run realism of the four CONUS projector arms by variant.
##   (a) share of stands whose basal area exceeds the FIA observed 99.9th percentile by year 50
##   (b) share of stands whose peak SDI exceeds the FIA observed 99.9th percentile SDI by year 50
##   (c) within-stand Reineke slope (median, interquartile range) on RD-gated trajectories
## Inputs: the catalog gate table (projector rows) and the four by-trajectory scorecards from the October 27 rescore.
## Output: paper2_fig2_bakuzis.png, 180 x 200 mm, 300 dpi; paper2_fig2_data.csv.
suppressMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE)
ga <- function(k, d = NA) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v else d }
G <- fread(ga("gates")); TRAJ <- ga("traj"); OUT <- ga("out", "paper2_fig2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
p <- G[grepl("^conus", arm) & !is.na(ba_p999)]
p[, arm_lab := factor(arm, levels = c("conus_b1", "conus_b1_gompit", "conus_b2", "conus_b2_gompit"),
                      labels = c("b1 native", "b1 GOMPIT", "b2 native", "b2 GOMPIT"))]
p[, survival := ifelse(grepl("gompit", arm), "GOMPIT", "native")]; p[, growth := ifelse(grepl("b1", arm), "b1", "b2")]
ord <- p[arm == "conus_b1"][order(-share_ba999_50), group]; p[, group := factor(group, levels = ord)]
T <- rbindlist(lapply(strsplit(TRAJ, ",")[[1]], fread), fill = TRUE)
T <- T[grepl("^T3 Reineke", relation) & status %in% c("PASS", "FLAG") & grepl("^conus", arm)]
T[, value := as.numeric(value)]; T[, group := toupper(group)]
r <- T[, .(slope_med = median(value), q25 = quantile(value, .25), q75 = quantile(value, .75), n = .N), by = .(arm, group)]
r[, arm_lab := factor(arm, levels = c("conus_b1", "conus_b1_gompit", "conus_b2", "conus_b2_gompit"),
                      labels = c("b1 native", "b1 GOMPIT", "b2 native", "b2 GOMPIT"))]
r <- r[group %in% ord]; r[, group := factor(group, levels = ord)]
th <- theme_bw(base_size = 9) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                       legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 0))
cols <- c("b1 native" = "#9e2a2b", "b1 GOMPIT" = "#e09f3e", "b2 native" = "#335c67", "b2 GOMPIT" = "#540b0e")
pa <- ggplot(p, aes(group, share_ba999_50, fill = arm_lab)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0.05, linetype = 2, linewidth = 0.3) + scale_fill_manual(values = cols, name = "Projector arm") +
  scale_y_continuous(labels = function(x) x * 100) + labs(x = NULL, y = expression("Stands above BA " * p[99.9] * " (%)")) + th
pb <- ggplot(p, aes(group, share_sdi_50, fill = arm_lab)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0.05, linetype = 2, linewidth = 0.3) + scale_fill_manual(values = cols, name = "Projector arm") +
  scale_y_continuous(labels = function(x) x * 100) + labs(x = NULL, y = expression("Stands above SDI " * p[99.9] * " (%)")) + th
pc <- ggplot(r, aes(group, slope_med, colour = arm_lab)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2.2, ymax = -1.2, alpha = 0.12, fill = "grey40") +
  geom_hline(yintercept = -1.605, linetype = 3, linewidth = 0.3) +
  geom_pointrange(aes(ymin = q25, ymax = q75), position = position_dodge(width = 0.7), size = 0.25, linewidth = 0.4) + scale_colour_manual(values = cols, name = "Projector arm", guide = "none") +
  labs(x = "FVS variant", y = "Reineke slope (median, IQR)") + th
png(file.path(OUT, "paper2_fig2_bakuzis.png"), width = 180, height = 200, units = "mm", res = 300)
if (requireNamespace("patchwork", quietly = TRUE)) { library(patchwork); print((pa / pb / pc) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")) } else {
  if (requireNamespace("gridExtra", quietly = TRUE)) gridExtra::grid.arrange(pa, pb, pc, ncol = 1) else { print(pa) } }
dev.off()
fwrite(merge(p[, .(arm, group, n_traj, share_ba999_50, share_ba999_100, share_sdi_50, share_sdi_100, reineke_flag_rate, reineke_eval)], r[, .(arm, group, slope_med, q25, q75, n_slope = n)], by = c("arm", "group"), all.x = TRUE), file.path(OUT, "paper2_fig2_data.csv"))
cat("wrote", file.path(OUT, "paper2_fig2_bakuzis.png"), "\n")
