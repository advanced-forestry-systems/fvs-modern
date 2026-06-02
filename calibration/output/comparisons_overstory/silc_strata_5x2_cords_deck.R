#!/usr/bin/env Rscript
# silc_strata_5x2_cords_deck.R
# =====================================================================
# Deck-ready merch-cords trajectory figure: one panel per 5x2 cell,
# y-axis fixed across the grid so cross-stratum yield potential is
# visually comparable. CRSF brand colours.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"
traj <- read.csv(file.path(od, "silc_strata_5x2_AGM_trajectories.csv"),
                 stringsAsFactors = FALSE)
type_order <- c("Cedar","Hardwood","Mixedwood",
                 "Commercial Softwood","Other Softwood")
dens_order <- c("A+B (high)","C+D (low)")
traj$forest_type <- factor(traj$forest_type, levels = type_order)
traj$density_class <- factor(traj$density_class, levels = dens_order)

n_lookup <- setNames(unique(traj[, c("forest_type","density_class","n_stands")])$n_stands,
                     paste(unique(traj[, c("forest_type","density_class","n_stands")])$forest_type,
                            unique(traj[, c("forest_type","density_class","n_stands")])$density_class,
                            sep="|"))

ymax_g <- max(traj$NetCords, na.rm = TRUE) * 1.05

png(file.path(od, "silc_strata_5x2_cords_deck.png"),
    width = 2700, height = 1100, res = 170)
par(mfrow = c(2, 5), mar = c(3.6, 4.2, 2.6, 0.8),
    mgp = c(2.4, 0.6, 0), oma = c(2.0, 2.5, 3.5, 0.5))

CRSF_GREEN <- "#1A3D28"; CRSF_ACCENT <- "#88A47A"

for (dc in dens_order) {
  for (ft in type_order) {
    sub <- traj[traj$forest_type == ft & traj$density_class == dc, ]
    if (nrow(sub) == 0) {
      plot.new(); box(col = "#cccccc", lty = 2)
      text(0.5, 0.55, ft, cex = 1.1, font = 2)
      text(0.5, 0.42, dc, cex = 0.9, col = "#666")
      text(0.5, 0.25, "no AGM data", cex = 0.95, col = "#aa3333", font = 3)
      next
    }
    col_line <- if (dc == "A+B (high)") CRSF_GREEN else CRSF_ACCENT
    plot(sub$Year, sub$NetCords, type = "l", col = col_line, lwd = 3,
         ylim = c(0, ymax_g),
         xlab = "Year", ylab = "Net merch cords/ac",
         main = sprintf("%s   (n=%d)", ft, n_lookup[paste(ft, dc, sep="|")]),
         cex.main = 1.0, font.main = 2, las = 1)
    points(sub$Year, sub$NetCords, pch = 19, col = col_line, cex = 0.8)
    grid(col = "#eeeeee", lty = 1)
    # year-100 value annotation
    yend <- sub$NetCords[nrow(sub)]
    abline(h = yend, col = col_line, lty = 3)
    text(sub$Year[1] + 5, yend, sprintf("%.0f cords", yend),
         pos = 3, col = col_line, cex = 0.85, font = 2)
  }
}
mtext("Top row: A+B (high) density          Bottom row: C+D (low) density",
      outer = TRUE, side = 3, line = 0.5, cex = 0.95, col = "#444")
dev.off()
cat("wrote silc_strata_5x2_cords_deck.png\n")
