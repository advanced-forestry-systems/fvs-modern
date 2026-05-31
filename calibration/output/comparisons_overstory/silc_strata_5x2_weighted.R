#!/usr/bin/env Rscript
# silc_strata_5x2_weighted.R
# =====================================================================
# Plot-count weighted AGM strata trajectories with uncertainty bands.
#
# Two upgrades vs the unweighted version:
#  1. Each byStrata stand contributes proportionally to its underlying
#     SILC plot count (NUM_PLOTS from StandInit), so a stand backing
#     1063 Hardwood C+D plots is weighted 13x more than the stand
#     backing 28 Cedar A+B plots.
#  2. Between-stand 5th-95th percentile ribbons within each cell.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"

# AGM byStrata trajectory data (per-stand)
gdb <- read.csv("/sessions/friendly-compassionate-rubin/mnt/outputs/silc_extracted/GrownDB_byStrata_ALL.csv")
si  <- read.csv("/sessions/friendly-compassionate-rubin/mnt/outputs/silc_extracted/Acadian_Matrix_StandInit_2023.csv")
mapping <- read.csv(file.path(od, "silc_strata_5x2_mapping.csv"))

# Per-stand year totals (DBH >= 4.6 in)
g <- gdb[gdb$DBH >= 4.6 & is.finite(gdb$DBH), ]
for (col in c("TPA","BA","Cords","NetCords")) {
  g[[col]] <- suppressWarnings(as.numeric(g[[col]]))
  g[[col]][!is.finite(g[[col]])] <- 0
}
sy <- aggregate(cbind(TPA, BA, Cords, NetCords) ~ StandID + Year,
                data = g, FUN = sum)
sy$QMD <- sqrt((sy$BA / sy$TPA) / 0.005454)

# Attach strata + per-stand plot counts
si$plot_count <- si$NUM_PLOTS
stand_plots <- setNames(si$plot_count, si$Stand)
sy$num_plots <- stand_plots[sy$StandID]
m <- merge(sy, mapping[, c("stand_id","forest_type","density_class")],
           by.x = "StandID", by.y = "stand_id", all.x = TRUE)

# Weighted mean + percentile band per (cell, year)
type_order <- c("Cedar","Hardwood","Mixedwood",
                 "Commercial Softwood","Other Softwood")
dens_order <- c("A+B (high)","C+D (low)")

build_cell_traj <- function(metric) {
  out <- data.frame()
  for (ft in type_order) {
    for (dc in dens_order) {
      sub <- m[m$forest_type == ft & m$density_class == dc, ]
      if (nrow(sub) == 0) next
      for (y in sort(unique(sub$Year))) {
        sy <- sub[sub$Year == y, ]
        w  <- sy$num_plots
        v  <- sy[[metric]]
        mu <- sum(w * v) / sum(w)
        # weighted 5-95 percentile: convert to weighted ECDF
        ord <- order(v); v <- v[ord]; w <- w[ord]
        cw <- cumsum(w) / sum(w)
        q05 <- v[which(cw >= 0.05)[1]]
        q95 <- v[which(cw >= 0.95)[1]]
        if (length(q05) == 0 || is.na(q05)) q05 <- min(v)
        if (length(q95) == 0 || is.na(q95)) q95 <- max(v)
        out <- rbind(out, data.frame(
          forest_type = ft, density_class = dc, Year = y,
          mean = mu, q05 = q05, q95 = q95, n_stands = nrow(sy),
          total_plots = sum(w)
        ))
      }
    }
  }
  out
}

ba_cell    <- build_cell_traj("BA")
cords_cell <- build_cell_traj("NetCords")

write.csv(ba_cell,    file.path(od, "silc_strata_5x2_BA_weighted.csv"),
          row.names = FALSE)
write.csv(cords_cell, file.path(od, "silc_strata_5x2_NetCords_weighted.csv"),
          row.names = FALSE)

# === Figure: 5x2 grid with uncertainty ribbon ===
draw_with_ribbon <- function(d, ylab, fname) {
  png(file.path(od, fname),
      width = 2700, height = 1100, res = 170)
  par(mfrow = c(2, 5), mar = c(3.6, 4.2, 2.6, 0.8),
      mgp = c(2.4, 0.6, 0), oma = c(2.0, 2.5, 3.0, 0.5))
  CRSF_GREEN <- "#1A3D28"; CRSF_ACCENT <- "#88A47A"
  ymax_g <- max(d$q95, na.rm = TRUE) * 1.05
  for (dc in dens_order) {
    for (ft in type_order) {
      sub <- d[d$forest_type == ft & d$density_class == dc, ]
      if (nrow(sub) == 0) {
        plot.new(); box(col = "#cccccc", lty = 2)
        text(0.5, 0.55, ft, cex = 1.0, font = 2)
        text(0.5, 0.42, dc, cex = 0.85, col = "#666")
        text(0.5, 0.25, "no AGM data", cex = 0.9, col = "#aa3333", font = 3)
        next
      }
      col_line <- if (dc == "A+B (high)") CRSF_GREEN else CRSF_ACCENT
      col_ribbon <- adjustcolor(col_line, alpha.f = 0.20)
      plot(NA, xlim = range(sub$Year), ylim = c(0, ymax_g),
           xlab = "Year", ylab = ylab,
           main = sprintf("%s   (n=%d stands, %d plots)",
                          ft, sub$n_stands[1], sub$total_plots[1]),
           cex.main = 0.95, font.main = 2, las = 1)
      polygon(c(sub$Year, rev(sub$Year)), c(sub$q05, rev(sub$q95)),
              col = col_ribbon, border = NA)
      lines(sub$Year, sub$mean, col = col_line, lwd = 3)
      points(sub$Year, sub$mean, pch = 19, col = col_line, cex = 0.7)
      grid(col = "#eeeeee", lty = 1)
    }
  }
  mtext("Plot-count weighted mean (line) and between-stand 5 to 95 percentile band (shaded)",
        outer = TRUE, side = 3, line = 0.5, cex = 0.95, col = "#444")
  dev.off()
  cat("wrote", fname, "\n")
}

draw_with_ribbon(cords_cell, "Net merch cords/ac",
                 "silc_strata_5x2_cords_weighted.png")
draw_with_ribbon(ba_cell,    "Basal area (ft^2/ac)",
                 "silc_strata_5x2_BA_weighted.png")
cat("\ndone\n")
