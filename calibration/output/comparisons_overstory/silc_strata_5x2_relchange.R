#!/usr/bin/env Rscript
# silc_strata_5x2_relchange.R
# =====================================================================
# Relative-change versions of the AGM 5x2 trajectory figures, plus
# uniform y-axis per density row so the cross-forest-type comparison
# is visually readable within each row.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"
traj <- read.csv(file.path(od, "silc_strata_5x2_AGM_trajectories.csv"),
                 stringsAsFactors = FALSE)
type_order <- c("Cedar","Hardwood","Mixedwood",
                 "Commercial Softwood","Other Softwood")
dens_order <- c("A+B (high)","C+D (low)")
traj$forest_type   <- factor(traj$forest_type, levels = type_order)
traj$density_class <- factor(traj$density_class, levels = dens_order)

# Relative-to-year-0 transform per cell
traj$rel_BA <- NA; traj$rel_TPA <- NA
traj$rel_QMD <- NA; traj$rel_NetCords <- NA
for (k in unique(paste(traj$forest_type, traj$density_class))) {
  parts <- unlist(strsplit(k, " ", fixed=TRUE))
  n <- length(parts)
  dc <- paste(parts[(n-1):n], collapse=" ")
  ft <- paste(parts[seq_len(n-2)], collapse=" ")
  idx <- traj$forest_type == ft & traj$density_class == dc
  sub <- traj[idx, ]
  sub <- sub[order(sub$Year), ]
  traj$rel_BA[idx]    <- sub$BA       / sub$BA[1]
  traj$rel_TPA[idx]   <- sub$TPA      / sub$TPA[1]
  traj$rel_QMD[idx]   <- sub$QMD      / sub$QMD[1]
  traj$rel_NetCords[idx] <- sub$NetCords / sub$NetCords[1]
}

# n_stands lookup
n_stands <- unique(traj[, c("forest_type","density_class","n_stands")])
n_lookup <- setNames(n_stands$n_stands,
                     paste(n_stands$forest_type, n_stands$density_class, sep="|"))
panel_label <- function(ft, dc) {
  k <- paste(ft, dc, sep = "|")
  n <- n_lookup[k]
  if (is.na(n)) "no AGM data" else
    sprintf("n=%d stand%s", n, ifelse(n == 1, "", "s"))
}

draw <- function(metric, ylab, fname, fixed_ymin = NULL, fixed_ymax = NULL,
                 horiz_line_at = NULL) {
  png(file.path(od, fname),
      width = 2400, height = 1100, res = 165)
  par(mfrow = c(2, 5), mar = c(3.6, 3.8, 2.4, 0.8),
      mgp = c(2.3, 0.6, 0), oma = c(2.0, 2.5, 3.5, 0.5))

  cols_dens <- c("A+B (high)" = "#2c7fb8",
                 "C+D (low)"  = "#d95f02")

  for (dc in dens_order) {
    # row y-range -- max across cells in this row, fixed across the row
    row_vals <- traj[traj$density_class == dc, metric]
    if (is.null(fixed_ymin)) ymn <- min(row_vals, na.rm = TRUE) * 0.97
    else ymn <- fixed_ymin
    if (is.null(fixed_ymax)) ymx <- max(row_vals, na.rm = TRUE) * 1.05
    else ymx <- fixed_ymax

    for (ft in type_order) {
      sub <- traj[traj$forest_type == ft &
                  traj$density_class == dc, ]
      if (nrow(sub) == 0) {
        plot.new()
        box(col = "#cccccc", lty = 2)
        text(0.5, 0.55, ft, cex = 1.0, font = 2)
        text(0.5, 0.42, dc, cex = 0.85, col = "#666")
        text(0.5, 0.25, "no AGM data", cex = 0.9,
             col = "#aa3333", font = 3)
        next
      }
      plot(sub$Year, sub[[metric]], type = "l",
           col = cols_dens[dc], lwd = 2.5,
           ylim = c(ymn, ymx),
           xlab = "Year", ylab = ylab,
           main = paste(ft, "  (", panel_label(ft, dc), ")", sep = ""),
           cex.main = 0.95, font.main = 2, las = 1)
      points(sub$Year, sub[[metric]], pch = 19,
             col = cols_dens[dc], cex = 0.7)
      grid(col = "#eeeeee", lty = 1)
      if (!is.null(horiz_line_at))
        abline(h = horiz_line_at, col = "#888", lty = 2)
    }
  }
  mtext(sprintf("AGM (AcadianGY) projection on SILC byStrata stands -- %s",
                ylab),
        outer = TRUE, line = 1.4, cex = 1.05, font = 2)
  mtext("Top row: A+B (high) density   |   Bottom row: C+D (low) density   |   y-axis fixed within each row",
        outer = TRUE, line = 0.0, cex = 0.85, col = "#444")
  mtext("Left -> Right: Cedar, Hardwood, Mixedwood, Commercial Softwood, Other Softwood",
        outer = TRUE, side = 1, line = 0.5, cex = 0.82, col = "#444")
  dev.off()
  cat("wrote", fname, "\n")
}

# Relative-change figures: y in growth-factor units
draw("rel_BA",    "BA / BA(year 0)",
     "silc_strata_5x2_AGM_relBA.png",
     horiz_line_at = 1.0)
draw("rel_NetCords","Net merch cords / cords(year 0)",
     "silc_strata_5x2_AGM_relNetCords.png",
     horiz_line_at = 1.0)

# Absolute-units versions with row-uniform y-axes
draw("BA",       "Basal area (ft^2/ac)",
     "silc_strata_5x2_AGM_BA_unifrow.png")
draw("NetCords", "Net merchantable cords / ac",
     "silc_strata_5x2_AGM_NetCords_unifrow.png")

cat("\nDone.\n")
