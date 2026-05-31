#!/usr/bin/env Rscript
# silc_strata_5x2_regen.R
# =====================================================================
# Re-aggregate SILC strata into a 5 forest-type by 2 density-class
# break (10 cells) per SILC's revised stratification request:
#   Forest type:  Cedar, Hardwood, Mixedwood, Commercial Softwood,
#                 Other Softwood
#   Density:      A+B (high), C+D (low)
#
# SILC stand-ID encoding (e.g. "S3B-N"):
#   C   = Cedar
#   H   = Hardwood
#   HS  = Hardwood-leading Mixedwood   } folded to Mixedwood
#   SH  = Softwood-leading Mixedwood   }
#   S   = Commercial Softwood
#   OS  = Other Softwood
#   3rd letter A/B/C/D = density class
#
# Data source for trajectories: GrownDB_byStrata_ALL.csv (AGM /
# AcadianGY output on the 11 SILC byStrata stands, year 2023->2123,
# 5 yr cycles). Multi-model trajectories from prior FVS-NE / FVS-ACD /
# OSM-ACD SILC runs are not in this local checkout -- they were in the
# silc_v25..v28 analysis tree that no longer mounts. This script
# produces the AGM figures on the new break and writes the
# stratification scaffold so the other-model data can be folded in
# with the same key once their projection outputs are restored.
# =====================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("Falling back to base R plotting (ggplot2 not available)\n")
    USE_GGPLOT <- FALSE
  } else {
    library(ggplot2)
    USE_GGPLOT <- TRUE
  }
})

base_dir <- "/sessions/friendly-compassionate-rubin/mnt/outputs/silc_extracted"
out_dir  <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"

# ----- 1. Strata mapping ---------------------------------------------
parse_stand_id <- function(stand_id) {
  # Returns list with forest_type, density_class, density_letter
  s <- toupper(as.character(stand_id))
  # extract the alphabetic prefix
  pre <- regmatches(s, regexpr("^[A-Z]+", s))
  # the density letter is the first A-D appearing after the prefix
  rest <- substring(s, nchar(pre) + 1)
  dens_letter <- regmatches(rest, regexpr("[A-D]", rest))
  if (length(dens_letter) == 0) dens_letter <- NA_character_

  ftype <- switch(pre,
    "C"  = "Cedar",
    "H"  = "Hardwood",
    "HS" = "Mixedwood",
    "SH" = "Mixedwood",
    "S"  = "Commercial Softwood",
    "OS" = "Other Softwood",
    NA_character_
  )
  dclass <- if (is.na(dens_letter)) NA_character_
            else if (dens_letter %in% c("A", "B")) "A+B (high)"
            else if (dens_letter %in% c("C", "D")) "C+D (low)"
            else NA_character_
  list(forest_type = ftype, density_class = dclass,
       density_letter = dens_letter, prefix = pre)
}

# Apply to all 79 matrix strata to produce the full reference table
all_strata <- read.csv(file.path(base_dir,
                                  "Acadian_Matrix_StandInit_2023.csv"),
                       stringsAsFactors = FALSE)
m <- t(sapply(all_strata$Stand, function(s) {
  p <- parse_stand_id(s)
  c(s, p$prefix, p$density_letter, p$forest_type, p$density_class)
}))
strata_map <- data.frame(
  stand_id      = m[, 1],
  prefix        = m[, 2],
  density_letter= m[, 3],
  forest_type   = m[, 4],
  density_class = m[, 5],
  num_plots     = all_strata$NUM_PLOTS,
  stringsAsFactors = FALSE
)

# Plot-count rollup on the new break
roll_plots <- aggregate(num_plots ~ forest_type + density_class,
                        data = strata_map, FUN = sum)
roll_plots <- roll_plots[order(roll_plots$forest_type,
                                roll_plots$density_class), ]
cat("=== Plot counts by new 5x2 stratification (n strata plots) ===\n")
print(roll_plots, row.names = FALSE)
write.csv(strata_map,
          file.path(out_dir, "silc_strata_5x2_mapping.csv"),
          row.names = FALSE)
write.csv(roll_plots,
          file.path(out_dir, "silc_strata_5x2_rollup.csv"),
          row.names = FALSE)

# ----- 2. AGM trajectories on the new break --------------------------
cat("\nReading GrownDB (this can take a minute)...\n")
grown <- read.csv(file.path(base_dir, "GrownDB_byStrata_ALL.csv"),
                  stringsAsFactors = FALSE)
cat(sprintf("GrownDB: %d rows, %d stands, %d years\n",
            nrow(grown), length(unique(grown$StandID)),
            length(unique(grown$Year))))

# Filter to commercial merch trees (DBH >= 4.6 in == BRK_DBH used
# in the StandInit). TPA in GrownDB is per-tree expansion.
g <- grown[grown$DBH >= 4.6 & is.finite(grown$DBH) & is.finite(grown$TPA), ]

# Stand-year aggregates (per acre)
# Force numeric and substitute zero for NA so they don't drop entire rows
for (col in c("TPA", "BA", "Cords", "NetCords")) {
  g[[col]] <- suppressWarnings(as.numeric(g[[col]]))
  g[[col]][!is.finite(g[[col]])] <- 0
}
stand_year <- aggregate(
  cbind(TPA, BA, Cords, NetCords) ~ StandID + Year,
  data = g, FUN = sum
)
# Use IT's QMD identity: QMD = sqrt(BA/TPA / 0.005454)  (inches)
stand_year$QMD <- with(stand_year, sqrt((BA / TPA) / 0.005454))

# Attach strata mapping
m2 <- do.call(rbind, lapply(stand_year$StandID, function(s) {
  p <- parse_stand_id(s)
  data.frame(forest_type   = p$forest_type,
             density_class = p$density_class,
             stringsAsFactors = FALSE)
}))
sy <- cbind(stand_year, m2)

# Trajectory: mean across stands in each (type, density) cell
traj <- aggregate(
  cbind(TPA, BA, QMD, Cords, NetCords) ~
    forest_type + density_class + Year,
  data = sy, FUN = mean
)
n_stands <- aggregate(StandID ~ forest_type + density_class,
                      data = sy[!duplicated(sy[, c("StandID",
                                                    "forest_type",
                                                    "density_class")]), ],
                      FUN = length)
names(n_stands)[3] <- "n_stands"
traj <- merge(traj, n_stands,
              by = c("forest_type", "density_class"))
traj <- traj[order(traj$forest_type, traj$density_class, traj$Year), ]
write.csv(traj,
          file.path(out_dir,
                    "silc_strata_5x2_AGM_trajectories.csv"),
          row.names = FALSE)
cat(sprintf("\nWrote AGM strata trajectories: %d cells, n_stands range %d-%d\n",
            length(unique(paste(traj$forest_type, traj$density_class))),
            min(n_stands$n_stands), max(n_stands$n_stands)))

# ----- 3. Figures ----------------------------------------------------
type_order <- c("Cedar", "Hardwood", "Mixedwood",
                "Commercial Softwood", "Other Softwood")
dens_order <- c("A+B (high)", "C+D (low)")
traj$forest_type   <- factor(traj$forest_type,   levels = type_order)
traj$density_class <- factor(traj$density_class, levels = dens_order)

# Build a label for each panel showing n stands (and "no data" if cell empty)
n_lookup <- setNames(n_stands$n_stands,
                     paste(n_stands$forest_type,
                           n_stands$density_class, sep = "|"))
panel_label <- function(ft, dc) {
  k <- paste(ft, dc, sep = "|")
  n <- n_lookup[k]
  if (is.na(n)) "no AGM data" else sprintf("n=%d stand%s", n,
                                            ifelse(n == 1, "", "s"))
}

# Figure helper using base R so we don't depend on ggplot
draw_metric <- function(metric, ylab, fname, ymin = NULL, ymax = NULL) {
  png(file.path(out_dir, fname),
      width = 2400, height = 1100, res = 165)
  par(mfrow = c(2, 5), mar = c(3.6, 3.8, 2.4, 0.8),
      mgp = c(2.3, 0.6, 0), oma = c(2.0, 2.5, 3.5, 0.5))

  cols_dens <- c("A+B (high)" = "#2c7fb8",
                 "C+D (low)"  = "#d95f02")

  # 2 rows (density), 5 cols (forest type)
  for (dc in dens_order) {
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
      if (is.null(ymin) || is.null(ymax)) {
        yr <- range(sub[[metric]], na.rm = TRUE)
        ymn <- yr[1]; ymx <- yr[2] * 1.05
      } else {
        ymn <- ymin; ymx <- ymax
      }
      plot(sub$Year, sub[[metric]], type = "l",
           col = cols_dens[dc], lwd = 2.5,
           ylim = c(ymn, ymx),
           xlab = "Year", ylab = ylab,
           main = paste(ft, "  (", panel_label(ft, dc), ")",
                        sep = ""),
           cex.main = 0.95, font.main = 2, las = 1)
      points(sub$Year, sub[[metric]], pch = 19,
             col = cols_dens[dc], cex = 0.7)
      grid(col = "#eeeeee", lty = 1)
    }
  }
  mtext(sprintf("AGM (AcadianGY) projection on SILC byStrata stands -- %s",
                ylab),
        outer = TRUE, line = 1.4, cex = 1.05, font = 2)
  mtext("Top row: A+B (high) density   |   Bottom row: C+D (low) density",
        outer = TRUE, line = 0.0, cex = 0.85, col = "#444")
  mtext("Left -> Right: Cedar, Hardwood, Mixedwood, Commercial Softwood, Other Softwood",
        outer = TRUE, side = 1, line = 0.5, cex = 0.82, col = "#444")
  dev.off()
  cat("wrote", file.path(out_dir, fname), "\n")
}

draw_metric("BA",       "Basal area (ft^2/ac)",
            "silc_strata_5x2_AGM_BA.png")
draw_metric("TPA",      "TPA (trees/ac, DBH >= 4.6 in)",
            "silc_strata_5x2_AGM_TPA.png")
draw_metric("QMD",      "QMD (in)",
            "silc_strata_5x2_AGM_QMD.png")
draw_metric("NetCords", "Net merchantable cords / ac",
            "silc_strata_5x2_AGM_NetCords.png")

# ----- 4. Stratification reference figure (plot counts by cell) ------
png(file.path(out_dir, "silc_strata_5x2_plot_counts.png"),
    width = 2000, height = 950, res = 175)
par(mar = c(5.5, 5.0, 4.0, 2.0), mgp = c(3.0, 0.8, 0))

# Build a 5x2 matrix of plot counts
M <- matrix(0, nrow = 2, ncol = 5,
            dimnames = list(dens_order, type_order))
for (i in seq_len(nrow(roll_plots))) {
  r <- roll_plots[i, ]
  if (!is.na(r$forest_type) && !is.na(r$density_class))
    M[r$density_class, r$forest_type] <- r$num_plots
}
bp <- barplot(M, beside = TRUE,
              col = c("#2c7fb8", "#d95f02"),
              border = NA,
              ylim = c(0, max(M) * 1.18),
              ylab = "Number of plots in SILC inventory",
              cex.names = 1.05, las = 1, font = 2)
for (j in seq_len(ncol(M))) {
  for (i in seq_len(nrow(M))) {
    v <- M[i, j]
    if (v > 0)
      text(bp[i, j], v + max(M) * 0.025,
           v, cex = 0.95, font = 2)
  }
}
legend("topright", legend = dens_order,
       fill = c("#2c7fb8", "#d95f02"),
       border = NA, bty = "n", cex = 1.0)
title(main = "SILC 5-type x 2-density stratification: plot counts",
      cex.main = 1.15, font.main = 2)
mtext("Density: A+B = high (BA index A/B); C+D = low (BA index C/D)",
      side = 1, line = 3.8, cex = 0.85, col = "#444")
dev.off()
cat("wrote silc_strata_5x2_plot_counts.png\n")

cat("\n=== Done. Artifacts in", out_dir, "===\n")
