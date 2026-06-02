#!/usr/bin/env Rscript
# silc_strata_5x2_year100.R
# =====================================================================
# Year-100 outcomes table for the SILC 5-type x 2-density break,
# using AGM (AcadianGY) byStrata trajectories. Reports per-cell year-0
# and year-100 BA / TPA / QMD / NetCords with 100-yr growth factor
# and n_stands.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"
traj <- read.csv(file.path(od, "silc_strata_5x2_AGM_trajectories.csv"),
                 stringsAsFactors = FALSE)

# Earliest and latest years per cell
out <- do.call(rbind, lapply(unique(paste(traj$forest_type, traj$density_class)),
  function(k) {
    parts <- unlist(strsplit(k, " ", fixed=TRUE))
    # forest_type can have a space (Commercial Softwood / Other Softwood)
    n_words <- length(parts)
    dc <- paste(parts[(n_words-1):n_words], collapse=" ")
    ft <- paste(parts[seq_len(n_words-2)], collapse=" ")
    sub <- traj[traj$forest_type == ft & traj$density_class == dc, ]
    sub <- sub[order(sub$Year), ]
    y0  <- sub[1, ]
    y100<- sub[nrow(sub), ]
    data.frame(
      forest_type   = ft,
      density_class = dc,
      year_start    = y0$Year,
      year_end      = y100$Year,
      BA_0    = round(y0$BA, 1),
      BA_100  = round(y100$BA, 1),
      TPA_0   = round(y0$TPA, 0),
      TPA_100 = round(y100$TPA, 0),
      QMD_0   = round(y0$QMD, 2),
      QMD_100 = round(y100$QMD, 2),
      NetCords_0  = round(y0$NetCords, 1),
      NetCords_100= round(y100$NetCords, 1),
      growth_factor_BA      = round(y100$BA / y0$BA, 2),
      growth_factor_NetCords= round(y100$NetCords / y0$NetCords, 2),
      n_stands = y100$n_stands,
      stringsAsFactors = FALSE
    )
  }))

# Add empty cells with NA
type_order <- c("Cedar", "Hardwood", "Mixedwood",
                "Commercial Softwood", "Other Softwood")
dens_order <- c("A+B (high)", "C+D (low)")
full_grid <- expand.grid(forest_type = type_order,
                          density_class = dens_order,
                          stringsAsFactors = FALSE)
out <- merge(full_grid, out,
             by = c("forest_type","density_class"), all.x = TRUE)
out$forest_type   <- factor(out$forest_type, levels = type_order)
out$density_class <- factor(out$density_class, levels = dens_order)
out <- out[order(out$forest_type, out$density_class), ]

# Mark empty cells
out$n_stands[is.na(out$n_stands)] <- 0

write.csv(out, file.path(od, "silc_strata_5x2_year100_outcomes.csv"),
          row.names = FALSE)

cat("=== Year-100 outcomes by 5x2 cell (AGM only) ===\n")
print(out[, c("forest_type","density_class","n_stands",
              "BA_0","BA_100","growth_factor_BA",
              "NetCords_0","NetCords_100","growth_factor_NetCords")],
      row.names = FALSE)

# Figure: BA growth factor + NetCords growth factor side by side
png(file.path(od, "silc_strata_5x2_year100_growth.png"),
    width = 2400, height = 950, res = 165)
par(mfrow = c(1, 2), mar = c(5.5, 5.0, 3.5, 1.0), mgp = c(2.8, 0.6, 0))

type_abbr <- c("Cedar" = "Cedar", "Hardwood" = "Hardwood",
                "Mixedwood" = "Mixedwood",
                "Commercial Softwood" = "Comm. SW",
                "Other Softwood" = "Other SW")
draw_panel <- function(metric, ylab, title_main) {
  M <- matrix(0, nrow = 2, ncol = 5,
              dimnames = list(dens_order, type_abbr[type_order]))
  for (i in seq_len(nrow(out))) {
    v <- out[[metric]][i]
    if (!is.na(v))
      M[as.character(out$density_class[i]),
        type_abbr[as.character(out$forest_type[i])]] <- v
  }
  bp <- barplot(M, beside = TRUE,
                col = c("#2c7fb8", "#d95f02"),
                border = NA,
                ylim = c(0, max(M, na.rm=TRUE) * 1.25),
                ylab = ylab, las = 1, cex.names = 1.0, font = 2)
  for (j in seq_len(ncol(M))) {
    for (i in seq_len(nrow(M))) {
      v <- M[i, j]
      if (v > 0)
        text(bp[i, j], v + max(M)*0.03, sprintf("%.2f", v),
             cex = 0.85, font = 2)
      else
        text(bp[i, j], max(M)*0.05, "n/a",
             cex = 0.8, col = "#aa3333", font = 3)
    }
  }
  legend("top", legend = dens_order, horiz = TRUE,
         fill = c("#2c7fb8", "#d95f02"),
         border = NA, bty = "n", cex = 0.95, inset = -0.02)
  title(main = title_main, cex.main = 1.1, font.main = 2)
  mtext("Growth factor = year-100 value / year-0 value",
        side = 1, line = 4.0, cex = 0.78, col = "#444")
}

draw_panel("growth_factor_BA",
           "BA growth factor (year100 / year0)",
           "BA growth factor by 5x2 cell")
draw_panel("growth_factor_NetCords",
           "Net merch cords growth factor",
           "Merchantable cords growth factor by 5x2 cell")

dev.off()
cat("\nwrote silc_strata_5x2_year100_outcomes.csv\n")
cat("wrote silc_strata_5x2_year100_growth.png\n")
