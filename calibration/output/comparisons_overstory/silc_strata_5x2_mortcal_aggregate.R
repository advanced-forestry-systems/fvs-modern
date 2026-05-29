#!/usr/bin/env Rscript
# silc_strata_5x2_mortcal_aggregate.R
# =====================================================================
# Re-aggregate the AcadianGY MORTCAL=TRUE 100-yr trajectories from
# Cardinal into the 5 forest-type x 2 density-class break, then
# compute side-by-side AGM default vs MORTCAL trajectory tables and
# year-100 outcomes.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"

mc <- read.csv(file.path(od, "silc_strata_100yr_mortcal_trajectories.csv"),
               stringsAsFactors=FALSE)
def<- read.csv(file.path(od, "silc_strata_5x2_AGM_trajectories.csv"),
               stringsAsFactors=FALSE)
sm <- read.csv(file.path(od, "silc_strata_5x2_mapping.csv"),
               stringsAsFactors=FALSE)

# Attach strata mapping to MORTCAL output
mc$forest_type   <- sm$forest_type[match(mc$StandID, sm$stand_id)]
mc$density_class <- sm$density_class[match(mc$StandID, sm$stand_id)]

# Aggregate to (forest_type x density_class x Year) means
mc_traj <- aggregate(
  cbind(TPA, BA, QMD, Cords, NetCords, BdFt) ~
    forest_type + density_class + Year,
  data = mc, FUN = mean
)
n_stands_mc <- aggregate(
  StandID ~ forest_type + density_class,
  data = mc[!duplicated(mc[, c("StandID","forest_type","density_class")]), ],
  FUN = length)
names(n_stands_mc)[3] <- "n_stands"
mc_traj <- merge(mc_traj, n_stands_mc, by = c("forest_type","density_class"))

write.csv(mc_traj, file.path(od, "silc_strata_5x2_AGM_MORTCAL_trajectories.csv"),
          row.names=FALSE)
cat(sprintf("Wrote MORTCAL strata trajectories: %d cells, n_stands range %d-%d\n",
            length(unique(paste(mc_traj$forest_type, mc_traj$density_class))),
            min(n_stands_mc$n_stands), max(n_stands_mc$n_stands)))

# Year-100 outcomes (default vs MORTCAL side by side)
def100 <- def[def$Year == 2123, c("forest_type","density_class","TPA","BA","QMD","Cords","NetCords","n_stands")]
names(def100) <- c("forest_type","density_class",
                    "TPA_def","BA_def","QMD_def","Cords_def","NetCords_def","n_def")
mc100 <- mc_traj[mc_traj$Year == 2123,
                  c("forest_type","density_class","TPA","BA","QMD","Cords","NetCords","BdFt","n_stands")]
names(mc100) <- c("forest_type","density_class",
                   "TPA_mc","BA_mc","QMD_mc","Cords_mc","NetCords_mc","BdFt_mc","n_mc")
y100 <- merge(def100, mc100, by=c("forest_type","density_class"), all=TRUE)
y100$dBA_pct    <- 100 * (y100$BA_mc / y100$BA_def - 1)
# Default AGM trajectory has Cords col stuck at 0; NetCords is the populated one
y100$dCords_pct <- 100 * (y100$NetCords_mc / y100$NetCords_def - 1)
write.csv(y100, file.path(od, "silc_strata_5x2_year100_mortcal_vs_default.csv"),
          row.names=FALSE)
cat("\n=== Year 2123 outcomes: AGM default vs MORTCAL ===\n")
print(y100[order(y100$forest_type, y100$density_class),
            c("forest_type","density_class","BA_def","BA_mc","dBA_pct",
              "NetCords_def","NetCords_mc","dCords_pct","n_mc")],
      row.names=FALSE, digits=3)

# === Figures: trajectory side-by-side per cell ===
CRSF_GREEN <- "#1A3D28"; MORT_GOLD <- "#B8860B"
type_order <- c("Cedar","Hardwood","Mixedwood","Commercial Softwood","Other Softwood")
dens_order <- c("A+B (high)","C+D (low)")

png(file.path(od, "silc_strata_5x2_AGM_MORTCAL_BA.png"),
    width=2400, height=1400, res=170)
par(mfrow=c(2,5), mar=c(3.5,4,3,1), mgp=c(2.4,0.6,0), oma=c(0,0,2.5,0))
for (dc in dens_order) {
  for (ft in type_order) {
    d_def <- def[def$forest_type==ft & def$density_class==dc, ]
    d_mc  <- mc_traj[mc_traj$forest_type==ft & mc_traj$density_class==dc, ]
    if (nrow(d_def) == 0 && nrow(d_mc) == 0) {
      plot.new(); title(main=sprintf("%s\n%s", ft, dc), cex.main=0.9)
      text(0.5, 0.5, "no data", col="grey50")
      next
    }
    ymax <- max(c(d_def$BA, d_mc$BA), na.rm=TRUE) * 1.1
    plot(d_def$Year, d_def$BA, type="l", lwd=2.5, col=CRSF_GREEN,
         ylim=c(0, ymax), xlab="", ylab="BA (ft^2/ac)",
         main=sprintf("%s\n%s (n=%d)", ft, dc, max(d_mc$n_stands, 0)),
         cex.main=0.9, las=1)
    if (nrow(d_mc) > 0)
      lines(d_mc$Year, d_mc$BA, lwd=2.5, col=MORT_GOLD)
    grid(col="grey90")
  }
}
mtext("100-yr AGM trajectories: default (green) vs MORTCAL=TRUE (gold)",
      side=3, line=0, outer=TRUE, font=2, cex=1.1)
dev.off()
cat("\nwrote silc_strata_5x2_AGM_MORTCAL_BA.png\n")

png(file.path(od, "silc_strata_5x2_AGM_MORTCAL_Cords.png"),
    width=2400, height=1400, res=170)
par(mfrow=c(2,5), mar=c(3.5,4,3,1), mgp=c(2.4,0.6,0), oma=c(0,0,2.5,0))
for (dc in dens_order) {
  for (ft in type_order) {
    d_def <- def[def$forest_type==ft & def$density_class==dc, ]
    d_mc  <- mc_traj[mc_traj$forest_type==ft & mc_traj$density_class==dc, ]
    if (nrow(d_def) == 0 && nrow(d_mc) == 0) {
      plot.new(); title(main=sprintf("%s\n%s", ft, dc), cex.main=0.9)
      text(0.5, 0.5, "no data", col="grey50")
      next
    }
    ymax <- max(c(d_def$NetCords, d_mc$NetCords), na.rm=TRUE) * 1.15
    plot(d_def$Year, d_def$NetCords, type="l", lwd=2.5, col=CRSF_GREEN,
         ylim=c(0, ymax), xlab="", ylab="Net Cords/ac",
         main=sprintf("%s\n%s (n=%d)", ft, dc, max(d_mc$n_stands, 0)),
         cex.main=0.9, las=1)
    if (nrow(d_mc) > 0)
      lines(d_mc$Year, d_mc$NetCords, lwd=2.5, col=MORT_GOLD)
    grid(col="grey90")
  }
}
mtext("100-yr AGM NetCords trajectories: default (green) vs MORTCAL=TRUE (gold)",
      side=3, line=0, outer=TRUE, font=2, cex=1.1)
dev.off()
cat("wrote silc_strata_5x2_AGM_MORTCAL_Cords.png\n")
