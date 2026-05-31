#!/usr/bin/env Rscript
# fia_overstory_perplot_dist.R   (base R)
# =====================================================================
# Per-plot bias distribution for the apples-to-apples FIA overstory
# recompute. Complements the bias-of-means scorecard with the spread of
# individual-plot errors, so the "essentially unbiased" headline reads
# correctly as landscape-mean, not per-stand.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"
d  <- read.csv(file.path(od, "validation_data_overstory.csv"))
pb <- function(p, o) 100*(p/o - 1)

acd <- d[d$VARIANT=="ACD" & is.finite(d$BA_t2) & d$BA_t2>0, ]
ne  <- d[d$VARIANT=="NE"  & is.finite(d$BA_t2) & d$BA_t2>0, ]

acd_cal <- pb(acd$BA_pred_calib,   acd$BA_t2); acd_def <- pb(acd$BA_pred_default, acd$BA_t2)
ne_cal  <- pb(ne$BA_pred_calib,    ne$BA_t2);  ne_def  <- pb(ne$BA_pred_default,  ne$BA_t2)

# trim extreme tails for plotting (keep distribution visible)
clip <- function(x) pmax(pmin(x, 80), -60)

png(file.path(od, "fia_overstory_perplot_dist.png"), width = 2300, height = 950, res = 170)
par(mfrow = c(1, 2), mar = c(4.8, 4.6, 3.6, 1.0), mgp = c(2.7, 0.7, 0))

GREEN <- "#2ca02c"; GREY <- "#5A5A5A"; BLUE <- "#7DB5D5"

## Panel A: ACD per-plot bias density (calibrated vs default)
plot(density(clip(acd_cal), bw=2), col=GREEN, lwd=2.5, xlim=c(-60,80),
     xlab="Per-plot BA bias (%)", main=sprintf("FVS-ACD per-plot bias  (n=%d)", nrow(acd)),
     las=1, font.main=2, cex.main=1.1)
lines(density(clip(acd_def), bw=2), col=GREY, lwd=2.5)
abline(v=0, col="#888", lty=2)
abline(v=median(acd_cal,na.rm=TRUE), col=GREEN, lty=3)
abline(v=median(acd_def,na.rm=TRUE), col=GREY, lty=3)
acd_bom_cal <- 100*(mean(acd$BA_pred_calib)/mean(acd$BA_t2)-1)
acd_bom_def <- 100*(mean(acd$BA_pred_default)/mean(acd$BA_t2)-1)
legend("topright", legend=c(
  sprintf("calibrated  bias-of-means %+.2f%%  median plot %+.1f%%", acd_bom_cal, median(acd_cal,na.rm=TRUE)),
  sprintf("default     bias-of-means %+.2f%%  median plot %+.1f%%", acd_bom_def, median(acd_def,na.rm=TRUE))
), col=c(GREEN, GREY), lwd=2.5, bty="n", cex=0.82)
mtext("dashed line = 0 ; dotted = median per-plot ; per-plot spread is wide (the calibration is landscape-mean accurate, not per-stand exact)",
      side=1, line=3.4, cex=0.72, col="#444")

## Panel B: NE
plot(density(clip(ne_cal), bw=2), col=BLUE, lwd=2.5, xlim=c(-60,80),
     xlab="Per-plot BA bias (%)", main=sprintf("FVS-NE per-plot bias  (n=%d)", nrow(ne)),
     las=1, font.main=2, cex.main=1.1)
lines(density(clip(ne_def), bw=2), col=GREY, lwd=2.5)
abline(v=0, col="#888", lty=2)
abline(v=median(ne_cal,na.rm=TRUE), col=BLUE, lty=3)
abline(v=median(ne_def,na.rm=TRUE), col=GREY, lty=3)
ne_bom_cal <- 100*(mean(ne$BA_pred_calib)/mean(ne$BA_t2)-1)
ne_bom_def <- 100*(mean(ne$BA_pred_default)/mean(ne$BA_t2)-1)
legend("topright", legend=c(
  sprintf("calibrated  bias-of-means %+.2f%%  median plot %+.1f%%", ne_bom_cal, median(ne_cal,na.rm=TRUE)),
  sprintf("default     bias-of-means %+.2f%%  median plot %+.1f%%", ne_bom_def, median(ne_def,na.rm=TRUE))
), col=c(BLUE, GREY), lwd=2.5, bty="n", cex=0.82)
mtext("bias-of-means = sum(pred)/sum(obs)-1 ; calibration narrows both stats toward zero",
      side=1, line=3.4, cex=0.72, col="#444")

dev.off()
cat("wrote", file.path(od, "fia_overstory_perplot_dist.png"), "\n")
