#!/usr/bin/env Rscript
# silc_cfi_threeway_panel.R
# =====================================================================
# Deck headline figure: three metrics (BA / Cords / BdFt) side by side,
# showing predicted vs observed for AcadianGY 12.3.9 on the routine-
# growth CFI subset, with key stats annotated.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
ag <- read.csv(file.path(od, "silc_cfi_merch_predictions.csv"))
core <- ag[!ag$establishment, ]

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse     <- function(p, o) sqrt(mean((p - o)^2, na.rm=TRUE))
r2       <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); if (sum(ok)<3) return(NA_real_)
  1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2)
}

png(file.path(od, "silc_cfi_threeway_panel.png"),
    width = 2700, height = 1050, res = 175)
par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.4, 1.0), mgp = c(2.8, 0.7, 0))

draw_metric <- function(pred, obs, title_label, unit, col_main) {
  ok <- is.finite(pred) & is.finite(obs)
  lim <- c(0, max(c(pred[ok], obs[ok])) * 1.06)
  plot(obs, pred, pch = 19, col = col_main, cex = 1.6,
       xlim = lim, ylim = lim,
       xlab = sprintf("Observed (%s)", unit),
       ylab = sprintf("Predicted (%s)", unit),
       main = title_label, las = 1, font.main = 2, cex.main = 1.25)
  abline(0, 1, lty = 2, col = "#888", lwd = 1.2)
  if (sum(ok) >= 3) abline(lm(pred[ok] ~ obs[ok]),
                            col = col_main, lwd = 2.2)
  bp <- bias_pct(pred, obs); rm <- rmse(pred, obs); r2v <- r2(pred, obs)
  legend("topleft", legend = c(
    sprintf("Bias    %+.1f%%", bp),
    sprintf("RMSE  %.2f", rm),
    sprintf("R^2    %.2f", r2v),
    sprintf("n        %d", sum(ok))
    ), bty = "n", cex = 1.05)
}

CRSF_GREEN <- "#1A3D28"
draw_metric(core$BA_PRED_ft2ac, core$BA_OBS_CURR,
            "Basal area",      "ft^2/ac", CRSF_GREEN)
draw_metric(core$Cords_PRED_ac, core$Cords_OBS_CURR,
            "Merchantable cords","cords/ac", CRSF_GREEN)
draw_metric(core$BdFt_PRED_ac,  core$BdFt_OBS_CURR,
            "Sawlog volume",    "BdFt/ac", CRSF_GREEN)
dev.off()
cat("wrote silc_cfi_threeway_panel.png\n")
