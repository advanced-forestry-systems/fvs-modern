#!/usr/bin/env Rscript
# silc_cfi_scorecard.R
# =====================================================================
# Build the SILC CFI predicted-vs-observed scorecard combining:
#   * AcadianGY (from Cardinal run silc_cfi_acadiangy_pred.csv)
#   * naive baselines (zero-growth and FIA-prior, from pair_summary)
# Exclude establishment pairs (BA_prev < 10 ft^2/ac) where the
# observed change is dominated by ingrowth crossing the 4.5 in
# measurement threshold, not regular stand growth.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
ps <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"),
               stringsAsFactors = FALSE)
ag <- read.csv(file.path(od, "silc_cfi_acadiangy_pred.csv"),
               stringsAsFactors = FALSE)

m <- merge(ps, ag, by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"),
           suffixes = c("","_agy"))

# Tag establishment pairs (BA_prev < 10 ft^2/ac OR PAI_NET_OBS > 5 ft^2/ac/yr)
m$establishment <- m$BA_PREV_FT2AC < 10 | abs(m$PAI_NET_OBS) > 5
cat(sprintf("Establishment pairs (excluded from routine scorecard): %d / %d\n",
            sum(m$establishment), nrow(m)))
core <- m[!m$establishment, ]

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse     <- function(p, o) sqrt(mean((p - o)^2, na.rm=TRUE))
r2       <- function(p, o) {
  ok <- is.finite(p) & is.finite(o)
  if (sum(ok) < 3) return(NA_real_)
  1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2)
}

# Scorecard table
make_row <- function(label, p, o, df) {
  data.frame(
    predictor      = label,
    n              = sum(is.finite(p) & is.finite(o)),
    bias_pct       = round(bias_pct(p, o), 2),
    RMSE_ft2_ac    = round(rmse(p, o), 2),
    R2             = round(r2(p, o), 3)
  )
}

cat("\n=== Full sample (n=", nrow(m), ") ===\n", sep="")
full_sc <- rbind(
  make_row("zero growth",     m$BA_pred_zero_growth,  m$BA_CURR_FT2AC, m),
  make_row("FIA prior PAI",   m$BA_pred_FIA_prior,    m$BA_CURR_FT2AC, m),
  make_row("AcadianGY 12.3.9",m$BA_PRED_ft2ac,        m$BA_CURR_FT2AC, m)
)
print(full_sc, row.names=FALSE)

cat(sprintf("\n=== Routine-growth subset (n=%d, excluding %d establishment pairs) ===\n",
            nrow(core), sum(m$establishment)))
core_sc <- rbind(
  make_row("zero growth",     core$BA_pred_zero_growth,  core$BA_CURR_FT2AC),
  make_row("FIA prior PAI",   core$BA_pred_FIA_prior,    core$BA_CURR_FT2AC),
  make_row("AcadianGY 12.3.9",core$BA_PRED_ft2ac,        core$BA_CURR_FT2AC)
)
print(core_sc, row.names=FALSE)

write.csv(full_sc,
          file.path(od, "silc_cfi_scorecard_full.csv"),
          row.names = FALSE)
write.csv(core_sc,
          file.path(od, "silc_cfi_scorecard_routine.csv"),
          row.names = FALSE)
write.csv(m,
          file.path(od, "silc_cfi_pair_with_pred.csv"),
          row.names = FALSE)

# Figure: predicted vs observed BA scatter (routine sample), three models
png(file.path(od, "silc_cfi_scorecard_scatter.png"),
    width = 2300, height = 950, res = 165)
par(mfrow = c(1, 3), mar = c(4.5, 4.6, 3.4, 1.0),
    mgp = c(2.7, 0.6, 0))

draw_scatter <- function(pred, obs, label, col_main) {
  lim <- c(0, max(c(pred, obs), na.rm=TRUE) * 1.05)
  plot(obs, pred, pch = 19, col = col_main, cex = 1.4,
       xlim = lim, ylim = lim,
       xlab = "Observed BA at year_curr (ft^2/ac)",
       ylab = "Predicted BA (ft^2/ac)",
       main = label, las = 1, font.main = 2)
  abline(0, 1, lty = 2, col = "#888")
  abline(lm(pred ~ obs), col = col_main, lwd = 1.5)
  bp <- bias_pct(pred, obs); rm <- rmse(pred, obs); r2v <- r2(pred, obs)
  legend("topleft", legend = c(
    sprintf("bias %+.2f%%", bp),
    sprintf("RMSE %.1f ft^2/ac", rm),
    sprintf("R^2 %.2f", r2v),
    sprintf("n = %d", sum(is.finite(pred)&is.finite(obs)))
    ), bty = "n", cex = 0.95)
}
draw_scatter(core$BA_pred_zero_growth, core$BA_CURR_FT2AC,
             "zero growth", "#888888")
draw_scatter(core$BA_pred_FIA_prior,   core$BA_CURR_FT2AC,
             "FIA prior PAI", "#7DB5D5")
draw_scatter(core$BA_PRED_ft2ac,       core$BA_CURR_FT2AC,
             "AcadianGY 12.3.9", "#2ca02c")
mtext("SILC CFI: predicted vs observed BA at year_curr  (routine-growth subset; establishment pairs excluded)",
      side = 3, line = -1.4, outer = TRUE, cex = 1.0, font = 2)
dev.off()

cat("\nwrote silc_cfi_scorecard_full.csv\n")
cat("wrote silc_cfi_scorecard_routine.csv\n")
cat("wrote silc_cfi_pair_with_pred.csv\n")
cat("wrote silc_cfi_scorecard_scatter.png\n")
