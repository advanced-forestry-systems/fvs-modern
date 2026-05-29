#!/usr/bin/env Rscript
# silc_cfi_merch_scorecard.R
# =====================================================================
# Compute CFI observed cords + BdFt at year_curr (apples-to-apples
# with the AcadianGY v3 predictions), then build the merch-volume
# scorecard.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
EXPF <- 5.0   # CFI 1/5-acre plots: each tree -> 5 trees/ac

tr <- read.csv(file.path(od, "TREE.csv"))
ag <- read.csv(file.path(od, "silc_cfi_acadiangy_pred_v3.csv"))
ps <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"))

# Observed metrics at year_curr per pair
obs_curr <- function(p, y) {
  t <- tr[tr$PLOT == p & tr$MEASYEAR == y &
          tr$STATUSCD == 1 & is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN
  h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  # imputation matches the driver: HT_m = pmax(2, 1.3 + 25*(1-exp(-0.04*DBH_cm)))
  if (any(miss))
    h[miss] <- pmax(6, 4.27 + 82 * (1 - exp(-0.04 * (d[miss] * 2.54))))
  e <- rep(EXPF, length(d))
  tcuft <- 0.0025 * d^2 * h
  merch <- tcuft * 0.90
  saw   <- d >= 9.0
  saw_cuft <- ifelse(saw, tcuft * 0.55, 0)
  list(
    BA      = sum(0.005454 * d^2 * e),
    Cords   = sum(merch * e) / 79,
    BdFt    = sum(saw_cuft * e) * 6.0,
    n_live  = nrow(t)
  )
}

ag$Cords_OBS_CURR <- NA_real_
ag$BdFt_OBS_CURR  <- NA_real_
for (i in seq_len(nrow(ag))) {
  o <- obs_curr(ag$PLOT[i], ag$YEAR_CURR[i])
  ag$Cords_OBS_CURR[i] <- o$Cords
  ag$BdFt_OBS_CURR[i]  <- o$BdFt
}

# Baselines: zero growth (Cords/BdFt at year_prev) and FIA-prior PAI
ag$Cords_zero <- ag$Cords_OBS_PREV
ag$BdFt_zero  <- ag$BdFt_OBS_PREV

# FIA prior PAI for cords: typical Acadian PAI ~ 0.20 cords/ac/yr
prior_cords_PAI <- c("Cedar"="0.10","Hardwood"="0.18",
                      "Mixedwood"="0.22","Commercial Softwood"="0.25",
                      "Other Softwood"="0.18","Unclassifiable"=NA)
strata_map <- read.csv(file.path(od, "silc_cfi_plot_strata_map.csv"))
ft <- setNames(strata_map$forest_type, strata_map$PLOT)
ag$forest_type <- ft[as.character(ag$PLOT)]
ag$prior_cords_PAI <- as.numeric(prior_cords_PAI[ag$forest_type])
ag$Cords_FIAprior  <- ag$Cords_OBS_PREV +
                       ag$prior_cords_PAI * ag$PERIOD_YR

# Establishment exclusion (consistent with BA scorecard)
ag$establishment <- ag$BA_OBS_PREV < 10 | abs(ag$PAI_NET_OBS) > 5

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse     <- function(p, o) sqrt(mean((p - o)^2, na.rm=TRUE))
r2       <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); if (sum(ok)<3) return(NA_real_)
  1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2)
}

make_row <- function(label, p, o, unit) {
  data.frame(predictor = label, n = sum(is.finite(p)&is.finite(o)),
             bias_pct = round(bias_pct(p, o), 2),
             RMSE = round(rmse(p, o), 2),
             RMSE_unit = unit,
             R2 = round(r2(p, o), 3))
}

cat("=== CFI scorecard at year_curr (full sample n=", nrow(ag), ") ===\n",
    sep="")
core <- ag[!ag$establishment, ]

cat("\n--- BA (ft^2/ac) ---\n")
print(rbind(
  make_row("zero growth",     ag$BA_OBS_PREV, ag$BA_OBS_CURR, "ft^2/ac"),
  make_row("AcadianGY 12.3.9",ag$BA_PRED_ft2ac, ag$BA_OBS_CURR, "ft^2/ac")
), row.names=FALSE)

cat("\n--- Cords/ac ---\n")
print(rbind(
  make_row("zero growth",     ag$Cords_zero,     ag$Cords_OBS_CURR, "cords/ac"),
  make_row("FIA prior PAI",   ag$Cords_FIAprior, ag$Cords_OBS_CURR, "cords/ac"),
  make_row("AcadianGY 12.3.9",ag$Cords_PRED_ac,  ag$Cords_OBS_CURR, "cords/ac")
), row.names=FALSE)

cat("\n--- BdFt/ac ---\n")
print(rbind(
  make_row("zero growth",     ag$BdFt_zero,    ag$BdFt_OBS_CURR, "bd ft/ac"),
  make_row("AcadianGY 12.3.9",ag$BdFt_PRED_ac, ag$BdFt_OBS_CURR, "bd ft/ac")
), row.names=FALSE)

cat(sprintf("\n=== Routine-growth subset (n=%d, excluding %d establishment pairs) ===\n",
            nrow(core), sum(ag$establishment)))

cat("\n--- BA (ft^2/ac) ---\n")
print(rbind(
  make_row("zero growth",     core$BA_OBS_PREV, core$BA_OBS_CURR, "ft^2/ac"),
  make_row("AcadianGY 12.3.9",core$BA_PRED_ft2ac, core$BA_OBS_CURR, "ft^2/ac")
), row.names=FALSE)

cat("\n--- Cords/ac ---\n")
core_cords <- rbind(
  make_row("zero growth",     core$Cords_zero,     core$Cords_OBS_CURR, "cords/ac"),
  make_row("FIA prior PAI",   core$Cords_FIAprior, core$Cords_OBS_CURR, "cords/ac"),
  make_row("AcadianGY 12.3.9",core$Cords_PRED_ac,  core$Cords_OBS_CURR, "cords/ac")
)
print(core_cords, row.names=FALSE)

cat("\n--- BdFt/ac ---\n")
core_bdft <- rbind(
  make_row("zero growth",     core$BdFt_zero,    core$BdFt_OBS_CURR, "bd ft/ac"),
  make_row("AcadianGY 12.3.9",core$BdFt_PRED_ac, core$BdFt_OBS_CURR, "bd ft/ac")
)
print(core_bdft, row.names=FALSE)

write.csv(ag, file.path(od, "silc_cfi_merch_predictions.csv"),
          row.names = FALSE)
write.csv(core_cords, file.path(od, "silc_cfi_merch_scorecard_cords.csv"),
          row.names = FALSE)
write.csv(core_bdft, file.path(od, "silc_cfi_merch_scorecard_bdft.csv"),
          row.names = FALSE)

# Figure: predicted vs observed cords scatter
png(file.path(od, "silc_cfi_merch_scatter.png"),
    width = 2300, height = 950, res = 165)
par(mfrow = c(1, 3), mar = c(4.5, 4.6, 3.4, 1.0), mgp = c(2.7, 0.6, 0))

draw <- function(pred, obs, label, col_main, ulab) {
  ok <- is.finite(pred) & is.finite(obs)
  lim <- c(0, max(c(pred[ok], obs[ok])) * 1.05)
  plot(obs, pred, pch = 19, col = col_main, cex = 1.4,
       xlim = lim, ylim = lim,
       xlab = sprintf("Observed %s", ulab),
       ylab = sprintf("Predicted %s", ulab),
       main = label, las = 1, font.main = 2)
  abline(0, 1, lty = 2, col = "#888")
  if (sum(ok) >= 3) abline(lm(pred[ok] ~ obs[ok]), col = col_main, lwd = 1.5)
  bp <- bias_pct(pred, obs); rm <- rmse(pred, obs); r2v <- r2(pred, obs)
  legend("topleft", legend = c(
    sprintf("bias %+.1f%%", bp),
    sprintf("RMSE %.2f", rm),
    sprintf("R^2 %.2f", r2v),
    sprintf("n = %d", sum(ok))), bty = "n", cex = 0.95)
}
draw(core$Cords_zero,     core$Cords_OBS_CURR,
     "zero growth (cords)", "#888888", "cords/ac")
draw(core$Cords_FIAprior, core$Cords_OBS_CURR,
     "FIA prior PAI (cords)", "#7DB5D5", "cords/ac")
draw(core$Cords_PRED_ac,  core$Cords_OBS_CURR,
     "AcadianGY 12.3.9 (cords)", "#2ca02c", "cords/ac")
mtext("SILC CFI: predicted vs observed merch cords at year_curr  (routine-growth subset)",
      side = 3, line = -1.5, outer = TRUE, cex = 1.0, font = 2)
dev.off()

cat("\nwrote silc_cfi_merch_predictions.csv, silc_cfi_merch_scorecard_*.csv\n")
cat("wrote silc_cfi_merch_scatter.png\n")
