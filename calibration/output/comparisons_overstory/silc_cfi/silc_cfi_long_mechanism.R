#!/usr/bin/env Rscript
# silc_cfi_long_mechanism.R
# =====================================================================
# Mechanism figure: decomposes AGM's +18.5% long-horizon BA bias into
# (1) under-mortality and (2) over-growth on diameter. Per-plot bars
# show observed vs AGM TPA loss and QMD growth side by side.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"

ps <- read.csv(file.path(od, "silc_cfi_longhorizon_pairs.csv"))
ag <- read.csv(file.path(od, "silc_cfi_long_agy_results.csv"))
m  <- merge(ps, ag[, c("PLOT","BA_PRED_ft2ac","TPA_PRED","QMD_PRED_in")], by = "PLOT")
m$obs_TPA_chg_pct <- 100 * (m$TPA_CURR / m$TPA_PREV - 1)
m$agm_TPA_chg_pct <- 100 * (m$TPA_PRED / m$TPA_PREV - 1)
m$obs_QMD_chg     <- m$QMD_CURR - m$QMD_PREV
m$agm_QMD_chg     <- m$QMD_PRED_in - m$QMD_PREV
m$BA_bias_pct     <- 100 * (m$BA_PRED_ft2ac / m$BA_CURR_FT2AC - 1)
m$estbl           <- (m$BA_CURR_FT2AC / m$BA_PREV_FT2AC) > 2.0

core <- m[!m$estbl, ]
core <- core[order(core$BA_bias_pct), ]
core$lab <- sprintf("Plot %d", core$PLOT)

CRSF_GREEN <- "#1A3D28"; OBS_GREY <- "#666666"
png(file.path(od, "silc_cfi_long_mechanism.png"),
    width = 2400, height = 1100, res = 170)
par(mfrow = c(1, 3), mar = c(4.5, 4.8, 3.2, 1.0), mgp = c(2.7, 0.6, 0))

# --- Panel 1: TPA change ---
ymat <- t(as.matrix(core[, c("obs_TPA_chg_pct", "agm_TPA_chg_pct")]))
barplot(ymat, beside = TRUE, names.arg = core$lab,
        col = c(OBS_GREY, CRSF_GREEN), las = 1,
        main = "TPA change over 10 years (%)",
        ylim = c(min(ymat, 0) * 1.1, max(ymat, 10) * 1.2),
        cex.main = 1.25, font.main = 2, ylab = "%")
abline(h = 0, lty = 1, col = "#333")
legend("topleft", legend = c("Observed", "AGM predicted"),
       fill = c(OBS_GREY, CRSF_GREEN), bty = "n", cex = 1.05)

# --- Panel 2: QMD growth ---
qmat <- t(as.matrix(core[, c("obs_QMD_chg", "agm_QMD_chg")]))
barplot(qmat, beside = TRUE, names.arg = core$lab,
        col = c(OBS_GREY, CRSF_GREEN), las = 1,
        main = "QMD growth over 10 years (in)",
        ylim = c(0, max(qmat) * 1.2),
        cex.main = 1.25, font.main = 2, ylab = "inches")
legend("topleft", legend = c("Observed", "AGM predicted"),
       fill = c(OBS_GREY, CRSF_GREEN), bty = "n", cex = 1.05)

# --- Panel 3: per-plot BA bias ---
bp <- barplot(core$BA_bias_pct, names.arg = core$lab, las = 1,
              col = ifelse(core$BA_bias_pct > 15, "#A03A2A", CRSF_GREEN),
              main = "AGM BA bias by plot (%)",
              ylab = "% over observed", cex.main = 1.25, font.main = 2,
              ylim = c(0, max(core$BA_bias_pct) * 1.2))
abline(h = mean(core$BA_bias_pct), lty = 2, col = "#A03A2A", lwd = 2)
text(bp, core$BA_bias_pct + 1.8,
     sprintf("%.0f%%", core$BA_bias_pct), cex = 0.95, font = 2)
text(par("usr")[2] * 0.96, mean(core$BA_bias_pct) + 2.5,
     sprintf("mean %.1f%%", mean(core$BA_bias_pct)),
     col = "#A03A2A", font = 2, cex = 0.95, adj = 1)

mtext("AGM long-horizon bias mechanism: under-mortality on disturbed plots + over-growth on slow plots",
      side = 3, line = -1.0, outer = TRUE, font = 2, cex = 1.1)
dev.off()
cat("wrote silc_cfi_long_mechanism.png\n")

# Companion CSV with per-plot decomposition
out <- core[, c("PLOT","PERIOD_YR","BA_PREV_FT2AC","BA_CURR_FT2AC",
                "BA_PRED_ft2ac","BA_bias_pct",
                "TPA_PREV","TPA_CURR","TPA_PRED",
                "obs_TPA_chg_pct","agm_TPA_chg_pct",
                "QMD_PREV","QMD_CURR","QMD_PRED_in",
                "obs_QMD_chg","agm_QMD_chg")]
write.csv(out, file.path(od, "silc_cfi_long_mechanism.csv"),
          row.names = FALSE)
cat("wrote silc_cfi_long_mechanism.csv\n")
