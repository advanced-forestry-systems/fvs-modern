vd <- read.csv("validation_data_acd_post.csv", stringsAsFactors = FALSE)
cat("ACD calibrated vs observed — all four core attrs (n =", nrow(vd), ")\n")
cat(strrep("=", 80), "\n\n")

f <- function(pred, obs, label) {
  v <- !is.na(pred) & !is.na(obs) & obs > 0
  p <- pred[v]; o <- obs[v]
  resid <- p - o
  cat(sprintf("%-7s  n=%5d  mean_obs=%7.2f  mean_pred=%7.2f  pred/obs=%.3f\n",
              label, length(p), mean(o), mean(p), mean(p)/mean(o)))
  cat(sprintf("         raw bias = %+7.2f   pop bias%% = %+6.2f%%\n",
              mean(resid), 100*mean(resid)/mean(o)))
  cat(sprintf("         RMSE     = %7.2f   RMSE%%    = %6.2f%%   R2 = %.3f\n",
              sqrt(mean(resid^2)),
              100*sqrt(mean(resid^2))/mean(o),
              max(1 - sum(resid^2)/sum((o-mean(o))^2), 0)))
  cat(sprintf("         median raw resid = %+.2f   median %%resid = %+.1f%%\n\n",
              median(resid), 100 * median(resid/o)))
}

f(vd$BA_pred_calib,         vd$BA_t2,         "BA")
f(vd$TPA_pred_calib,        vd$TPA_t2,        "TPA")
f(vd$QMD_pred_calib,        vd$QMD_t2,        "QMD")
f(vd$VOL_CFGRS_pred_calib,  vd$VOL_CFGRS_t2,  "VOL")
f(vd$HT_top_calib,          vd$HT_top_t2,     "HT_top")

cat("Same comparison vs DEFAULT arm (to gauge what calibration is buying us)\n")
cat(strrep("-", 80), "\n\n")
f(vd$BA_pred_default,        vd$BA_t2,         "BA-def")
f(vd$VOL_CFGRS_pred_default, vd$VOL_CFGRS_t2,  "VOL-def")
f(vd$HT_top_default,         vd$HT_top_t2,     "HT-def")

cat("Where does the BA error concentrate? (component decomposition)\n")
cat(strrep("-", 80), "\n\n")
# BA = sum(TPA * pi*(DIA/24)^2). If we missed in TPA, that pulls BA. If we missed
# in QMD (which derives from BA/TPA), that's a different signal.
cat("  Are TPA and QMD biases consistent with the BA bias?\n")
cat(sprintf("    TPA  pred/obs ratio = %.3f\n", mean(vd$TPA_pred_calib)/mean(vd$TPA_t2)))
cat(sprintf("    QMD  pred/obs ratio = %.3f\n", mean(vd$QMD_pred_calib)/mean(vd$QMD_t2)))
cat(sprintf("    BA   pred/obs ratio = %.3f\n", mean(vd$BA_pred_calib)/mean(vd$BA_t2)))
cat("\n")

# Restrict to mature stands (BA > 50) to see how it performs where it should
mature <- vd[vd$BA_t2 > 50, ]
cat(sprintf("On mature stands (BA > 50, n=%d, %.1f%% of sample):\n",
            nrow(mature), 100*nrow(mature)/nrow(vd)))
f(mature$BA_pred_calib, mature$BA_t2, "BA")
f(mature$VOL_CFGRS_pred_calib, mature$VOL_CFGRS_t2, "VOL")
