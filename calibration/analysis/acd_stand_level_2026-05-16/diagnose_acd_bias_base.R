#!/usr/bin/env Rscript
vd <- read.csv("validation_data_acd_post.csv", stringsAsFactors = FALSE)
cat("N =", nrow(vd), "\n\n")

resid <- vd$BA_pred_calib - vd$BA_t2
pct   <- 100 * resid / vd$BA_t2

cat("Headline restatements (calibrated arm)\n")
cat("--------------------------------------\n")
cat(sprintf("  mean raw residual (pred-obs) : %+.2f sq ft/ac\n", mean(resid, na.rm = TRUE)))
cat(sprintf("  median raw residual          : %+.2f sq ft/ac\n", median(resid, na.rm = TRUE)))
cat(sprintf("  mean of per-record pct err   : %+.2f%%\n", mean(pct, na.rm = TRUE)))
cat(sprintf("  median of per-record pct err : %+.2f%%\n", median(pct, na.rm = TRUE)))
cat(sprintf("  mean obs BA                  : %.1f\n", mean(vd$BA_t2, na.rm = TRUE)))
cat(sprintf("  mean pred BA                 : %.1f\n", mean(vd$BA_pred_calib, na.rm = TRUE)))
cat(sprintf("  pred / obs ratio             : %.3f\n\n",
            mean(vd$BA_pred_calib, na.rm = TRUE) / mean(vd$BA_t2, na.rm = TRUE)))

cat("Bias by BA strata\n")
cat("-----------------\n")
vd$BA_stratum <- cut(vd$BA_t2, breaks = c(0, 25, 50, 75, 100, 150, 200, 1000), right = FALSE)
out <- do.call(rbind, lapply(split(vd, vd$BA_stratum), function(d) {
  data.frame(stratum  = d$BA_stratum[1],
             n        = nrow(d),
             mean_obs = round(mean(d$BA_t2, na.rm = TRUE), 1),
             mean_pred = round(mean(d$BA_pred_calib, na.rm = TRUE), 1),
             bias_raw  = round(mean(d$BA_pred_calib - d$BA_t2, na.rm = TRUE), 2),
             bias_pct  = round(100 * mean((d$BA_pred_calib - d$BA_t2) / d$BA_t2, na.rm = TRUE), 1),
             rmse_pct  = round(100 * sqrt(mean((d$BA_pred_calib - d$BA_t2)^2, na.rm = TRUE)) /
                                mean(d$BA_t2, na.rm = TRUE), 1))
}))
print(out, row.names = FALSE)
cat("\n")

cat("Bias by interval_years\n")
cat("----------------------\n")
yr <- aggregate(cbind(
  bias_raw = vd$BA_pred_calib - vd$BA_t2,
  bias_rec_pct = 100 * (vd$BA_pred_calib - vd$BA_t2) / vd$BA_t2),
  by = list(yr = vd$interval_years), FUN = function(x) mean(x, na.rm = TRUE))
yr$n <- aggregate(vd$BA_t2, by = list(yr = vd$interval_years), length)$x
yr$bias_raw <- round(yr$bias_raw, 2)
yr$bias_rec_pct <- round(yr$bias_rec_pct, 1)
print(yr[order(yr$yr), ], row.names = FALSE)
cat("\n")

cat("BA<50 vs BA>=100 split\n")
cat("----------------------\n")
u50 <- vd[vd$BA_t2 < 50, ]
o100 <- vd[vd$BA_t2 >= 100, ]
cat(sprintf("  BA<50 (n=%d): bias_pct = %+.1f%%  raw bias = %+.1f\n",
            nrow(u50),
            100 * mean((u50$BA_pred_calib - u50$BA_t2) / u50$BA_t2, na.rm = TRUE),
            mean(u50$BA_pred_calib - u50$BA_t2, na.rm = TRUE)))
cat(sprintf("  BA>=100 (n=%d): bias_pct = %+.1f%%  raw bias = %+.1f\n",
            nrow(o100),
            100 * mean((o100$BA_pred_calib - o100$BA_t2) / o100$BA_t2, na.rm = TRUE),
            mean(o100$BA_pred_calib - o100$BA_t2, na.rm = TRUE)))
cat("\n")

cat("Calibrated vs Default arm comparison\n")
cat("------------------------------------\n")
cat(sprintf("  calibrated  raw bias = %+.2f   bias_pct = %+.2f%%\n",
            mean(vd$BA_pred_calib  - vd$BA_t2, na.rm = TRUE),
            100 * mean((vd$BA_pred_calib  - vd$BA_t2) / vd$BA_t2, na.rm = TRUE)))
cat(sprintf("  default     raw bias = %+.2f   bias_pct = %+.2f%%\n",
            mean(vd$BA_pred_default - vd$BA_t2, na.rm = TRUE),
            100 * mean((vd$BA_pred_default - vd$BA_t2) / vd$BA_t2, na.rm = TRUE)))
