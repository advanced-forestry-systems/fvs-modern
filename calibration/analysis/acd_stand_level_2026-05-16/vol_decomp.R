# Volume bias decomposition: is the -4.7% VOL bias explainable by DIA/HT bias,
# or is it residual to the volume equation itself?
#
# VOL ~ TPA * (a * DIA^b * HT^c)   (Honer-style; NVEL uses regional variants)
# If pred VOL/obs VOL ratio = (pred TPA/obs TPA) * (pred DIA/obs DIA)^b * (pred HT/obs HT)^c
# at population level, then the residual is the volume-equation issue.

vd <- read.csv("validation_data_acd_post.csv", stringsAsFactors = FALSE)

# Use only rows with non-NA volumes
v <- !is.na(vd$VOL_CFGRS_pred_calib) & !is.na(vd$VOL_CFGRS_t2) &
     vd$VOL_CFGRS_t2 > 0 & vd$QMD_t2 > 0 & vd$HT_top_t2 > 0 &
     !is.na(vd$QMD_pred_calib) & !is.na(vd$HT_top_calib)
d <- vd[v, ]

# Population-level ratios
r_TPA <- mean(d$TPA_pred_calib) / mean(d$TPA_t2)
r_QMD <- mean(d$QMD_pred_calib) / mean(d$QMD_t2)
r_HT  <- mean(d$HT_top_calib)   / mean(d$HT_top_t2)
r_VOL <- mean(d$VOL_CFGRS_pred_calib) / mean(d$VOL_CFGRS_t2)

cat("Population pred/obs ratios (calibrated arm, n =", nrow(d), ")\n")
cat("--------------------------------------------------------\n")
cat(sprintf("  TPA   ratio = %.4f   bias = %+.2f%%\n", r_TPA, 100*(r_TPA-1)))
cat(sprintf("  QMD   ratio = %.4f   bias = %+.2f%%\n", r_QMD, 100*(r_QMD-1)))
cat(sprintf("  HT    ratio = %.4f   bias = %+.2f%%\n", r_HT,  100*(r_HT-1)))
cat(sprintf("  VOL   ratio = %.4f   bias = %+.2f%%\n\n", r_VOL, 100*(r_VOL-1)))

# What would VOL ratio be if it equaled the Honer prediction
# VOL ~ TPA * DIA^2 * HT  (b=2, c=1 — common cubic-foot form)
predicted_VOL_ratio_b2c1 <- r_TPA * r_QMD^2 * r_HT
cat("Honer/NVEL approximations of VOL ratio implied by marginal ratios:\n")
cat(sprintf("  TPA * QMD^2 * HT          (cubic form, b=2,c=1)   = %.4f  (%+.2f%%)\n",
            predicted_VOL_ratio_b2c1, 100*(predicted_VOL_ratio_b2c1-1)))

predicted_VOL_ratio_b18c11 <- r_TPA * r_QMD^1.8 * r_HT^1.1
cat(sprintf("  TPA * QMD^1.8 * HT^1.1    (Schumacher-Hall-like)   = %.4f  (%+.2f%%)\n",
            predicted_VOL_ratio_b18c11, 100*(predicted_VOL_ratio_b18c11-1)))

cat("\n")
cat("If actual VOL bias != approximation, the gap is the volume-equation residual.\n")
cat(sprintf("  observed VOL bias                   : %+.2f%%\n", 100*(r_VOL-1)))
cat(sprintf("  approximated bias (cubic form)      : %+.2f%%\n", 100*(predicted_VOL_ratio_b2c1-1)))
cat(sprintf("  unexplained residual (vol equation) : %+.2f%%\n",
            100*(r_VOL/predicted_VOL_ratio_b2c1 - 1)))
cat("\n")

# Per-record decomposition: where DIA and HT predictions are individually accurate,
# does the volume prediction agree with what DIA*HT would suggest?
# Compute per-tree-ish proxy: BA * HT_top (a rough basis for volume)
# and compare to VOL.
d$prod_obs   <- d$BA_t2          * d$HT_top_t2
d$prod_pred  <- d$BA_pred_calib  * d$HT_top_calib
d$resid_prod <- d$prod_pred - d$prod_obs
d$resid_VOL  <- d$VOL_CFGRS_pred_calib - d$VOL_CFGRS_t2

cat("Correlation of VOL residuals with BA*HT residuals\n")
cat(sprintf("  cor(VOL_resid, (BA*HT)_resid) = %.3f\n",
            cor(d$resid_VOL, d$resid_prod, use = "pairwise.complete")))

# Fit a quick OLS: VOL_resid = a + b*(BA*HT)_resid + e
m <- lm(resid_VOL ~ resid_prod, data = d)
cat(sprintf("  Slope from VOL_resid ~ (BA*HT)_resid: %.3f\n",
            coef(m)["resid_prod"]))
cat(sprintf("  Intercept (mean residual not explained by BA*HT): %.2f\n",
            coef(m)["(Intercept)"]))
cat(sprintf("  R^2 of decomposition fit: %.3f\n\n", summary(m)$r.squared))

# Per-record VOL bias when QMD and HT predictions are individually good (within +-5%)
ok_qmd <- abs(d$QMD_pred_calib / d$QMD_t2 - 1) < 0.05
ok_ht  <- abs(d$HT_top_calib   / d$HT_top_t2 - 1) < 0.05
both_ok <- ok_qmd & ok_ht
cat(sprintf("On the %d records where pred QMD AND HT are within +-5%% of obs:\n",
            sum(both_ok)))
g <- d[both_ok, ]
cat(sprintf("  VOL pred/obs ratio = %.4f   bias = %+.2f%%\n",
            mean(g$VOL_CFGRS_pred_calib) / mean(g$VOL_CFGRS_t2),
            100*(mean(g$VOL_CFGRS_pred_calib)/mean(g$VOL_CFGRS_t2) - 1)))
cat("If this bias is still substantial when DIA + HT predictions are perfect,\n")
cat("the volume equation itself (NVEL, not the Bayesian posteriors) is responsible.\n")
