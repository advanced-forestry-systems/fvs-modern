#!/usr/bin/env Rscript
# silc_cfi_long_scorecard_v2.R
# =====================================================================
# Extends the long-horizon scorecard with the MORTCAL=TRUE AGM run.
# Produces:
#   silc_cfi_long_BA_v2.csv (with AGM MORTCAL row)
#   silc_cfi_long_Cords_v2.csv
#   silc_cfi_long_BdFt_v2.csv
#   silc_cfi_long_mortcal_compare.png (AGM default vs MORTCAL scatter)
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
SCRIB_TO_INTL <- 1.18

ps   <- read.csv(file.path(od, "silc_cfi_longhorizon_pairs.csv"))
ag   <- read.csv(file.path(od, "silc_cfi_long_agy_results.csv"))
agmc <- read.csv(file.path(od, "silc_cfi_long_agy_mortcal_results.csv"))
fvs  <- read.csv(file.path(od, "silc_cfi_long_fvs_results.csv"))
osm  <- read.csv(file.path(od, "silc_cfi_long_osm_results.csv"))
tr   <- read.csv(file.path(od, "TREE.csv"))
EXPF_CFI <- 5.0

obs_curr <- function(p, y) {
  t <- tr[tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN; h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  if (any(miss)) h[miss] <- pmax(6, 4.27 + 82*(1 - exp(-0.04 * (d[miss] * 2.54))))
  list(BA   = sum(0.005454 * d^2 * EXPF_CFI),
       Cords= sum(0.0025 * d^2 * h * 0.90 * EXPF_CFI) / 79,
       BdFt = sum(0.01 * d[d>=9]^2 * h[d>=9] * EXPF_CFI))
}
ps$BA_o    <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$BA)
ps$Cords_o <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$Cords)
ps$BdFt_o  <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$BdFt)

# Merge predictions
ag2  <- ag[, c("PLOT","BA_PRED_ft2ac","Cords_PRED_ac","BdFt_intl_PRED_ac")]
names(ag2) <- c("PLOT","BA_agm","Cords_agm","BdFt_agm")
ag3 <- agmc[, c("PLOT","BA_PRED_ft2ac","Cords_PRED_ac","BdFt_intl_PRED_ac")]
names(ag3) <- c("PLOT","BA_agmmc","Cords_agmmc","BdFt_agmmc")
fvs$model <- paste(fvs$variant, fvs$config, sep="_")
fvs$BdFt_intl <- fvs$BdFt_PRED * SCRIB_TO_INTL
fvs$Cords <- fvs$MCuFt_PRED / 79
fvs_w <- reshape(fvs[, c("PLOT","YEAR_PREV","model","BA_PRED_ft2ac","Cords","BdFt_intl")],
                 idvar=c("PLOT","YEAR_PREV"), timevar="model", direction="wide")
osm2 <- osm[, c("PLOT","BA_PRED_ft2ac")]; names(osm2) <- c("PLOT","BA_osm")

m <- merge(ps, ag2,  by="PLOT", all.x=TRUE)
m <- merge(m,  ag3,  by="PLOT", all.x=TRUE)
m <- merge(m,  fvs_w, by=c("PLOT","YEAR_PREV"), all.x=TRUE)
m <- merge(m,  osm2, by="PLOT", all.x=TRUE)

m$estbl <- (m$BA_CURR_FT2AC / pmax(m$BA_PREV_FT2AC, 1)) > 2.0
core <- m[!m$estbl, ]

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse <- function(p, o) sqrt(mean((p-o)^2, na.rm=TRUE))
r2 <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); if (sum(ok)<3) return(NA_real_)
  1 - sum((o[ok]-p[ok])^2) / sum((o[ok]-mean(o[ok]))^2)
}
mk <- function(lab, p, o) data.frame(predictor=lab, n=sum(is.finite(p)&is.finite(o)),
                                      bias_pct=round(bias_pct(p,o),2),
                                      RMSE=round(rmse(p,o),2), R2=round(r2(p,o),3))

cat("=== ROUTINE n=", nrow(core), " (mean ", round(mean(core$PERIOD_YR),1), " yr) ===\n", sep="")
cat("\n--- BA (ft^2/ac) ---\n")
ba <- rbind(
  mk("zero growth",       core$BA_PREV_FT2AC,            core$BA_CURR_FT2AC),
  mk("AGM default",       core$BA_agm,                   core$BA_CURR_FT2AC),
  mk("AGM MORTCAL",       core$BA_agmmc,                 core$BA_CURR_FT2AC),
  mk("FVS-NE calibrated", core$BA_PRED_ft2ac.NE_calibrated, core$BA_CURR_FT2AC),
  mk("FVS-ACD default",   core$BA_PRED_ft2ac.ACD_default,   core$BA_CURR_FT2AC),
  mk("OSM-ACD",           core$BA_osm,                   core$BA_CURR_FT2AC)
)
print(ba, row.names=FALSE)

cat("\n--- Cords/ac ---\n")
cs <- rbind(
  mk("zero growth",       core$Cords_o*0+core$Cords_o,    core$Cords_o), # zero growth proxy
  mk("AGM default",       core$Cords_agm,                 core$Cords_o),
  mk("AGM MORTCAL",       core$Cords_agmmc,               core$Cords_o),
  mk("FVS-NE calibrated", core$Cords.NE_calibrated,       core$Cords_o),
  mk("FVS-ACD default",   core$Cords.ACD_default,         core$Cords_o)
)
print(cs, row.names=FALSE)

cat("\n--- BdFt Intl 1/4 /ac ---\n")
bd <- rbind(
  mk("AGM default",       core$BdFt_agm,                  core$BdFt_o),
  mk("AGM MORTCAL",       core$BdFt_agmmc,                core$BdFt_o),
  mk("FVS-NE calibrated", core$BdFt_intl.NE_calibrated,   core$BdFt_o),
  mk("FVS-ACD default",   core$BdFt_intl.ACD_default,     core$BdFt_o)
)
print(bd, row.names=FALSE)

write.csv(ba, file.path(od, "silc_cfi_long_BA_v2.csv"), row.names=FALSE)
write.csv(cs, file.path(od, "silc_cfi_long_Cords_v2.csv"), row.names=FALSE)
write.csv(bd, file.path(od, "silc_cfi_long_BdFt_v2.csv"), row.names=FALSE)

# AGM default vs MORTCAL comparison figure
png(file.path(od, "silc_cfi_long_mortcal_compare.png"),
    width=2400, height=900, res=170)
par(mfrow=c(1,3), mar=c(4.5, 4.6, 3.4, 1.0), mgp=c(2.7, 0.6, 0))
CRSF_GREEN <- "#1A3D28"; MORT_GOLD <- "#B8860B"
draw <- function(o, def_p, mc_p, lab) {
  lim <- c(0, max(c(o, def_p, mc_p), na.rm=TRUE) * 1.05)
  plot(NA, xlim=lim, ylim=lim, xlab="Observed", ylab="Predicted",
       main=lab, las=1, font.main=2, cex.main=1.2)
  abline(0, 1, lty=2, col="#888")
  points(o, def_p, col=CRSF_GREEN, pch=19, cex=1.5)
  points(o, mc_p,  col=MORT_GOLD,  pch=17, cex=1.5)
  legend("topleft", legend=c(
    sprintf("AGM default (%.1f%%, R^2 %.2f)", bias_pct(def_p, o), r2(def_p, o)),
    sprintf("AGM MORTCAL (%.1f%%, R^2 %.2f)", bias_pct(mc_p,  o), r2(mc_p,  o))
  ), col=c(CRSF_GREEN, MORT_GOLD), pch=c(19,17), bty="n", cex=0.95)
}
draw(core$BA_CURR_FT2AC, core$BA_agm,   core$BA_agmmc,   "BA (ft^2/ac)")
draw(core$Cords_o,       core$Cords_agm, core$Cords_agmmc, "Cords/ac")
draw(core$BdFt_o,        core$BdFt_agm,  core$BdFt_agmmc,  "BdFt Intl 1/4 /ac")
mtext(sprintf("AGM default vs MORTCAL=TRUE (#126b in-source size-dep mortality, n=%d routine)", nrow(core)),
      side=3, line=-1.2, outer=TRUE, font=2, cex=1.05)
dev.off()
cat("\nwrote silc_cfi_long_mortcal_compare.png\n")
