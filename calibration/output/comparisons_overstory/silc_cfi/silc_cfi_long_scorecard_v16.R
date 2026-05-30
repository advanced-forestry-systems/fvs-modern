#!/usr/bin/env Rscript
# silc_cfi_long_scorecard_v16.R
# Build scorecard CSVs and the long-horizon scatter figure for the v16
# AcadianGY MORTCAL run (n=73 pairs, 65 routine).
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
ps <- read.csv(file.path(od, "silc_cfi_longhorizon_pairs_v16.csv"))
mc <- read.csv(file.path(od, "silc_cfi_long_agy_mortcal_v16_results.csv"))
tr <- read.csv(file.path(od, "v16/TREE.csv"))
EXPF <- 5.0

obs_curr <- function(p, y) {
  t <- tr[tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN; h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  if (any(miss)) h[miss] <- pmax(6, 4.27 + 82*(1 - exp(-0.04 * (d[miss] * 2.54))))
  list(BA=sum(0.005454*d^2*EXPF), Cords=sum(0.0025*d^2*h*0.90*EXPF)/79,
       BdFt=sum(0.01*d[d>=9]^2*h[d>=9]*EXPF))
}
ps$BA_o    <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$BA)
ps$Cords_o <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$Cords)
ps$BdFt_o  <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$BdFt)

m <- merge(mc, ps[, c("PLOT","YEAR_PREV","BA_o","Cords_o","BdFt_o")],
           by=c("PLOT","YEAR_PREV"), all.x=TRUE)
m$estbl <- (m$BA_OBS_CURR / pmax(m$BA_OBS_PREV, 1)) > 2.0
core <- m[!m$estbl, ]

bias_pct <- function(p,o) 100*(mean(p,na.rm=TRUE)/mean(o,na.rm=TRUE)-1)
rmse <- function(p,o) sqrt(mean((p-o)^2,na.rm=TRUE))
r2 <- function(p,o) { ok <- is.finite(p)&is.finite(o); if(sum(ok)<3) return(NA_real_)
  1 - sum((o[ok]-p[ok])^2)/sum((o[ok]-mean(o[ok]))^2) }
mk <- function(lab,p,o) data.frame(predictor=lab, n=sum(is.finite(p)&is.finite(o)),
  bias_pct=round(bias_pct(p,o),2), RMSE=round(rmse(p,o),2), R2=round(r2(p,o),3))

cat(sprintf("=== v16 ROUTINE n=%d (mean horizon %.1f yr) ===\n", nrow(core), mean(core$PERIOD_YR)))
ba <- rbind(mk("zero growth", core$BA_OBS_PREV, core$BA_OBS_CURR),
            mk("AcadianGY MORTCAL", core$BA_PRED_ft2ac, core$BA_OBS_CURR))
cs <- rbind(mk("zero growth", core$Cords_o*0+core$Cords_o, core$Cords_o),
            mk("AcadianGY MORTCAL", core$Cords_PRED_ac, core$Cords_o))
bd <- rbind(mk("zero growth", core$BdFt_o, core$BdFt_o),
            mk("AcadianGY MORTCAL", core$BdFt_intl_PRED_ac, core$BdFt_o))
cat("--- BA ---\n"); print(ba, row.names=FALSE)
cat("--- Cords ---\n"); print(cs, row.names=FALSE)
cat("--- BdFt Intl ¼ ---\n"); print(bd, row.names=FALSE)
write.csv(ba, file.path(od, "silc_cfi_long_BA_v16_routine.csv"), row.names=FALSE)
write.csv(cs, file.path(od, "silc_cfi_long_Cords_v16_routine.csv"), row.names=FALSE)
write.csv(bd, file.path(od, "silc_cfi_long_BdFt_v16_routine.csv"), row.names=FALSE)

# Scatter figure
png(file.path(od, "silc_cfi_long_v16_scatter.png"),
    width=2700, height=900, res=170)
par(mfrow=c(1,3), mar=c(4.5,4.6,3.4,1.0), mgp=c(2.7,0.6,0))
GOLD <- "#B8860B"
draw <- function(o, p, lab) {
  lim <- c(0, max(c(o,p),na.rm=TRUE)*1.05)
  plot(NA, xlim=lim, ylim=lim, xlab="Observed", ylab="AcadianGY MORTCAL pred",
       main=lab, las=1, font.main=2, cex.main=1.2)
  abline(0,1,lty=2,col="#888")
  points(o, p, col=GOLD, pch=17, cex=1.3)
  legend("topleft", legend=sprintf("n=%d, bias=%.1f%%, R²=%.2f",
         sum(is.finite(p)&is.finite(o)), bias_pct(p,o), r2(p,o)),
         bty="n", cex=1.05)
}
draw(core$BA_OBS_CURR,   core$BA_PRED_ft2ac,     "BA (ft²/ac)")
draw(core$Cords_o,       core$Cords_PRED_ac,     "Cords/ac")
draw(core$BdFt_o,        core$BdFt_intl_PRED_ac, "BdFt Intl ¼ /ac")
mtext(sprintf("SILC CFI v16 long-horizon scorecard (n=%d routine pairs, mean %.1f yr)",
              nrow(core), mean(core$PERIOD_YR)),
      side=3, line=-1.2, outer=TRUE, font=2, cex=1.05)
dev.off()
cat("\nwrote silc_cfi_long_v16_scatter.png\n")
