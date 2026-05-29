#!/usr/bin/env Rscript
# silc_cfi_long_scorecard.R
# =====================================================================
# Long-horizon multi-model scorecard. Compares each CFI plot's
# earliest-measurement to latest-measurement projection across
# all available models.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
SCRIB_TO_INTL <- 1.18

# === Load ===
ps  <- read.csv(file.path(od, "silc_cfi_longhorizon_pairs.csv"))
ag  <- read.csv(file.path(od, "silc_cfi_long_agy_results.csv"))
fvs <- read.csv(file.path(od, "silc_cfi_long_fvs_results.csv"))
osm <- read.csv(file.path(od, "silc_cfi_long_osm_results.csv"))
tr  <- read.csv(file.path(od, "TREE.csv"))
EXPF_CFI <- 5.0

# === Compute observed Cords and BdFt at year_curr (apples-to-apples) ===
calc_cords <- function(d, h, e) {
  ok <- d >= 4.5 & is.finite(d) & is.finite(h)
  sum(0.0025 * d[ok]^2 * h[ok] * 0.90 * e[ok]) / 79
}
calc_intl_bdft <- function(d, h, e) {
  ok <- d >= 9.0 & is.finite(d) & is.finite(h)
  sum(0.01 * d[ok]^2 * h[ok] * e[ok])
}
obs_curr <- function(p, y) {
  t <- tr[tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN; h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  if (any(miss))
    h[miss] <- pmax(6, 4.27 + 82 * (1 - exp(-0.04 * (d[miss] * 2.54))))
  list(BA = sum(0.005454 * d^2 * EXPF_CFI),
       Cords = calc_cords(d, h, rep(EXPF_CFI, length(d))),
       BdFt  = calc_intl_bdft(d, h, rep(EXPF_CFI, length(d))))
}
ps$BA_obs_curr_calc    <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$BA)
ps$Cords_obs_curr_calc <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$Cords)
ps$BdFt_obs_curr_calc  <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_CURR[i])$BdFt)
ps$BA_obs_prev_calc    <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_PREV[i])$BA)
ps$Cords_obs_prev_calc <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_PREV[i])$Cords)
ps$BdFt_obs_prev_calc  <- sapply(seq_len(nrow(ps)), function(i) obs_curr(ps$PLOT[i], ps$YEAR_PREV[i])$BdFt)

# FIA-prior PAI for long horizon: scale by years
prior_PAI <- c("Cedar"=0.40,"Hardwood"=0.65,"Mixedwood"=0.85,
                "Commercial Softwood"=0.93,"Other Softwood"=0.70,
                "Unclassifiable"=NA)
prior_BdFt_PAI <- c("Cedar"=15,"Hardwood"=35,"Mixedwood"=42,
                     "Commercial Softwood"=45,"Other Softwood"=30)
prior_cords_PAI <- c("Cedar"=0.10,"Hardwood"=0.18,"Mixedwood"=0.22,
                      "Commercial Softwood"=0.25,"Other Softwood"=0.18)
strata_map <- read.csv(file.path(od, "silc_cfi_plot_strata_map.csv"))
ft <- setNames(strata_map$forest_type, strata_map$PLOT)
ps$forest_type <- ft[as.character(ps$PLOT)]

ps$BA_pred_zero    <- ps$BA_obs_prev_calc
ps$BA_pred_prior   <- ps$BA_obs_prev_calc + prior_PAI[ps$forest_type] * ps$PERIOD_YR
ps$Cords_pred_zero <- ps$Cords_obs_prev_calc
ps$Cords_pred_prior<- ps$Cords_obs_prev_calc + prior_cords_PAI[ps$forest_type] * ps$PERIOD_YR
ps$BdFt_pred_zero  <- ps$BdFt_obs_prev_calc

# === Combine model predictions ===
ag2 <- ag[, c("PLOT","YEAR_PREV","BA_PRED_ft2ac","Cords_PRED_ac","BdFt_intl_PRED_ac")]
names(ag2) <- c("PLOT","YEAR_PREV","BA_agm","Cords_agm","BdFt_agm")
fvs$model <- paste(fvs$variant, fvs$config, sep="_")
fvs$BdFt_intl <- fvs$BdFt_PRED * SCRIB_TO_INTL
fvs$Cords <- fvs$MCuFt_PRED / 79
fvs_w <- reshape(fvs[, c("PLOT","YEAR_PREV","model","BA_PRED_ft2ac","Cords","BdFt_intl")],
                 idvar=c("PLOT","YEAR_PREV"), timevar="model", direction="wide")

# OSM-ACD predicted tree list for Cords/BdFt — use Stand Summary only here
# (predicted BA only; cords/BdFt from same formula on osm tree list later)
osm2 <- osm[, c("PLOT","YEAR_PREV","BA_PRED_ft2ac")]
names(osm2) <- c("PLOT","YEAR_PREV","BA_osm")

m <- merge(ps, ag2, by=c("PLOT","YEAR_PREV"), all.x=TRUE)
m <- merge(m, fvs_w, by=c("PLOT","YEAR_PREV"), all.x=TRUE)
m <- merge(m, osm2,  by=c("PLOT","YEAR_PREV"), all.x=TRUE)
write.csv(m, file.path(od, "silc_cfi_long_predictions.csv"), row.names=FALSE)

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse <- function(p, o) sqrt(mean((p-o)^2, na.rm=TRUE))
r2 <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); if (sum(ok)<3) return(NA_real_)
  1 - sum((o[ok]-p[ok])^2) / sum((o[ok]-mean(o[ok]))^2)
}
mk <- function(lab, p, o) data.frame(predictor=lab, n=sum(is.finite(p)&is.finite(o)),
                                      bias_pct=round(bias_pct(p,o),2),
                                      RMSE=round(rmse(p,o),2),
                                      R2=round(r2(p,o),3))

# Routine subset: exclude establishment-driven plots (BA growth > 2x = ingrowth surge)
m$establishment <- (m$BA_CURR_FT2AC / pmax(m$BA_PREV_FT2AC, 1)) > 2.0
core <- m[!m$establishment, ]

cat("=== Long-horizon CFI scorecard FULL (n=", nrow(m), ", mean horizon ",
    round(mean(m$PERIOD_YR),1), " yr) ===\n", sep="")
cat("    Establishment pairs excluded from routine subset: ", sum(m$establishment),
    " (plots ", paste(m$PLOT[m$establishment], collapse=", "), ")\n", sep="")
cat("\n--- BA (ft^2/ac) ---\n")
ba <- rbind(
  mk("zero growth",       m$BA_pred_zero,                 m$BA_CURR_FT2AC),
  mk("FIA prior PAI",     m$BA_pred_prior,                m$BA_CURR_FT2AC),
  mk("AGM / AcadianGY",   m$BA_agm,                       m$BA_CURR_FT2AC),
  mk("FVS-NE default",    m$BA_PRED_ft2ac.NE_default,     m$BA_CURR_FT2AC),
  mk("FVS-NE calibrated", m$BA_PRED_ft2ac.NE_calibrated,  m$BA_CURR_FT2AC),
  mk("FVS-ACD default",   m$BA_PRED_ft2ac.ACD_default,    m$BA_CURR_FT2AC),
  mk("OSM-ACD",           m$BA_osm,                       m$BA_CURR_FT2AC)
)
print(ba, row.names=FALSE)

cat("\n--- Cords/ac ---\n")
cs <- rbind(
  mk("zero growth",       m$Cords_pred_zero,              m$Cords_obs_curr_calc),
  mk("FIA prior PAI",     m$Cords_pred_prior,             m$Cords_obs_curr_calc),
  mk("AGM / AcadianGY",   m$Cords_agm,                    m$Cords_obs_curr_calc),
  mk("FVS-NE default",    m$Cords.NE_default,             m$Cords_obs_curr_calc),
  mk("FVS-NE calibrated", m$Cords.NE_calibrated,          m$Cords_obs_curr_calc),
  mk("FVS-ACD default",   m$Cords.ACD_default,            m$Cords_obs_curr_calc)
)
print(cs, row.names=FALSE)

cat("\n--- BdFt Intl 1/4 /ac ---\n")
bd <- rbind(
  mk("zero growth",       m$BdFt_pred_zero,                m$BdFt_obs_curr_calc),
  mk("AGM / AcadianGY",   m$BdFt_agm,                      m$BdFt_obs_curr_calc),
  mk("FVS-NE default",    m$BdFt_intl.NE_default,          m$BdFt_obs_curr_calc),
  mk("FVS-NE calibrated", m$BdFt_intl.NE_calibrated,       m$BdFt_obs_curr_calc),
  mk("FVS-ACD default",   m$BdFt_intl.ACD_default,         m$BdFt_obs_curr_calc)
)
print(bd, row.names=FALSE)

write.csv(ba, file.path(od, "silc_cfi_long_BA.csv"), row.names=FALSE)
write.csv(cs, file.path(od, "silc_cfi_long_Cords.csv"), row.names=FALSE)
write.csv(bd, file.path(od, "silc_cfi_long_BdFt.csv"), row.names=FALSE)

cat("\n\n=== Long-horizon ROUTINE-GROWTH subset (n=", nrow(core), ") ===\n", sep="")
cat("\n--- BA (ft^2/ac) ---\n")
ba_c <- rbind(
  mk("zero growth",       core$BA_pred_zero,                 core$BA_CURR_FT2AC),
  mk("FIA prior PAI",     core$BA_pred_prior,                core$BA_CURR_FT2AC),
  mk("AGM / AcadianGY",   core$BA_agm,                       core$BA_CURR_FT2AC),
  mk("FVS-NE default",    core$BA_PRED_ft2ac.NE_default,     core$BA_CURR_FT2AC),
  mk("FVS-NE calibrated", core$BA_PRED_ft2ac.NE_calibrated,  core$BA_CURR_FT2AC),
  mk("FVS-ACD default",   core$BA_PRED_ft2ac.ACD_default,    core$BA_CURR_FT2AC),
  mk("OSM-ACD",           core$BA_osm,                       core$BA_CURR_FT2AC)
)
print(ba_c, row.names=FALSE)

cat("\n--- Cords/ac ---\n")
cs_c <- rbind(
  mk("zero growth",       core$Cords_pred_zero,              core$Cords_obs_curr_calc),
  mk("FIA prior PAI",     core$Cords_pred_prior,             core$Cords_obs_curr_calc),
  mk("AGM / AcadianGY",   core$Cords_agm,                    core$Cords_obs_curr_calc),
  mk("FVS-NE default",    core$Cords.NE_default,             core$Cords_obs_curr_calc),
  mk("FVS-NE calibrated", core$Cords.NE_calibrated,          core$Cords_obs_curr_calc),
  mk("FVS-ACD default",   core$Cords.ACD_default,            core$Cords_obs_curr_calc)
)
print(cs_c, row.names=FALSE)

cat("\n--- BdFt Intl 1/4 /ac ---\n")
bd_c <- rbind(
  mk("zero growth",       core$BdFt_pred_zero,                core$BdFt_obs_curr_calc),
  mk("AGM / AcadianGY",   core$BdFt_agm,                      core$BdFt_obs_curr_calc),
  mk("FVS-NE default",    core$BdFt_intl.NE_default,          core$BdFt_obs_curr_calc),
  mk("FVS-NE calibrated", core$BdFt_intl.NE_calibrated,       core$BdFt_obs_curr_calc),
  mk("FVS-ACD default",   core$BdFt_intl.ACD_default,         core$BdFt_obs_curr_calc)
)
print(bd_c, row.names=FALSE)

write.csv(ba_c, file.path(od, "silc_cfi_long_BA_routine.csv"), row.names=FALSE)
write.csv(cs_c, file.path(od, "silc_cfi_long_Cords_routine.csv"), row.names=FALSE)
write.csv(bd_c, file.path(od, "silc_cfi_long_BdFt_routine.csv"), row.names=FALSE)

# === Figure: BA / Cords / BdFt 3-panel scatter for AGM, FVS-NE cal, OSM-ACD ===
png(file.path(od, "silc_cfi_long_scatter.png"),
    width=2700, height=1100, res=170)
par(mfrow=c(1,3), mar=c(4.5, 4.6, 3.4, 1.0), mgp=c(2.7, 0.6, 0))
CRSF_GREEN <- "#1A3D28"; OSM_BLUE <- "#2c7fb8"; FVS_BLUE <- "#7DB5D5"

draw_multi <- function(o_col, models, ylab) {
  obs <- o_col
  lim <- c(0, max(c(obs, unlist(lapply(models, function(x) x$y))),
                  na.rm=TRUE) * 1.05)
  plot(NA, xlim=lim, ylim=lim, xlab="Observed", ylab="Predicted",
       main=ylab, las=1, font.main=2, cex.main=1.2)
  abline(0, 1, lty=2, col="#888")
  for (md in models) {
    points(obs, md$y, col=md$col, pch=md$pch, cex=1.5)
  }
  legend("topleft",
         legend = sapply(models, function(x) sprintf("%s (%.0f%%, R^2 %.2f)",
                                                       x$name, bias_pct(x$y, obs),
                                                       r2(x$y, obs))),
         col = sapply(models, function(x) x$col),
         pch = sapply(models, function(x) x$pch), bty="n", cex=0.95)
}

# Use routine subset for the scatter
draw_multi(core$BA_CURR_FT2AC, list(
  list(y=core$BA_agm, col=CRSF_GREEN, pch=19, name="AGM"),
  list(y=core$BA_PRED_ft2ac.NE_calibrated, col=FVS_BLUE, pch=17, name="FVS-NE cal"),
  list(y=core$BA_osm, col=OSM_BLUE, pch=15, name="OSM-ACD")), "BA (ft^2/ac)")
draw_multi(core$Cords_obs_curr_calc, list(
  list(y=core$Cords_agm, col=CRSF_GREEN, pch=19, name="AGM"),
  list(y=core$Cords.NE_calibrated, col=FVS_BLUE, pch=17, name="FVS-NE cal"),
  list(y=core$Cords.ACD_default, col=OSM_BLUE, pch=15, name="FVS-ACD def")), "Cords/ac")
draw_multi(core$BdFt_obs_curr_calc, list(
  list(y=core$BdFt_agm, col=CRSF_GREEN, pch=19, name="AGM"),
  list(y=core$BdFt_intl.NE_calibrated, col=FVS_BLUE, pch=17, name="FVS-NE cal"),
  list(y=core$BdFt_intl.ACD_default, col=OSM_BLUE, pch=15, name="FVS-ACD def")), "BdFt Intl 1/4/ac")
mtext(sprintf("SILC CFI long-horizon scorecard (n=%d routine pairs, %d-%d yr horizons, mean %.1f yr)",
              nrow(core), min(core$PERIOD_YR), max(core$PERIOD_YR), mean(core$PERIOD_YR)),
      side=3, line=-1.2, outer=TRUE, cex=1.0, font=2)
dev.off()
cat("\nwrote silc_cfi_long_scatter.png\n")
