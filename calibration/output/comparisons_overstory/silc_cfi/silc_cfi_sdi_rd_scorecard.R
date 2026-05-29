#!/usr/bin/env Rscript
# silc_cfi_sdi_rd_scorecard.R
# =====================================================================
# Adds SDI (Reineke) and Curtis RD scorecards to the multi-model CFI
# benchmark. Both are computable directly from each model's predicted
# TPA and QMD at year_curr -- no new runs needed.
#
# Formulas:
#   SDI       = TPA * (QMD / 10)^1.605     (Reineke 1933)
#   Curtis RD = BA  / sqrt(QMD)             (NE convention, m2/ha not converted)
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"

sdi <- function(tpa, qmd) tpa * (qmd / 10)^1.605
curtis_rd <- function(ba, qmd) ba / sqrt(qmd)

# === Load ===
ps  <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"))
sm  <- read.csv(file.path(od, "STAND_METRICS.csv"))   # has observed REL_DENSITY + CURTIS_RD
ag  <- read.csv(file.path(od, "silc_cfi_acadiangy_pred_v3.csv"))
os  <- read.csv(file.path(od, "silc_cfi_osm_results.csv"))
fvs <- read.csv(file.path(od, "silc_cfi_fvs_results.csv"))

# Compute SDI + RD for each model from TPA / QMD
ag$SDI_PRED       <- sdi(ag$TPA_PRED,       ag$QMD_PRED_in)
ag$CurtisRD_PRED  <- curtis_rd(ag$BA_PRED_ft2ac, ag$QMD_PRED_in)
os$SDI_PRED       <- sdi(os$TPA_PRED,       os$QMD_PRED_in)
os$CurtisRD_PRED  <- curtis_rd(os$BA_PRED_ft2ac, os$QMD_PRED_in)
fvs$SDI_PRED      <- sdi(fvs$TPA_PRED,      fvs$QMD_PRED_in)
fvs$CurtisRD_PRED <- curtis_rd(fvs$BA_PRED_ft2ac, fvs$QMD_PRED_in)
fvs$model <- paste(fvs$variant, fvs$config, sep = "_")

fvs_w <- reshape(fvs[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","model",
                          "SDI_PRED","CurtisRD_PRED","TPA_PRED","QMD_PRED_in",
                          "BA_PRED_ft2ac")],
                  idvar = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"),
                  timevar = "model", direction = "wide")

# === Observed SDI and Curtis RD at year_curr (from STAND_METRICS) ===
sm_obs <- sm[sm$METRICS_RELIABLE == "Y", ]
sm_obs$SDI_obs <- sm_obs$SDI                              # already computed
sm_obs$CurtisRD_obs <- sm_obs$CURTIS_RD                  # already computed
sm_obs_curr <- sm_obs[, c("PLOT","MEASYEAR","SDI_obs","CurtisRD_obs")]

m <- merge(ps, ag, by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"),
           suffixes = c("","_agy"))
m <- merge(m, os[, c("PLOT","YEAR_PREV","YEAR_CURR","SDI_PRED","CurtisRD_PRED",
                     "BA_PRED_ft2ac")],
           by = c("PLOT","YEAR_PREV","YEAR_CURR"), all.x = TRUE,
           suffixes = c("","_osm"))
m <- merge(m, fvs_w,
           by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"), all.x = TRUE)
m <- merge(m, sm_obs_curr,
           by.x = c("PLOT","YEAR_CURR"), by.y = c("PLOT","MEASYEAR"),
           all.x = TRUE)

# Naive baselines: SDI / RD at year_prev
sm_prev <- sm_obs[, c("PLOT","MEASYEAR","SDI_obs","CurtisRD_obs")]
names(sm_prev)[3:4] <- c("SDI_prev","CurtisRD_prev")
m <- merge(m, sm_prev,
           by.x = c("PLOT","YEAR_PREV"), by.y = c("PLOT","MEASYEAR"),
           all.x = TRUE)

m$establishment <- m$BA_PREV_FT2AC < 10 | abs(m$PAI_NET_OBS) > 5
core <- m[!m$establishment, ]
write.csv(core, file.path(od, "silc_cfi_sdi_rd_predictions.csv"),
          row.names = FALSE)

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse     <- function(p, o) sqrt(mean((p - o)^2, na.rm=TRUE))
r2       <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); if (sum(ok) < 3) return(NA_real_)
  1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2)
}
mk <- function(label, p, o) {
  data.frame(predictor=label, n=sum(is.finite(p)&is.finite(o)),
             bias_pct=round(bias_pct(p,o),2),
             RMSE=round(rmse(p,o),2),
             R2=round(r2(p,o),3))
}

cat("=== SDI scorecard (n=", nrow(core), ") ===\n", sep="")
sdi_sc <- rbind(
  mk("zero growth",       core$SDI_prev,            core$SDI_obs),
  mk("AGM / AcadianGY",   core$SDI_PRED,            core$SDI_obs),
  mk("FVS-NE default",    core$SDI_PRED.NE_default,    core$SDI_obs),
  mk("FVS-NE calibrated", core$SDI_PRED.NE_calibrated, core$SDI_obs),
  mk("FVS-ACD default",   core$SDI_PRED.ACD_default,   core$SDI_obs),
  mk("FVS-ACD calibrated",core$SDI_PRED.ACD_calibrated,core$SDI_obs),
  mk("OSM-ACD",           core$SDI_PRED_osm,        core$SDI_obs)
)
print(sdi_sc, row.names = FALSE)

cat("\n=== Curtis RD scorecard (n=", nrow(core), ") ===\n", sep="")
rd_sc <- rbind(
  mk("zero growth",       core$CurtisRD_prev,             core$CurtisRD_obs),
  mk("AGM / AcadianGY",   core$CurtisRD_PRED,             core$CurtisRD_obs),
  mk("FVS-NE default",    core$CurtisRD_PRED.NE_default,    core$CurtisRD_obs),
  mk("FVS-NE calibrated", core$CurtisRD_PRED.NE_calibrated, core$CurtisRD_obs),
  mk("FVS-ACD default",   core$CurtisRD_PRED.ACD_default,   core$CurtisRD_obs),
  mk("FVS-ACD calibrated",core$CurtisRD_PRED.ACD_calibrated,core$CurtisRD_obs),
  mk("OSM-ACD",           core$CurtisRD_PRED_osm,         core$CurtisRD_obs)
)
print(rd_sc, row.names = FALSE)

write.csv(sdi_sc, file.path(od, "silc_cfi_sdi_scorecard.csv"), row.names = FALSE)
write.csv(rd_sc,  file.path(od, "silc_cfi_rd_scorecard.csv"),  row.names = FALSE)

cat("\n--- Formulas ---\n")
cat("SDI        = TPA * (QMD / 10)^1.605     (Reineke 1933)\n")
cat("Curtis RD  = BA  / sqrt(QMD)             (Curtis 1982 NE convention)\n")
