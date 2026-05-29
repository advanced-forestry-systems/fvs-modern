#!/usr/bin/env Rscript
# silc_cfi_multimodel_scorecard.R
# =====================================================================
# Combine all available CFI model predictions into one scorecard:
#   * naive baselines (zero growth, FIA prior PAI)  -- from pair_summary
#   * AGM / AcadianGY 12.3.9 (in-source mortality + ingrowth fix)
#   * OSM-ACD                (Open Stand Model Acadian v2.26.1)
#
# Apples-to-apples volume formula across models:
#   per tree:   merch_cuft = 0.0025 * DBH_in^2 * HT_ft  if DBH >= 4.5
#               saw_cuft   = 0.55 * total_cuft         if DBH >= 9.0
#   cords/ac  = sum(merch_cuft * EXPF_per_ac) / 79
#   bdft/ac   = sum(saw_cuft   * EXPF_per_ac) * 6.0
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"

# === Load ===
ps <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"))
ag <- read.csv(file.path(od, "silc_cfi_acadiangy_pred_v3.csv"))
os <- read.csv(file.path(od, "silc_cfi_osm_results.csv"))
osm_tree <- read.csv(file.path(od, "silc_cfi_TreeListProjections.csv"))
fvs <- read.csv(file.path(od, "silc_cfi_fvs_results.csv"))

# Pivot FVS to wide: one row per pair, columns per variant/config
fvs$model <- paste(fvs$variant, fvs$config, sep="_")
fvs_keep <- fvs[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","model",
                    "BA_PRED_ft2ac","TPA_PRED","QMD_PRED_in","MCuFt_PRED")]
fvs_wide <- reshape(fvs_keep, idvar=c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"),
                    timevar="model", direction="wide")
# FVS Cords ~ MCuFt / 79
for (m in unique(fvs$model)) {
  c0 <- paste0("MCuFt_PRED.", m)
  fvs_wide[[paste0("Cords_", m)]] <- fvs_wide[[c0]] / 79
}

# Merch volume helper (imperial, per tree DBH in inches, HT in feet)
calc_vol <- function(dbh_in, ht_ft, expf_per_ac) {
  live <- dbh_in >= 4.5 & is.finite(dbh_in) & is.finite(ht_ft)
  d <- dbh_in[live]; h <- ht_ft[live]; e <- expf_per_ac[live]
  if (length(d) == 0) return(list(BA=0, Cords=0, BdFt=0, TPA=0))
  tcuft <- 0.0025 * d^2 * h
  merch <- tcuft * 0.90
  saw_ok<- d >= 9.0
  saw   <- ifelse(saw_ok, tcuft * 0.55, 0)
  list(BA    = sum(0.005454 * d^2 * e),
       TPA   = sum(e),
       Cords = sum(merch * e) / 79,
       BdFt  = sum(saw   * e) * 6.0)
}

# === Observed at year_curr (from TREE.csv, same uniform formula) ===
tr <- read.csv(file.path(od, "TREE.csv"))
EXPF_CFI <- 5.0
obs_curr <- function(p, y) {
  t <- tr[tr$PLOT == p & tr$MEASYEAR == y &
          tr$STATUSCD == 1 & is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN; h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  if (any(miss))
    h[miss] <- pmax(6, 4.27 + 82 * (1 - exp(-0.04 * (d[miss] * 2.54))))
  calc_vol(d, h, rep(EXPF_CFI, length(d)))
}

# === Compute OSM cords + BdFt from TreeListProjections ===
# OSM is metric: DBH in cm, HT in m, Stems in trees/ha. Convert per tree.
# Match each (SurveyID, year_off) to the manifest interval.
ACRES_PER_HA <- 2.4710538147
osm_tree$DBH_in <- osm_tree$DBH / 2.54
osm_tree$HT_ft  <- osm_tree$HT  / 0.3048
osm_tree$EXPF_ac<- osm_tree$Stems / ACRES_PER_HA

# Match SurveyID to (PLOT, YEAR_PREV) via sequential mapping (same order as build_inputs)
sid_to_pair <- ps[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR")]
sid_to_pair$SurveyID <- seq_len(nrow(sid_to_pair))

# Each pair: take trees at the row where (Year - first_year_in_pair) == PERIOD_YR
osm_pair_vol <- data.frame()
for (i in seq_len(nrow(sid_to_pair))) {
  pp <- sid_to_pair[i, ]
  sub <- osm_tree[osm_tree$SurveyID == pp$SurveyID, ]
  if (nrow(sub) == 0) next
  y0 <- min(sub$Year)
  sub$yr_off <- sub$Year - y0
  # Live trees only: Died is numeric 0 (not yet dead), Cut is string "False"
  died_ok <- is.na(sub$Died) | as.numeric(sub$Died) == 0
  cut_ok  <- is.na(sub$Cut)  | sub$Cut %in% c("False","false","FALSE","0",FALSE)
  sub_t <- sub[sub$yr_off == pp$PERIOD_YR &
               died_ok & cut_ok &
               is.finite(sub$DBH_in) & sub$DBH_in >= 4.5, ]
  if (nrow(sub_t) == 0) next
  v <- calc_vol(sub_t$DBH_in, sub_t$HT_ft, sub_t$EXPF_ac)
  osm_pair_vol <- rbind(osm_pair_vol, data.frame(
    PLOT = pp$PLOT, YEAR_PREV = pp$YEAR_PREV, YEAR_CURR = pp$YEAR_CURR,
    OSM_BA_ft2ac = v$BA, OSM_TPA = v$TPA,
    OSM_Cords_ac = v$Cords, OSM_BdFt_ac = v$BdFt
  ))
}

# === Build the combined panel ===
m <- merge(ps, ag, by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"),
           suffixes = c("","_agy"))
m <- merge(m, osm_pair_vol,
           by = c("PLOT","YEAR_PREV","YEAR_CURR"), all.x = TRUE)
m <- merge(m, fvs_wide,
           by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"), all.x = TRUE)

# Observed cords/BdFt at year_curr
o_curr <- t(sapply(seq_len(nrow(m)),
                   function(i) {
                     oc <- obs_curr(m$PLOT[i], m$YEAR_CURR[i])
                     c(BA=oc$BA, Cords=oc$Cords, BdFt=oc$BdFt)
                   }))
m$BA_obs_curr_calc    <- o_curr[, "BA"]
m$Cords_obs_curr_calc <- o_curr[, "Cords"]
m$BdFt_obs_curr_calc  <- o_curr[, "BdFt"]

# Establishment exclusion
m$establishment <- m$BA_PREV_FT2AC < 10 | abs(m$PAI_NET_OBS) > 5
core <- m[!m$establishment, ]

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse     <- function(p, o) sqrt(mean((p - o)^2, na.rm=TRUE))
r2       <- function(p, o) {
  ok <- is.finite(p) & is.finite(o); if (sum(ok)<3) return(NA_real_)
  1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2)
}

mk <- function(label, p, o, unit) {
  data.frame(predictor = label,
             n         = sum(is.finite(p) & is.finite(o)),
             bias_pct  = round(bias_pct(p, o), 2),
             RMSE      = round(rmse(p, o), 2),
             unit      = unit,
             R2        = round(r2(p, o), 3))
}

cat("=== Routine-growth subset (n=", nrow(core), ") ===\n", sep = "")
cat("\n--- BA (ft^2/ac) ---\n")
ba_sc <- rbind(
  mk("zero growth",       core$BA_pred_zero_growth,       core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("FIA prior PAI",     core$BA_pred_FIA_prior,         core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("AGM / AcadianGY",   core$BA_PRED_ft2ac,             core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("FVS-NE default",    core$`BA_PRED_ft2ac.NE_default`,    core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("FVS-NE calibrated", core$`BA_PRED_ft2ac.NE_calibrated`, core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("FVS-ACD default",   core$`BA_PRED_ft2ac.ACD_default`,   core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("FVS-ACD calibrated",core$`BA_PRED_ft2ac.ACD_calibrated`,core$BA_CURR_FT2AC, "ft^2/ac"),
  mk("OSM-ACD",           core$OSM_BA_ft2ac,              core$BA_CURR_FT2AC, "ft^2/ac")
)
print(ba_sc, row.names = FALSE)

cp <- read.csv(file.path(od, "silc_cfi_merch_predictions.csv"))
cp_lookup <- setNames(cp$Cords_zero,     paste(cp$PLOT, cp$YEAR_PREV))
core$Cords_zero <- cp_lookup[paste(core$PLOT, core$YEAR_PREV)]
cp_lookup <- setNames(cp$Cords_FIAprior, paste(cp$PLOT, cp$YEAR_PREV))
core$Cords_FIAprior <- cp_lookup[paste(core$PLOT, core$YEAR_PREV)]
cp_lookup <- setNames(cp$BdFt_zero,      paste(cp$PLOT, cp$YEAR_PREV))
core$BdFt_zero <- cp_lookup[paste(core$PLOT, core$YEAR_PREV)]

cat("\n--- Merch cords/ac ---\n")
cords_sc <- rbind(
  mk("zero growth",       core$Cords_zero,                core$Cords_obs_curr_calc, "cords/ac"),
  mk("FIA prior PAI",     core$Cords_FIAprior,            core$Cords_obs_curr_calc, "cords/ac"),
  mk("AGM / AcadianGY",   core$Cords_PRED_ac,             core$Cords_obs_curr_calc, "cords/ac"),
  mk("FVS-NE default",    core$Cords_NE_default,          core$Cords_obs_curr_calc, "cords/ac"),
  mk("FVS-NE calibrated", core$Cords_NE_calibrated,       core$Cords_obs_curr_calc, "cords/ac"),
  mk("FVS-ACD default",   core$Cords_ACD_default,         core$Cords_obs_curr_calc, "cords/ac"),
  mk("FVS-ACD calibrated",core$Cords_ACD_calibrated,      core$Cords_obs_curr_calc, "cords/ac"),
  mk("OSM-ACD",           core$OSM_Cords_ac,              core$Cords_obs_curr_calc, "cords/ac")
)
print(cords_sc, row.names = FALSE)

cat("\n--- Sawlog BdFt/ac ---\n")
bdft_sc <- rbind(
  mk("zero growth",       core$BdFt_zero,    core$BdFt_obs_curr_calc, "BdFt/ac"),
  mk("AGM / AcadianGY",   core$BdFt_PRED_ac, core$BdFt_obs_curr_calc, "BdFt/ac"),
  mk("OSM-ACD",           core$OSM_BdFt_ac,  core$BdFt_obs_curr_calc, "BdFt/ac")
)
print(bdft_sc, row.names = FALSE)

# (Cords_zero / Cords_FIAprior columns already present via merge from ag)

write.csv(core, file.path(od, "silc_cfi_multimodel_predictions.csv"), row.names = FALSE)
write.csv(ba_sc,    file.path(od, "silc_cfi_multimodel_BA.csv"),    row.names = FALSE)
write.csv(cords_sc, file.path(od, "silc_cfi_multimodel_Cords.csv"), row.names = FALSE)
write.csv(bdft_sc,  file.path(od, "silc_cfi_multimodel_BdFt.csv"),  row.names = FALSE)

# === Figure: 3 metrics x 2 named models (AGM, OSM) scatter ===
png(file.path(od, "silc_cfi_multimodel_scatter.png"),
    width = 2700, height = 1600, res = 165)
par(mfrow = c(2, 3), mar = c(4.5, 4.6, 3.4, 1.0), mgp = c(2.7, 0.6, 0))

CRSF_GREEN <- "#1A3D28"; OSM_BLUE <- "#2c7fb8"

draw <- function(pred, obs, label, col_main, unit) {
  ok <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 2) {
    plot.new(); title(main = sprintf("%s\n(insufficient data)", label))
    return()
  }
  lim <- c(0, max(c(pred[ok], obs[ok])) * 1.05)
  plot(obs, pred, pch = 19, col = col_main, cex = 1.4,
       xlim = lim, ylim = lim,
       xlab = sprintf("Observed (%s)", unit),
       ylab = sprintf("Predicted (%s)", unit),
       main = label, las = 1, font.main = 2, cex.main = 1.15)
  abline(0, 1, lty = 2, col = "#888")
  abline(lm(pred[ok] ~ obs[ok]), col = col_main, lwd = 1.8)
  bp <- bias_pct(pred, obs); rm <- rmse(pred, obs); r2v <- r2(pred, obs)
  legend("topleft", legend = c(
    sprintf("bias %+.1f%%", bp),
    sprintf("RMSE %.2f", rm),
    sprintf("R^2 %.2f", r2v),
    sprintf("n %d", sum(ok))), bty = "n", cex = 1.0)
}

draw(core$`BA_PRED_ft2ac.NE_calibrated`, core$BA_CURR_FT2AC,
     "FVS-NE calibrated: BA",   "#7DB5D5", "ft^2/ac")
draw(core$Cords_NE_calibrated, core$Cords_obs_curr_calc,
     "FVS-NE calibrated: cords","#7DB5D5", "cords/ac")
draw(core$BA_PRED_ft2ac, core$BA_CURR_FT2AC,
     "AGM / AcadianGY: BA",      CRSF_GREEN, "ft^2/ac")
draw(core$Cords_PRED_ac, core$Cords_obs_curr_calc,
     "AGM / AcadianGY: cords",   CRSF_GREEN, "cords/ac")
draw(core$OSM_BA_ft2ac,  core$BA_CURR_FT2AC,
     "OSM-ACD: BA",              OSM_BLUE, "ft^2/ac")
draw(core$OSM_Cords_ac,  core$Cords_obs_curr_calc,
     "OSM-ACD: cords",           OSM_BLUE, "cords/ac")
dev.off()
cat("\nwrote silc_cfi_multimodel_scatter.png\n")
cat("wrote silc_cfi_multimodel_*.csv\n")
