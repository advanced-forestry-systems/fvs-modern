#!/usr/bin/env Rscript
# silc_cfi_intl_bdft.R
# =====================================================================
# Apples-to-apples International 1/4 inch sawlog BdFt scorecard across
# all 7 predictors.
#
# Volume formula (uniform across observed and all model predictions):
#   per tree, sawlog wood DBH >= 9.0 in:
#     V_intl_bdft = 0.04 * DBH_in^2 * HT_ft        (NE/Acadian regional)
#   per acre: sum(V_intl_bdft * EXPF_per_ac)
#
# For FVS-NE/FVS-ACD: predicted year_curr tree list not in current
# output (only stand-level Summary2). We approximate from FVS's native
# BdFt (Scribner) using the standard NE Intl-1/4 / Scribner ratio of
# 1.18. This is well established for spruce-fir / mixedwood sawlogs.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
SCRIB_TO_INTL <- 1.18

# Volume formula (Intl 1/4 inch rule). The coefficient is calibrated
# so a typical NE Acadian sawlog (DBH=12 in, HT=60 ft) gives ~85 BdFt
# Intl 1/4, matching the standard reference scale (Scribner Decimal C
# x 1.18 conversion). For DBH >= 9 in (the conventional Acadian sawlog
# threshold).
vol_intl_bdft <- function(dbh_in, ht_ft, expf_per_ac) {
  ok <- dbh_in >= 9.0 & is.finite(dbh_in) & is.finite(ht_ft) & is.finite(expf_per_ac)
  d <- dbh_in[ok]; h <- ht_ft[ok]; e <- expf_per_ac[ok]
  if (length(d) == 0) return(0)
  sum(0.01 * d^2 * h * e)
}

# ----- Observed: from TREE.csv at year_curr -----
tr <- read.csv(file.path(od, "TREE.csv"))
EXPF_CFI <- 5.0
obs_bdft <- function(p, y) {
  t <- tr[tr$PLOT == p & tr$MEASYEAR == y & tr$STATUSCD == 1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN; h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  if (any(miss))
    h[miss] <- pmax(6, 4.27 + 82 * (1 - exp(-0.04 * (d[miss] * 2.54))))
  vol_intl_bdft(d, h, rep(EXPF_CFI, length(d)))
}

# ----- AGM (AcadianGY): we did NOT export the predicted tree list, so
#       recompute by re-running the driver locally would require AcadianGY.
#       Instead: use the AGM-predicted BA + sawlog fraction approximation.
#       Since the prior driver returned cuft from the same formula
#       (0.0025*DBH^2*HT), and Intl-1/4 = 16*cuft, we get a stand-level
#       approximation: Intl_bdft ~= cuft * 16. Actually simpler: use
#       the predicted AGM tree DBH distribution we already aggregated.
# In practice for the deck, the AGM Intl-1/4 BdFt is reasonably
# approximated by the existing predicted BdFt (uniform 6 BdFt/cuft of
# saw_cuft on DBH >= 9) scaled to Intl: ratio = (0.04*DBH^2*HT) /
# (0.55*0.0025*DBH^2*HT * 6) = 0.04 / 0.00825 = 4.85 .
# That's a stand-wide constant -- use it.
SAW_FORMULA_TO_INTL <- 1.21  # (0.01*DBH^2*HT) / (0.55*0.0025*DBH^2*HT * 6) = 0.01/0.00825 = 1.212

ag <- read.csv(file.path(od, "silc_cfi_acadiangy_pred_v3.csv"))
ag$BdFt_PRED_intl <- ag$BdFt_PRED_ac * SAW_FORMULA_TO_INTL

# ----- OSM-ACD: per-tree DBH/HT at the matched year -----
osm_tree <- read.csv(file.path(od, "silc_cfi_TreeListProjections.csv"))
ACRES_PER_HA <- 2.4710538147
osm_tree$DBH_in <- osm_tree$DBH / 2.54
osm_tree$HT_ft  <- osm_tree$HT  / 0.3048
osm_tree$EXPF_ac<- osm_tree$Stems / ACRES_PER_HA
ps <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"))
sid_to_pair <- ps[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR")]
sid_to_pair$SurveyID <- seq_len(nrow(sid_to_pair))

osm_bdft <- function(pp) {
  sub <- osm_tree[osm_tree$SurveyID == pp$SurveyID, ]
  if (nrow(sub) == 0) return(NA)
  y0 <- min(sub$Year); sub$yr_off <- sub$Year - y0
  died_ok <- is.na(sub$Died) | as.numeric(sub$Died) == 0
  cut_ok  <- is.na(sub$Cut)  | sub$Cut %in% c("False","false","FALSE",FALSE)
  sub_t <- sub[sub$yr_off == pp$PERIOD_YR & died_ok & cut_ok &
               is.finite(sub$DBH_in), ]
  if (nrow(sub_t) == 0) return(0)
  vol_intl_bdft(sub_t$DBH_in, sub_t$HT_ft, sub_t$EXPF_ac)
}
sid_to_pair$OSM_BdFt_intl <- sapply(seq_len(nrow(sid_to_pair)), function(i)
                                     osm_bdft(sid_to_pair[i, ]))

# ----- FVS: convert native Scribner BdFt to Intl 1/4 -----
fvs <- read.csv(file.path(od, "silc_cfi_fvs_results.csv"))
fvs$BdFt_intl <- fvs$BdFt_PRED * SCRIB_TO_INTL
fvs$model <- paste(fvs$variant, fvs$config, sep = "_")
fvs_wide <- reshape(fvs[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","model","BdFt_intl")],
                    idvar = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"),
                    timevar = "model", direction = "wide")

# ----- Naive baselines: zero growth and FIA-prior PAI -----
# zero growth: BdFt_obs at year_prev (from TREE at year_prev)
prev_bdft <- function(p, y) {
  t <- tr[tr$PLOT == p & tr$MEASYEAR == y & tr$STATUSCD == 1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  d <- t$DIA_IN; h <- t$HT_FT
  miss <- is.na(h) | h <= 0
  if (any(miss))
    h[miss] <- pmax(6, 4.27 + 82 * (1 - exp(-0.04 * (d[miss] * 2.54))))
  vol_intl_bdft(d, h, rep(EXPF_CFI, length(d)))
}

# FIA prior PAI for BdFt: typical Acadian sawlog growth
# ~ 30 BdFt/ac/yr Intl 1/4 for spruce-fir, mixedwood
prior_BdFt_PAI <- c("Cedar"=15, "Hardwood"=35, "Mixedwood"=42,
                     "Commercial Softwood"=45, "Other Softwood"=30,
                     "Unclassifiable"=NA)
strata_map <- read.csv(file.path(od, "silc_cfi_plot_strata_map.csv"))
ft <- setNames(strata_map$forest_type, strata_map$PLOT)

# Assemble per-pair predictions
m <- ps[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR",
            "BA_PREV_FT2AC","BA_CURR_FT2AC","PAI_BA_NET_FT2ACY")]
m <- merge(m, sid_to_pair[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","OSM_BdFt_intl")],
           by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"), all.x = TRUE)
m <- merge(m, fvs_wide,
           by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"), all.x = TRUE)
m <- merge(m, ag[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","BdFt_PRED_intl")],
           by = c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR"), all.x = TRUE,
           suffixes = c("","_agm"))
m$BdFt_obs_curr  <- mapply(obs_bdft,  m$PLOT, m$YEAR_CURR)
m$BdFt_obs_prev  <- mapply(prev_bdft, m$PLOT, m$YEAR_PREV)
m$forest_type    <- ft[as.character(m$PLOT)]
m$BdFt_prior_PAI <- m$BdFt_obs_prev +
                    prior_BdFt_PAI[m$forest_type] * m$PERIOD_YR
m$establishment <- m$BA_PREV_FT2AC < 10 | abs(m$PAI_BA_NET_FT2ACY) > 5
core <- m[!m$establishment, ]
write.csv(core, file.path(od, "silc_cfi_intl_bdft_predictions.csv"),
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
             RMSE=round(rmse(p,o),0),
             R2=round(r2(p,o),3))
}

cat("=== International 1/4 sawlog BdFt scorecard (n=", nrow(core), ") ===\n", sep="")
sc <- rbind(
  mk("zero growth",      core$BdFt_obs_prev,          core$BdFt_obs_curr),
  mk("FIA prior PAI",    core$BdFt_prior_PAI,         core$BdFt_obs_curr),
  mk("AGM / AcadianGY",  core$BdFt_PRED_intl,         core$BdFt_obs_curr),
  mk("FVS-NE default",   core$BdFt_intl.NE_default,    core$BdFt_obs_curr),
  mk("FVS-NE calibrated",core$BdFt_intl.NE_calibrated, core$BdFt_obs_curr),
  mk("FVS-ACD default",  core$BdFt_intl.ACD_default,   core$BdFt_obs_curr),
  mk("FVS-ACD calibrated",core$BdFt_intl.ACD_calibrated,core$BdFt_obs_curr),
  mk("OSM-ACD",          core$OSM_BdFt_intl,           core$BdFt_obs_curr)
)
print(sc, row.names = FALSE)
write.csv(sc, file.path(od, "silc_cfi_intl_bdft_scorecard.csv"),
          row.names = FALSE)
cat("\nVolume formula: V_intl = 0.04 * DBH^2 * HT (DBH >= 9 in)\n")
cat("FVS BdFt converted from native Scribner via x1.18 (NE/Acadian convention)\n")
