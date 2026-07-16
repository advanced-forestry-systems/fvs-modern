#!/usr/bin/env Rscript
# dg_default_recal_per_variant.R
# Per-FVS-variant held-out DG RMSE/bias for the DG comparison arms that ARE
# reproducible offline from conus_remeasurement_pairs, using the SAME held-out
# split and filters as the reference six-arm harness (dg_site_driver_variants.R
# and dg_hook_ab_conus.R). Seed 20260703 to match the reference split.
#
# ARMS produced here (offline-reproducible):
#   (1) spp_specific : Greg-form DG kernel refit PER SPECIES on TRAIN, held-out
#                      TEST RMSE, grouped by fvs_variant. This is the offline
#                      "species-specific" benchmark (CONUS pool ~0.090 in/yr),
#                      now split by variant. Single-step annual form:
#                      dg = exp(B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3)
#                               + B2*bal^B4/log(dbh+2.7))
#   (2) spp_free     : ONE pooled Greg-form kernel fit on all TRAIN species,
#                      applied to all TEST, grouped by variant (CONUS pool ~0.098).
#   (3) null         : per-variant mean(dg_obs) predictor (skill floor).
#
# NOT reproducible offline (engine-native; reported as a gap, not fabricated):
#   default     : native FVS Wykoff ln(DDS) with the compiled per-variant/species
#                 coefficient vectors + internal SITE index/CCF state. Config JSON
#                 carries only growth$B1/B2 (a partial small-tree piece) + species
#                 crosswalk, NOT the full DDS coefficient set. Needs the engine.
#   recalibrated: native default DG * dds_multiplier (per-species Bayesian intercept
#                 shift). The dds_multiplier IS in config, but the native base it
#                 multiplies is engine-only, and dds_mapping is flagged "approximate
#                 (no crosswalk)". Needs the compiled engine.
suppressPackageStartupMessages({ library(data.table); library(minpack.lm) })
set.seed(20260703)
setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/dg_default_recal_per_variant"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rmse <- function(a,b) sqrt(mean((a-b)^2)); bias <- function(a,b) mean(a-b)

base_start <- c(B0=-1.326, B1=-0.57, B2=-0.05, B3=1.5, B4=0.8)
frm <- dg_obs ~ exp(B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7))
fit_kernel <- function(dat) tryCatch(
  nlsLM(frm, data=dat, start=as.list(base_start), control=nls.lm.control(maxiter=200)),
  error=function(e) NULL)
pred_kernel <- function(p, te)
  exp(p["B0"] + p["B1"]*log((te$dbh+1)^2/(te$cr*te$ht+1)^p["B3"]) +
      p["B2"]*te$bal^p["B4"]/log(te$dbh+2.7))

## ---- load + filter (exact mirror of the reference harness) ----
d <- as.data.table(readRDS("~/fvs-conus/calibration/data/conus_remeasurement_pairs.rds"))
d <- d[is.finite(DBH1)&is.finite(DBH2)&is.finite(CR1)&is.finite(HT1)&is.finite(BAL1)&
       is.finite(YEARS)&YEARS>0 & CR1>0 & CR1<1 & DBH1>=1 & STATUS1==1 & STATUS2==1]
d[, dg_obs := (DBH2-DBH1)/YEARS]
d <- d[dg_obs > -0.1 & dg_obs < 2]
setnames(d, c("DBH1","CR1","HT1","BAL1"), c("dbh","cr","ht","bal"))
d[, variant := toupper(fvs_variant)]
cat("rows:", format(nrow(d), big.mark=","), " variants:", uniqueN(d$variant),
    " species:", uniqueN(d$SPCD), "\n")

## ---- per-species 70/30 split (same seed => same split as reference) ----
d[, is_test := FALSE]
for (s in unique(d$SPCD)) {
  idx <- which(d$SPCD == s); if (!length(idx)) next
  te  <- sample(idx, length(idx) - floor(length(idx)*0.7))
  d[te, is_test := TRUE]
}
tr_all <- d[is_test == FALSE]; te_all <- d[is_test == TRUE]

## ---- (1) spp_specific : per-species kernel on TRAIN, predict TEST ----
sp_n <- d[, .N, by=SPCD][N >= 3000, SPCD]
te_all[, pred_spp := NA_real_]
for (s in sp_n) {
  tr <- tr_all[SPCD==s]; if (nrow(tr) < 1400) next
  m  <- fit_kernel(tr); if (is.null(m)) next
  ii <- which(te_all$SPCD==s); if (!length(ii)) next
  te_all[ii, pred_spp := pred_kernel(coef(m), te_all[ii])]
}

## ---- (2) spp_free : one pooled kernel on all TRAIN species with >=3000 obs ----
m_free <- fit_kernel(tr_all[SPCD %in% sp_n][sample(.N, min(.N, 3e5))])
if (!is.null(m_free)) te_all[SPCD %in% sp_n, pred_free := pred_kernel(coef(m_free), .SD)]

## ---- (3) per-variant scorecard over test rows with a fitted species ----
ev <- te_all[SPCD %in% sp_n & is.finite(pred_spp)]
sc <- ev[, .(
  n_test        = .N,
  n_species     = uniqueN(SPCD),
  RMSE_spp      = rmse(pred_spp,  dg_obs),
  bias_spp      = bias(pred_spp,  dg_obs),
  RMSE_sppfree  = if (all(is.finite(pred_free))) rmse(pred_free, dg_obs) else NA_real_,
  bias_sppfree  = if (all(is.finite(pred_free))) bias(pred_free, dg_obs) else NA_real_,
  RMSE_null     = rmse(rep(mean(dg_obs), .N), dg_obs)
), by = variant][order(-n_test)]

## CONUS pool (should reproduce the ~0.090 spp-specific / ~0.098 spp-free benchmark)
pool <- ev[, .(variant="CONUS_POOL", n_test=.N, n_species=uniqueN(SPCD),
               RMSE_spp=rmse(pred_spp,dg_obs), bias_spp=bias(pred_spp,dg_obs),
               RMSE_sppfree=rmse(pred_free,dg_obs), bias_sppfree=bias(pred_free,dg_obs),
               RMSE_null=rmse(rep(mean(dg_obs),.N),dg_obs))]
SC <- rbind(sc, pool)
cat("\n=== Per-variant offline-reproducible DG arms (held-out, in/yr) ===\n")
print(SC, nrow=100)
fwrite(SC, file.path(OUT, "dg_offline_arms_per_variant.csv"))
cat("\nBenchmarks (six-arm held-out, in/yr): default 0.269 | recalibrated 0.134 |",
    "spp-specific 0.090 | spp-free 0.098\n")
cat("wrote", file.path(OUT, "dg_offline_arms_per_variant.csv"), "\nDONE\n")
