#!/usr/bin/env Rscript
# bakuzis_assess_eichhorn.R
# Patched copy of bakuzis_assess.R. Replaces the hardcoded line-114 Eichhorn UNTESTABLE flag with
# an actual score. Does NOT regenerate any projection; reads existing AGB-bearing metrics and the
# dominant heights derived by derive_htdom.py from the on-disk treelists.
#
# Eichhorn convention (documented, defensible default):
#   Eichhorn's rule = total stand aboveground biomass is a function of dominant height, largely
#   independent of site class. Dominant height = TPA-weighted mean height of the ~100 largest-DBH
#   trees per acre (Assmann). For each arm we fit AGB as a function of HT_DOM alone (base) and then
#   AGB as a function of HT_DOM interacted with site class (full). If dominant height alone explains
#   the AGB-height relationship well and adding site class reduces residual RMSE by < 10 percent,
#   the arm PASSes Eichhorn; otherwise it is FLAGged.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260629)
D    <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
SBV  <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
HTD  <- file.path(D, "eichhorn_htdom")   # derive_htdom.py outputs live here
OUT  <- D
elog <- file.path(OUT, "error_log_eichhorn.txt")
logerr <- function(w, e) cat(sprintf("[%s] %s: %s\n", Sys.time(), w, conditionMessage(e)),
                             file = elog, append = TRUE)

arms <- list(
  b1_native = list(metrics = file.path(D, "out_conus_eq",          "*conus_b1_metrics.csv"),
                   htd     = "b1_native"),
  b2_native = list(metrics = file.path(D, "out_conus_eq",          "*conus_b2_metrics.csv"),
                   htd     = "b2_native"),
  b2_gompit = list(metrics = file.path(D, "out_conus_eq_gompit",   "*conus_b2_gompit_metrics.csv"),
                   htd     = "b2_gompit"),
  b1_gompit = list(metrics = file.path(D, "out_conus_eq_b1_gompit","*conus_b1_gompit_metrics.csv"),
                   htd     = "b1_gompit"))

## ---- per-stand SITE_INDEX map + global tertiles (identical to bakuzis_assess.R) ----
si_map <- rbindlist(lapply(Sys.glob(file.path(SBV, "standinit_*.csv")), function(f) tryCatch({
  d <- fread(f, select = c("STAND_CN","SITE_INDEX"), colClasses = list(character = "STAND_CN"))
  d[is.finite(SITE_INDEX) & SITE_INDEX > 0]
}, error = function(e) { logerr(paste("si", f), e); NULL })), fill = TRUE)
si_map <- unique(si_map, by = "STAND_CN")
qt <- quantile(si_map$SITE_INDEX, c(1/3, 2/3), na.rm = TRUE)
si_map[, site_class := fifelse(SITE_INDEX <= qt[1], "low",
                       fifelse(SITE_INDEX <= qt[2], "mid", "high"))]
setkey(si_map, STAND_CN)
cat(sprintf("SITE_INDEX tertile cutpoints: %.1f / %.1f ; n stands=%d\n", qt[1], qt[2], nrow(si_map)))

rmse <- function(r) sqrt(mean(r^2))

score_rows <- list()
for (a in names(arms)) {
  tryCatch({
    # AGB per stand-year
    mfs <- Sys.glob(arms[[a]]$metrics)
    M <- rbindlist(lapply(mfs, function(f) fread(f,
             select = c("STAND_CN","PROJ_YEAR","AGB_TONS_AC"),
             colClasses = list(character = "STAND_CN"))), fill = TRUE)
    M <- M[is.finite(AGB_TONS_AC)]
    # dominant heights per stand-year (concatenate the derive_htdom outputs for this arm)
    hfs <- Sys.glob(file.path(HTD, paste0("htdom_", arms[[a]]$htd, "_*.csv")))
    H <- rbindlist(lapply(hfs, function(f) fread(f,
             select = c("STAND_CN","PROJ_YEAR","HT_DOM_M","HT_MEAN_M"),
             colClasses = list(character = "STAND_CN"))), fill = TRUE)
    H <- H[is.finite(HT_DOM_M) & HT_DOM_M > 0]
    # join AGB <- HT_DOM on (STAND_CN, PROJ_YEAR), then site_class on STAND_CN
    DT <- merge(M, H, by = c("STAND_CN","PROJ_YEAR"))
    DT <- DT[si_map[, .(STAND_CN, site_class)], on = "STAND_CN", nomatch = 0L]
    n_raw <- nrow(DT[is.finite(AGB_TONS_AC) & is.finite(HT_DOM_M) & HT_DOM_M > 0])
    # Physical-plausibility guard. The un-thinned b1_native arm over-densifies and emits a small
    # number of nonphysical stand-years (runaway AGB / height-diameter extrapolation) that would
    # otherwise dominate an ordinary-least-squares fit. Bounds are generous: no real North American
    # stand has dominant height > 80 m or AGB > 2000 tons ac-1; these thresholds drop only corrupt
    # rows, not real forest. n_drop is reported for transparency.
    DT <- DT[is.finite(AGB_TONS_AC) & is.finite(HT_DOM_M) &
             HT_DOM_M > 0 & HT_DOM_M <= 80 & AGB_TONS_AC > 0 & AGB_TONS_AC <= 2000]
    DT[, site_class := factor(site_class, levels = c("low","mid","high"))]
    n <- nrow(DT); n_drop <- n_raw - n
    # base: AGB ~ f(HT_DOM) only ; full: AGB ~ f(HT_DOM) x site_class
    mb <- lm(AGB_TONS_AC ~ HT_DOM_M + I(HT_DOM_M^2), data = DT)
    mf <- lm(AGB_TONS_AC ~ (HT_DOM_M + I(HT_DOM_M^2)) * site_class, data = DT)
    rb <- rmse(residuals(mb)); rf <- rmse(residuals(mf))
    r2b <- summary(mb)$r.squared
    pct_red <- 100 * (rb - rf) / rb
    # site-class band: max spread of predicted AGB across site classes at fixed HT_DOM,
    # relative to overall AGB sd (a scale-free tightness check)
    hgrid <- quantile(DT$HT_DOM_M, seq(0.1, 0.9, 0.1), na.rm = TRUE)
    pg <- CJ(HT_DOM_M = hgrid, site_class = factor(c("low","mid","high"),
                                                   levels = c("low","mid","high")))
    pg[, pred := predict(mf, newdata = pg)]
    band <- pg[, .(spread = max(pred) - min(pred)), by = HT_DOM_M]
    band_frac <- mean(band$spread) / sd(DT$AGB_TONS_AC)
    verdict <- if (is.finite(pct_red) && pct_red < 10) "PASS" else "FLAG"
    score_rows[[a]] <- data.table(
      arm = a, n_stand_year = n, n_dropped_nonphysical = n_drop,
      R2_AGB_HTDOM = round(r2b, 4),
      rmse_base_tons_ac = round(rb, 3),
      rmse_full_tons_ac = round(rf, 3),
      site_class_rmse_reduction_pct = round(pct_red, 2),
      site_class_band_frac_of_sd = round(band_frac, 3),
      eichhorn_flag = verdict)
    cat(sprintf("arm %s: n=%d drop=%d R2=%.3f rmse_base=%.2f rmse_full=%.2f red=%.2f%% band=%.3f -> %s\n",
                a, n, n_drop, r2b, rb, rf, pct_red, band_frac, verdict))
    rm(M, H, DT); gc()
  }, error = function(e) logerr(paste("score", a), e))
}
scores <- rbindlist(score_rows, fill = TRUE)
fwrite(scores, file.path(OUT, "eichhorn_scores.csv"))
print(scores)
cat("bakuzis_assess_eichhorn done\n")
