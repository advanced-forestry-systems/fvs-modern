## Site productivity head-to-head: BGI vs climate_si (CSI) vs CSPI v4 as
## predictors of (a) diameter growth, (b) height growth, (c) HT-DBH allometry,
## (d) FIA SICOND. Controls: ln(DBH1), BAL, sqrt(BA), species fixed effect.
## Subsamples for tractability. Reports incremental R^2 and coefficient sign/magnitude.
suppressMessages({library(data.table); library(bit64)})
CAL <- "/users/PUOM0008/crsfaaron/fvs-modern/calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
CSPI <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration/cspi_v4_at_calib_plots.csv"
OUT  <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

dat <- as.data.table(readRDS(CAL))
cat("nrow:", nrow(dat), "\n")
cat("site/resp cols present:", paste(intersect(c("bgi","climate_si","cspi","SICOND","SICOND_FVS","DBH1","DBH2","HT1","HT2","YEARS","BAL_SW1","BAL_HW1","BA1","SPCD"), names(dat)), collapse=", "), "\n")

## merge CSPI v4
lk <- fread(CSPI); dat[, .pid := as.character(PLT_CN_cond1)]; lk[, .pid := as.character(PLT_CN_cond1)]
dat[lk, cspi_v4 := i.cspi_v4, on=".pid"]
cat("cspi_v4 merged for", sum(is.finite(dat$cspi_v4)), "of", nrow(dat), "rows\n")

## derive responses
dat[, ln_dbh := log(pmax(DBH1, 2.54))]
dat[, bal := (BAL_SW1 + BAL_HW1)]
dat[, sqrt_ba := sqrt(pmax(BA1, 0))]
if ("DBH2" %in% names(dat) && "YEARS" %in% names(dat)) dat[, dg_ann := (DBH2-DBH1)/YEARS]
if ("HT2"  %in% names(dat) && "YEARS" %in% names(dat)) dat[, hg_ann := (HT2-HT1)/YEARS]
dat[, ht_log_above_bh := log(pmax(HT1-1.37, 0.01))]

## subsample for tractability
set.seed(42); idx <- sample.int(nrow(dat), min(200000, nrow(dat)))
d <- dat[idx]
d <- d[is.finite(ln_dbh) & is.finite(bal) & is.finite(sqrt_ba) & !is.na(SPCD)]
cat("subsample after clean:", nrow(d), "\n")

## ---- helper: incremental R^2 + coefficient for a productivity var ----
fit_inc <- function(resp_name, var_name) {
  d2 <- d[is.finite(get(resp_name)) & is.finite(get(var_name))]
  if (nrow(d2) < 1000) return(data.table(response=resp_name, var=var_name, n=nrow(d2), r2_base=NA, r2_full=NA, dR2=NA, coef=NA, se=NA, t=NA))
  f0 <- lm(as.formula(paste(resp_name, "~ ln_dbh + bal + sqrt_ba + factor(SPCD)")), data=d2)
  f1 <- lm(as.formula(paste(resp_name, "~ ln_dbh + bal + sqrt_ba + factor(SPCD) +", var_name)), data=d2)
  cf <- summary(f1)$coefficients[var_name, ]
  data.table(response=resp_name, var=var_name, n=nrow(d2),
             r2_base=round(summary(f0)$r.squared,4),
             r2_full=round(summary(f1)$r.squared,4),
             dR2=round(summary(f1)$r.squared - summary(f0)$r.squared, 5),
             coef=signif(cf["Estimate"],3), se=signif(cf["Std. Error"],3),
             t=signif(cf["Estimate"]/cf["Std. Error"],3))
}

resps <- intersect(c("dg_ann","hg_ann","ht_log_above_bh","SICOND"), names(d))
vars  <- intersect(c("bgi","climate_si","cspi_v4"), names(d))
cat("responses:", paste(resps,collapse=", "), "| vars:", paste(vars,collapse=", "), "\n")

out <- rbindlist(lapply(resps, function(r) rbindlist(lapply(vars, function(v) fit_inc(r,v)))))
cat("\n=== HEAD-TO-HEAD: incremental R^2 and coefficient per (response, productivity var) ===\n")
print(out)
fwrite(out, file.path(OUT,"site_var_headtohead_growth.csv"))

## correlation matrix among productivity vars
co <- function(x,y){ok<-is.finite(x)&is.finite(y); if(sum(ok)<100) NA else round(cor(x[ok],y[ok]),3)}
cat("\n=== correlation matrix among productivity vars (incl SICOND) ===\n")
vmat <- intersect(c("bgi","climate_si","cspi_v4","SICOND","SICOND_FVS"), names(d))
M <- outer(vmat, vmat, Vectorize(function(a,b) co(d[[a]], d[[b]])))
dimnames(M) <- list(vmat, vmat); print(M)
fwrite(as.data.table(M, keep.rownames="var"), file.path(OUT,"site_var_correlations_growth.csv"))
cat("\nDONE\n")
