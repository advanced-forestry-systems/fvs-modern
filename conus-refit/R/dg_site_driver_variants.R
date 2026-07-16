#!/usr/bin/env Rscript
# dg_site_driver_variants.R
# Fit a FAMILY of Greg-Johnson-form CONUS diameter-growth (DG) equations that
# differ only in the site / climate driver term, so the greg arm can offer them
# as FVS-style keyword-selectable options (like an FVS variant keyword):
#
#   base size+competition kernel (all variants):
#     eta0 = B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7)
#   driver term (swappable):  + B5 * DRIVER
#   DRIVER in { none, elev, bgi, cspi, esi(=climate_si) }   [+ emt, DF-only below]
#   dg_annual = exp(eta0 [+ B5*DRIVER])
#
# Single-step annual approximation (compounding ignored) applied identically to
# every variant, so cross-variant RMSE is a fair comparison. Fit on our CONUS
# remeasurement pairs (bgi/cspi/climate_si present for all species; EMT is not,
# so EMT stays a Douglas-fir reference via Greg's df_dg_res in Part A).
# Emits greg_dg_coefficients_<driver>.csv per variant + a held-out scorecard.
# Author: Cowork autopilot 2026-07-03. Seed 20260703.

suppressPackageStartupMessages({ library(data.table); library(minpack.lm) })
set.seed(20260703)
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/dg_site_variants"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CFG <- "/users/PUOM0008/crsfaaron/fvs-modern/config"
rmse <- function(a,b) sqrt(mean((a-b)^2)); bias <- function(a,b) mean(a-b)

drivers <- c("none","elev","bgi","cspi","esi")   # esi := climate_si
drv_col <- c(elev="ELEV", bgi="bgi", cspi="cspi", esi="esi")   # driver -> data column
base_start <- c(B0=-1.326, B1=-0.57, B2=-0.05, B3=1.5, B4=0.8)
fit_driver <- function(dat, drv) {
  if (drv == "none") {
    frm <- dg_obs ~ exp(B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7))
    st  <- as.list(base_start)
  } else {
    frm <- dg_obs ~ exp(B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7) + B5*drv_x)
    st  <- as.list(c(base_start, B5 = -1e-4))
  }
  tryCatch(nlsLM(frm, data=dat, start=st, control=nls.lm.control(maxiter=200)),
           error=function(e) NULL)
}

## ---- load our CONUS pairs ---------------------------------------------------
d <- as.data.table(readRDS("~/fvs-conus/calibration/data/conus_remeasurement_pairs.rds"))
d <- d[is.finite(DBH1)&is.finite(DBH2)&is.finite(CR1)&is.finite(HT1)&is.finite(BAL1)&
       is.finite(YEARS)&YEARS>0 & CR1>0 & CR1<1 & DBH1>=1 & STATUS1==1 & STATUS2==1]
d[, dg_obs := (DBH2-DBH1)/YEARS]
d <- d[dg_obs > -0.1 & dg_obs < 2]
setnames(d, c("DBH1","CR1","HT1","BAL1"), c("dbh","cr","ht","bal"))
d[, esi := climate_si]
cat("coverage (finite):",
    paste(sprintf("%s %.0f%%", c("elev","bgi","cspi","esi"),
      100*c(mean(is.finite(d$ELEV)),mean(is.finite(d$bgi)),mean(is.finite(d$cspi)),mean(is.finite(d$esi)))),
      collapse="  "), "\n")

sp_n <- d[, .N, by=SPCD][N >= 3000]
score <- list()
for (drv in drivers) {
  co <- list(); te_all <- list()
  for (s in sp_n$SPCD) {
    ds <- d[SPCD==s]
    if (drv != "none") { ds[, drv_x := get(drv_col[[drv]])]; ds <- ds[is.finite(drv_x)] }
    if (nrow(ds) < 2000) next
    i <- sample(nrow(ds), floor(nrow(ds)*0.7)); tr <- ds[i]; te <- ds[-i]
    m <- fit_driver(tr, drv); if (is.null(m)) next
    p <- coef(m)
    rowc <- data.table(SPCD=s, n=nrow(ds), B0=p["B0"],B1=p["B1"],B2=p["B2"],B3=p["B3"],B4=p["B4"],
                       B5=if("B5"%in%names(p)) p["B5"] else NA_real_)
    co[[as.character(s)]] <- rowc
    eta <- p["B0"]+p["B1"]*log((te$dbh+1)^2/(te$cr*te$ht+1)^p["B3"])+p["B2"]*te$bal^p["B4"]/log(te$dbh+2.7)
    if (drv!="none") eta <- eta + p["B5"]*te$drv_x
    te_all[[as.character(s)]] <- data.table(SPCD=s, obs=te$dg_obs, pred=exp(eta))
  }
  CO <- rbindlist(co); TE <- rbindlist(te_all)
  fwrite(CO, file.path(CFG, sprintf("greg_dg_coefficients_%s.csv", drv)))
  score[[drv]] <- data.table(driver=drv, species=nrow(CO), n_test=nrow(TE),
                             RMSE=rmse(TE$pred,TE$obs), bias=bias(TE$pred,TE$obs),
                             r=suppressWarnings(cor(TE$pred,TE$obs)))
}
SC <- rbindlist(score)[order(RMSE)]
cat("\n=== DG site-driver scorecard (held-out, in/yr; lower RMSE better) ===\n"); print(SC)
fwrite(SC, file.path(OUT, "dg_site_variant_scorecard.csv"))
cat("\nBenchmarks (six-arm held-out, in/yr): default 0.269 | recalibrated 0.134 | spp-specific 0.090 | spp-free 0.098\n")
cat("coefficient sets written to", CFG, "as greg_dg_coefficients_{none,elev,bgi,cspi,esi}.csv\n")
