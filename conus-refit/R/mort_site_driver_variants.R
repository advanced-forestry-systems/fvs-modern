#!/usr/bin/env Rscript
# mort_site_driver_variants.R
# Site-driver family for CONUS mortality, to offer keyword-selectable options
# alongside the DG family. Greg's base gompit uses only CR and CCH (no site
# term), so this asks: does adding a productivity driver improve survival
# discrimination? Fit as discrete-time survival (cloglog = gompit family) with a
# log(YEARS) exposure offset, so the annual hazard carries the driver:
#   cloglog(P_die_interval) = log(YEARS) + b0 + b1*lnDBH + b2*CR + b3*CCH + b5*Zdriver
# Drivers standardized (mean0/sd1) -> bounded, comparable, and fixes the BGI
# overflow seen in the DG family. Driver in {none,elev,bgi,cspi,esi}.
# Metrics: held-out AUC (discrimination) + log-loss. Author: Cowork autopilot 2026-07-03.

suppressPackageStartupMessages({ library(data.table) })
set.seed(20260703)
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/mort_site_variants"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CFG <- "/users/PUOM0008/crsfaaron/fvs-modern/config"
drivers <- c("none","elev","bgi","cspi","esi","emt")
drv_col <- c(elev="ELEV", bgi="bgi", cspi="cspi", esi="esi", emt="EMT")

auc <- function(score, y) {                       # fast rank AUC (drop NA preds)
  ok <- is.finite(score) & !is.na(y); score <- score[ok]; y <- y[ok]
  r <- rank(score); n1 <- sum(y==1); n0 <- sum(y==0)
  if (n1==0 || n0==0) return(NA_real_)
  (sum(r[y==1]) - n1*(n1+1)/2) / (n1*n0)
}
logloss <- function(p, y){ p <- pmin(pmax(p,1e-6),1-1e-6); -mean(y*log(p)+(1-y)*log(1-p)) }

d <- as.data.table(readRDS("~/fvs-conus/calibration/data/conus_remeasurement_pairs.rds"))
# join EMT climate from the P11 crosswalk (100% coverage)
xw <- fread("/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/emt_crosswalk/plot_key_emt_td_crosswalk.csv")
d <- merge(d, xw[, .(plot_key, EMT)], by="plot_key", all.x=TRUE)
# mortality response: alive (STATUS2==1) vs dead (STATUS2==2); drop harvested/other
d <- d[STATUS1==1 & STATUS2 %in% c(1,2) & is.finite(DBH1)&DBH1>=1 & CR1>0 & CR1<1 &
       is.finite(CCH1) & is.finite(YEARS) & YEARS>0]
d[, died := as.integer(STATUS2==2)]
d[, lnDBH := log(DBH1)][, esi := climate_si]
cat(sprintf("trees: %s  overall annualized-ish mortality frac: %.3f\n",
            format(nrow(d),big.mark=","), mean(d$died)))

sp_n <- d[, .N, by=SPCD][N >= 5000]
score <- list(); coefs <- list()
for (drv in drivers) {
  te_all <- list(); co <- list()
  for (s in sp_n$SPCD) {
    ds <- d[SPCD==s]
    if (drv!="none") { ds[, Z := scale(get(drv_col[[drv]]))[,1]]; ds <- ds[is.finite(Z)] }
    if (sum(ds$died) < 50 || nrow(ds) < 3000) next    # need enough deaths
    i <- sample(nrow(ds), floor(nrow(ds)*0.7)); tr <- ds[i]; te <- ds[-i]
    form <- if (drv=="none") died ~ lnDBH + CR1 + CCH1 else died ~ lnDBH + CR1 + CCH1 + Z
    m <- tryCatch(glm(form, data=tr, family=binomial("cloglog"),
                      offset=log(tr$YEARS)), error=function(e) NULL)
    if (is.null(m)) next
    pr <- tryCatch(predict(m, newdata=te, type="response"), error=function(e) NULL)
    if (is.null(pr)) next
    te_all[[as.character(s)]] <- data.table(SPCD=s, died=te$died, p=pr)
    cf <- coef(m); co[[as.character(s)]] <- data.table(SPCD=s, n=nrow(ds), deaths=sum(ds$died),
        b0=cf[1], b_lnDBH=cf["lnDBH"], b_CR=cf["CR1"], b_CCH=cf["CCH1"],
        b_driver=if("Z"%in%names(cf)) cf["Z"] else NA_real_)
  }
  TE <- rbindlist(te_all); CO <- rbindlist(co)
  fwrite(CO, file.path(CFG, sprintf("greg_mort_coefficients_%s.csv", drv)))
  score[[drv]] <- data.table(driver=drv, species=nrow(CO), n_test=nrow(TE),
     AUC=auc(TE$p, TE$died), logloss=logloss(TE$p, TE$died),
     median_abs_bdriver=if(drv=="none") NA_real_ else median(abs(CO$b_driver),na.rm=TRUE))
}
SC <- rbindlist(score)[order(-AUC)]
cat("\n=== mortality site-driver scorecard (held-out; higher AUC / lower logloss better) ===\n"); print(SC)
fwrite(SC, file.path(OUT, "mort_site_variant_scorecard.csv"))
cat("\ncoefficient sets: greg_mort_coefficients_{none,elev,bgi,cspi,esi}.csv in", CFG, "\n")
