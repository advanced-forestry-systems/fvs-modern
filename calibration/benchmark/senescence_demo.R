#!/usr/bin/env Rscript
# senescence_demo.R -- does adding a relative-size term recover the observed U-shaped mortality?
# cloglog mortality with exposure: base (absolute size) vs +relsize (DBH/species max).
suppressPackageStartupMessages({ library(data.table) })
set.seed(1)
d <- as.data.table(readRDS("data/conus_remeasurement_pairs_metric_cond_v2.rds")); nm<-names(d)
tr <- as.data.table(readRDS("traits/species_traits_v2.rds"))
md <- tr[, .(SPCD, max_dbh=max_dbh_cm)]
st2<-intersect(c("TREESTATUS2","STATUS2"),nm)[1]; st1<-intersect(c("TREESTATUS1","STATUS1"),nm)[1]
d <- d[is.finite(YEARS)&YEARS>=1&YEARS<=20&is.finite(DBH1)&DBH1>=2.54&get(st1)==1&get(st2)%in%c(1,2)&is.finite(BA1)&is.finite(BAL_SW1)&is.finite(BAL_HW1)]
d[, died:=as.integer(get(st2)==2)]
d <- merge(d, md, by="SPCD", all.x=TRUE); d <- d[is.finite(max_dbh)&max_dbh>0]
d[, relsz:=DBH1/max_dbh]; d <- d[relsz>0&relsz<1.3]
if(nrow(d)>400000) d <- d[sample(.N,400000)]
d[, ln_dbh:=log(DBH1)]; d[, bal:=BAL_SW1+BAL_HW1]; d[, ba:=BA1*0.2296]; d[, logT:=log(YEARS)]
cat("n:", nrow(d), " overall dead frac:", round(mean(d$died),4), "\n\n")
# BASE: absolute size only (mirrors current model: ln_dbh + ln_dbh^2)
mB <- glm(died ~ ln_dbh + I(ln_dbh^2) + ba + bal, family=binomial("cloglog"), offset=logT, data=d)
# +RELSIZE senescence: add relsz + relsz^2
mR <- glm(died ~ ln_dbh + I(ln_dbh^2) + ba + bal + relsz + I(relsz^2), family=binomial("cloglog"), offset=logT, data=d)
cat(sprintf("relsz coef=%.3f (p=%.1e)  relsz^2 coef=%.3f (p=%.1e)  dAIC(base-rel)=%.0f\n",
    coef(mR)["relsz"], summary(mR)$coef["relsz",4], coef(mR)["I(relsz^2)"], summary(mR)$coef["I(relsz^2)",4], AIC(mB)-AIC(mR)))
# predicted ANNUAL mortality (set offset to 0 => 1yr) by relsize class, both models vs observed
d[, eta_B:=predict(mB, newdata=d) - logT]; d[, eta_R:=predict(mR, newdata=d) - logT]
d[, pB:=1-exp(-exp(eta_B))]; d[, pR:=1-exp(-exp(eta_R))]
d[, rel_cls:=cut(relsz, c(0,0.15,0.3,0.45,0.6,0.75,0.9,1.3), right=FALSE)]
tab <- d[, .(n=.N, obs=round(1-(mean(1-died))^(1/mean(YEARS)),4), base_pred=round(mean(pB),4), relsize_pred=round(mean(pR),4)), by=rel_cls][order(rel_cls)]
cat("\n== annual mortality by relative size: observed vs base vs +relsize ==\n"); print(tab)
cat("\nFix works if relsize_pred rises with the observed at high relative size and base_pred does not.\n")
cat("DONE_SENESCENCE\n")
