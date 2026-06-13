#!/usr/bin/env Rscript
# stand_state_v2.R -- Garcia-grounded refits: basal area carrying capacity from SDImax,
# density self-thinning as a power-form rate, coupling BA growth to height increment.
suppressPackageStartupMessages({ library(data.table) })
set.seed(1)
d <- as.data.table(readRDS("data/conus_remeasurement_pairs_metric_cond_v2.rds")); nm<-names(d)
if(!("TPH1" %in% nm) && "TPA1" %in% nm){ d[,TPH1:=TPA1*2.47105]; d[,TPH2:=TPA2*2.47105] }
sdimax_col <- intersect(c("SDImax_brms","SDImax","sdimax"), nm)[1]
d <- d[is.finite(YEARS)&YEARS>=1&YEARS<=20]
topht <- function(ht,dbh){ thr<-quantile(dbh,0.8,na.rm=TRUE); mean(ht[dbh>=thr],na.rm=TRUE) }
pl <- d[, .(H1=topht(HT1,DBH1),H2=topht(HT2,DBH2),dt=YEARS[1],N1=TPH1[1],N2=TPH2[1],
            BA1=BA1[1],BA2=BA2[1],SDI1=SDI1[1],SDImax=if(length(sdimax_col))get(sdimax_col)[1] else NA,
            QMD1=QMD1[1],QMD2=QMD2[1],cspi=cspi[1],ntree=.N), by=plot_key]
pl <- pl[ntree>=5 & is.finite(SDImax)&SDImax>0]
pl[, SDImax_z:=as.numeric(scale(SDImax))]
r2<-function(o,p)1-sum((o-p)^2)/sum((o-mean(o))^2); rmse<-function(o,p)sqrt(mean((o-p)^2))
cat("plot pairs:", nrow(pl), "\n\n")

## (C') BASAL AREA: Gmax driven by SDImax (carrying capacity), rate coupled to height increment
tryCatch({
  s<-pl[is.finite(BA1)&is.finite(BA2)&BA1>0&BA2>0&is.finite(H1)&is.finite(H2)]
  s[, dH:=pmax(H2-H1,0)]
  BA1<-s$BA1;BA2<-s$BA2;DT<-s$dt;SM<-s$SDImax;SMz<-s$SDImax_z;dH<-s$dH
  # Gmax = g0 + g1*SDImax ; BA2 = BA1 + (Gmax-BA1)*(1-exp(-k*dt))
  sse<-function(p){Gmax<-p[1]+p[2]*SM; k<-p[3]; pr<-BA1+(Gmax-BA1)*(1-exp(-k*DT)); if(any(!is.finite(pr)))return(1e12); sum((BA2-pr)^2)}
  f<-nlminb(c(10,0.05,0.05),sse,lower=c(0,0,0.001),upper=c(80,0.6,0.6))
  Gmax<-f$par[1]+f$par[2]*SM; pr<-BA1+(Gmax-BA1)*(1-exp(-f$par[3]*DT))
  cat(sprintf("(C') BA Gmax=%.1f + %.3f*SDImax  k=%.4f/yr | g1 sign %s | RMSE=%.3f R2_inc=%.3f\n",
      f$par[1],f$par[2],f$par[3], ifelse(f$par[2]>0,"POSITIVE (expected)","NEG flag"),
      rmse(BA2,pr), r2(BA2-BA1, pr-BA1)))
}, error=function(e) cat("(C') BA FAILED:",conditionMessage(e),"\n"))

## (B') DENSITY self-thinning, power form: lnN2 = lnN1 - exp(c0)*RD^c1 * dt
tryCatch({
  s<-pl[is.finite(N1)&is.finite(N2)&N1>0&N2>0&is.finite(SDI1)]
  s[, RD:=SDI1/SDImax]; s<-s[RD>0.02&RD<1.5]
  lnN1<-log(s$N1);lnN2<-log(s$N2);RD<-s$RD;DT<-s$dt
  sse<-function(p){ pr<-lnN1 - exp(p[1])*RD^p[2]*DT; sum((lnN2-pr)^2) }
  f<-nlminb(c(log(0.02),1.5),sse,lower=c(log(1e-4),0.3),upper=c(log(2),5))
  pr<-lnN1-exp(f$par[1])*RD^f$par[2]*DT
  cat(sprintf("(B') density: rate = %.4f * RD^%.2f /yr | RMSE(lnN)=%.3f R2=%.3f\n",
      exp(f$par[1]),f$par[2],rmse(lnN2,pr),r2(lnN2,pr)))
  for(rd in c(0.3,0.6,1.0)) cat(sprintf("    RD=%.1f -> %.2f%%/yr loss\n", rd, 100*exp(f$par[1])*rd^f$par[2]))
}, error=function(e) cat("(B') density FAILED:",conditionMessage(e),"\n"))
cat("\nDONE_STAND_V2\n")
