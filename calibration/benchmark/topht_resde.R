#!/usr/bin/env Rscript
# topht_resde.R -- Garcia reducible-SDE top-height fit WITH measurement error (his method).
# Bertalanffy-Richards SDE; compare measurement-error ON (mum=1) vs OFF (mum=0) to show how
# much apparent top-height variation is observation noise (the cause of negative increment R2).
suppressPackageStartupMessages({ library(resde); library(data.table) })
set.seed(1)
d <- as.data.table(readRDS("data/conus_remeasurement_pairs_metric_cond_v2.rds"))
d <- d[is.finite(YEARS)&YEARS>=1&YEARS<=20&is.finite(HT1)&is.finite(HT2)&is.finite(DBH1)]
topht <- function(ht,dbh){ thr<-quantile(dbh,0.8,na.rm=TRUE); mean(ht[dbh>=thr],na.rm=TRUE) }
pl <- d[, .(H1=topht(HT1,DBH1), H2=topht(HT2,DBH2), dt=YEARS[1], ntree=.N), by=plot_key]
pl <- pl[ntree>=5 & is.finite(H1)&is.finite(H2)&H1>1.5&H2>1.5&H1<90&H2<90]
pl[, hinc:=(H2-H1)/dt]; pl <- pl[hinc> -0.25 & hinc<2.0]
if(nrow(pl)>15000) pl <- pl[sample(.N,15000)]
long <- rbind(pl[,.(unit=plot_key, t=0, x=H1)], pl[,.(unit=plot_key, t=dt, x=H2)])[order(unit,t)]
cat("plots:", nrow(pl), " obs:", nrow(long), "\n\n")

rich <- function(mum) sdemodel(phi=~x^c, beta0=~b*A^c, beta1=~-b, mum=mum)
fitone <- function(mum,label){ tryCatch({
  f <- sdefit(rich(mum), x="x", t="t", unit="unit", data=long,
              global=c(A=38, b=0.04, c=-0.6), method="nls")
  p <- f$fit; m <- f$more
  cat(sprintf("[%s] params: %s | sigma_proc=%.4f sigma_meas=%.4f | logLik=%.1f AIC=%.1f\n",
      label, paste(names(coef(p)),round(coef(p),4),sep="=",collapse=" "),
      ifelse(is.null(m$sigma),NA,m$sigma[1]), ifelse(length(m$sigma)>1,m$sigma[2],0),
      as.numeric(m$logLik), as.numeric(m$AIC)))
}, error=function(e) cat(sprintf("[%s] FAILED: %s\n",label,conditionMessage(e)))) }

cat("==== Bertalanffy-Richards top-height SDE ====\n")
fitone(0, "no measurement error (mum=0)")
fitone(1, "WITH measurement error (mum=1)")
cat("\nKey: if sigma_meas is a large share of total, the negative increment skill of the\n")
cat("least-squares GADA was measurement noise, and the SDE recovers the real trajectory.\n")
cat("DONE_TOPHT_RESDE\n")
