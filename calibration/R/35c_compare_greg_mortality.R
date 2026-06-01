#!/usr/bin/env Rscript
# Model-level comparison of Greg's re-fit mortality (gompit on cr + cch) against
# the per-species base rate, on the held remeasurement panel. Quantifies what the
# new covariates buy: discrimination (AUC), calibration (predicted vs observed
# survival), and whether the model tracks the crown-closure-at-tip (cch) signal
# Greg's paper highlights. Pure prediction check -- no FVS engine needed.
#
# Usage:
#   Rscript compare_new_mortality.R --coeffs full_out/greg_mortality_coefficients.csv \
#       --panel mort_slim.rds --out full_out
suppressWarnings(suppressMessages(library(data.table)))
args <- commandArgs(trailingOnly = TRUE)
ga <- function(f,d=NULL){i<-grep(paste0("^--",f,"="),args,value=TRUE);if(!length(i))return(d);sub(paste0("^--",f,"="),"",i[1])}
COEF <- ga("coeffs"); PANEL <- ga("panel"); OUT <- ga("out","."); stopifnot(!is.null(COEF),!is.null(PANEL))

co <- fread(COEF); setkey(co, SPCD)
d  <- as.data.table(readRDS(PANEL))[SPCD %in% co$SPCD]
d  <- merge(d, co[,.(SPCD,b0,b1,b2,b3,b4)], by="SPCD")

# new-model annual hazard + period survival
d[, eta := b0 + b1*(cr+0.01)^b2 + b3*fifelse(cch>0, cch^b4, 0)]
d[, eta := pmin(pmax(eta,-30),30)]
d[, psurv_new := exp(-exp(eta)*years)]
# per-species base rate: annual survival ^ years
base <- d[, .(p_ann = mean(alive)^(1/mean(years))), by=SPCD]   # crude annualization
d <- merge(d, base, by="SPCD")
d[, psurv_base := p_ann^years]

# metrics (died = event for AUC so higher prob = more death)
d[, pdie_new := 1-psurv_new]; d[, pdie_base := 1-psurv_base]; d[, died := 1L-alive]
fast_auc <- function(score, y){           # Mann-Whitney AUC, sampled for speed
  i1<-which(y==1); i0<-which(y==0); if(!length(i1)||!length(i0)) return(NA_real_)
  n<-min(40000,length(i1),length(i0)); s1<-sample(score[i1],n); s0<-sample(score[i0],n)
  mean(outer(s1[1:min(2000,n)], s0[1:min(2000,n)], ">")) }
ll <- function(p,y){p<-pmin(pmax(p,1e-9),1-1e-9); -mean(y*log(p)+(1-y)*log(1-p))}

overall <- data.table(
  n=nrow(d), obs_surv=round(mean(d$alive),4),
  pred_surv_new=round(mean(d$psurv_new),4), pred_surv_base=round(mean(d$psurv_base),4),
  auc_new=round(fast_auc(d$pdie_new,d$died),4), auc_base=round(fast_auc(d$pdie_base,d$died),4),
  logloss_new=round(ll(d$pdie_new,d$died),4), logloss_base=round(ll(d$pdie_base,d$died),4))
fwrite(overall, file.path(OUT,"mortality_compare_overall.csv"))
print(overall)

# survival vs cch (Greg's key signal): observed vs new-predicted by cch quartile
d[, cch_q := cut(cch, quantile(cch, 0:4/4, na.rm=TRUE), include.lowest=TRUE)]
by_cch <- d[, .(n=.N, obs_surv=round(mean(alive),4),
                pred_new=round(mean(psurv_new),4), pred_base=round(mean(psurv_base),4)),
            by=cch_q][order(cch_q)]
fwrite(by_cch, file.path(OUT,"mortality_compare_by_cch.csv"))
cat("\n=== survival by cch quartile ===\n"); print(by_cch)
cat("\nwrote mortality_compare_overall.csv + mortality_compare_by_cch.csv\n")
