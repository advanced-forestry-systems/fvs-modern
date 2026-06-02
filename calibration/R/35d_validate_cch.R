#!/usr/bin/env Rscript
# Validate the cch_organon port against the panel's stored CCH1.
# Self-contained: uses the panel's own per-tree DBH1/HT1/CR1 grouped by PLT_CN as
# the stand (approximation: only remeasured trees, but usually most of the
# canopy), computes cch via the ORGANON crown profile, and correlates with the
# stored CCH1. High correlation => the port + species crosswalk reproduce cch.
#
# Panel is metric (DBH1 cm, HT1 m); CAL_CCH.for is imperial -> convert.
# Species crosswalk (coarse): softwood (SPCD<300) -> ORGANON group 1 (DF),
# hardwood (>=300) -> group 16 (RA). Refine later if the signal is there.

suppressWarnings(suppressMessages(library(data.table)))
args <- commandArgs(trailingOnly=TRUE)
ga <- function(f,d=NULL){i<-grep(paste0("^--",f,"="),args,value=TRUE);if(!length(i))return(d);sub(paste0("^--",f,"="),"",i[1])}
DATA <- ga("data","data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT  <- ga("out","/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out")
NPLOT<- as.integer(ga("nplot","4000"))

# ---- ORGANON SWO crown params (groups used: 1=DF softwood, 16=RA hardwood) ----
MCW <- list("1"=c(4.6366,1.6078,-0.009625,88.52), "16"=c(8.0,1.53,0.0,999.99))
LCWp<- list("1"=c(0.0,0.00371834,0.808121),       "16"=c(0.3227140,0.0,0.0))
CWA <- list("1"=c(0.929973,-0.135212,-0.0157579), "16"=c(0.5,0.0,0.0))
DACB<- list("1"=0.062, "16"=0.0)
grp <- function(spcd) ifelse(spcd<300L, "1", "16")
mcw_f <- function(g,D,H){p<-MCW[[g]];d<-pmin(D,p[4]);ifelse(H<4.501,H/4.5*p[1],p[1]+p[2]*d+p[3]*d*d)}
hlcw_f<- function(g,H,CR) H-(1-DACB[[g]])*CR*H
lcw_f <- function(g,M,CR,D,H){p<-LCWp[[g]];CL<-CR*H;M*CR^(p[1]+p[2]*CL+p[3]*(D/H))}
cw_f  <- function(g,HL,LC,H,D,XL){p<-CWA[[g]];rp<-(H-XL)/(H-HL);rp[rp<=0]<-NA;LC*rp^(p[1]+p[2]*sqrt(rp)+p[3]*(H/D))}

plot_cch <- function(dt){          # dt: DBH(in), HT(ft), CR(0-1), EXPAN, g
  n<-nrow(dt); if(n<1) return(rep(NA,n))
  top<-max(dt$HT); cch<-numeric(41); cch[41]<-top
  for(i in 1:n){
    g<-dt$g[i];D<-dt$DBH[i];H<-dt$HT[i];CR<-dt$CR[i];E<-dt$EXPAN[i]
    CL<-CR*H;HCB<-H-CL;M<-mcw_f(g,D,H);LC<-lcw_f(g,M,CR,D,H);HL<-hlcw_f(g,H,CR)
    if(!is.finite(LC)||!is.finite(HL)||!is.finite(HCB)||!is.finite(H)) next
    thr<-max(HCB,HL)
    for(ii in 40:1){xl<-(ii-1)*(top/40); cw<-0
      if(xl<=thr){cw<-if(HCB<=HL) LC else { v<-cw_f(g,HL,LC,H,D,max(xl,HCB)); if(is.finite(v)) v else LC }}
      else if(xl<H){v<-cw_f(g,HL,LC,H,D,xl); if(is.finite(v)) cw<-v}
      if(!is.finite(cw)) cw<-0
      cch[ii]<-cch[ii]+(cw^2)*(0.001803*E)}
  }
  # interpolate to each tree's tip (defensive against non-finite / edge indices)
  sapply(dt$HT,function(h){
    if(!is.finite(h)||!is.finite(top)||top<=0||h>=top) return(0)
    xi<-40*(h/top); idx<-as.integer(xi)+1
    if(idx>=40) return(cch[40]*(40-xi))
    if(idx<2) idx<-2
    xxi<-(idx+1)-1; v<-cch[idx+1]+(cch[idx]-cch[idx+1])*(xxi-xi)
    if(is.finite(v)) v else 0})
}

cat("loading panel...\n"); d<-as.data.table(readRDS(DATA))
if("PLT_CN_cond1" %in% names(d)) setnames(d,"PLT_CN_cond1","PLT_CN")
keep<-c("PLT_CN","SPCD","DBH1","HT1","CR1","CCH1")
stopifnot(all(keep %in% names(d)))
d<-d[, ..keep]
d[, PLT_CN := as.character(PLT_CN)]   # avoid integer64 join-type issues[is.finite(DBH1)&is.finite(HT1)&is.finite(CR1)&is.finite(CCH1)&DBH1>0&HT1>0&CR1>0&CR1<=1]
# convert metric -> imperial; FIA expansion: 6.018 TPA (DIA>=5in), 74.965 (<5in)
d[, DBH:=DBH1/2.54][, HT:=HT1/0.3048][, CR:=CR1][, g:=grp(SPCD)]
d[, EXPAN:=fifelse(DBH>=5,6.018,74.965)]
set.seed(1); plots<-sample(unique(d$PLT_CN), min(NPLOT,uniqueN(d$PLT_CN)))
d<-d[PLT_CN %in% plots]
cat("plots:",length(plots)," trees:",nrow(d),"\n")
safe_cch <- function(dt) tryCatch(plot_cch(dt), error=function(e) rep(NA_real_, nrow(dt)))
d[, cch_hat := safe_cch(.SD), by=PLT_CN, .SDcols=c("DBH","HT","CR","EXPAN","g")]
cat("plots with NA cch_hat:", d[, sum(is.na(cch_hat))], "of", nrow(d), "\n")
r<-cor(d$cch_hat, d$CCH1, use="complete.obs")
rs<-cor(d$cch_hat, d$CCH1, use="complete.obs", method="spearman")
cat(sprintf("\n=== cch_hat vs stored CCH1 ===\nPearson r = %.3f   Spearman = %.3f\n", r, rs))
cat(sprintf("CCH1   range [%.3f, %.3f] mean %.3f\n", min(d$CCH1),max(d$CCH1),mean(d$CCH1)))
cat(sprintf("cch_hat range [%.3f, %.3f] mean %.3f\n", min(d$cch_hat,na.rm=TRUE),max(d$cch_hat,na.rm=TRUE),mean(d$cch_hat,na.rm=TRUE)))
fit<-lm(CCH1~cch_hat, d); cat(sprintf("CCH1 ~ a + b*cch_hat: a=%.3f b=%.4f R2=%.3f\n",
  coef(fit)[1],coef(fit)[2],summary(fit)$r.squared))
fwrite(d[, .(PLT_CN,SPCD,DBH,HT,CR,CCH1,cch_hat)], file.path(OUT,"cch_validation_sample.csv"))
cat("wrote cch_validation_sample.csv\n")
