#!/usr/bin/env Rscript
## validate_ne_greg.R -- NE GO-bar validation for the faithful Greg arm.
## GO criteria: (1) year-0 BA/QMD/TPH match other arms (same seed); (2) DG+mort spot-check
## reproduces Greg's own est_dg/survival; (3) EMT/TD present for NE; (4) trajectories bounded;
## (5) fallback fraction reported. Writes validation_summary_NE_greg.csv + thumbnail.
suppressPackageStartupMessages({library(data.table)})
P<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
OUTD<-file.path(P,"out_conus_eq_greg")
GREG<-file.path(OUTD,"conus_eq_ne_greg_metrics.csv")
FB<-file.path(OUTD,"conus_eq_ne_greg_fallback.csv")
RDS<-"/users/PUOM0008/crsfaaron/fvs_remodeling/rds"
stopifnot(file.exists(GREG))
g<-fread(GREG); fb<-if(file.exists(FB)) fread(FB) else NULL

## ---- (1) year-0 seed match vs the v3 fvs-conus arm (same standinit+treeinit) ----
seed_ref<-NULL
for(cand in c(file.path(P,"out_conus_eq_v3","conus_eq_ne_conus_b2_metrics.csv"),
              file.path(P,"conus_eq_ne_conus_b2_metrics.csv"))){
  if(file.exists(cand)){ seed_ref<-fread(cand); break } }
y0g<-g[PROJ_YEAR==0,.(BA0=mean(BA_FT2AC,na.rm=TRUE),QMD0=mean(QMD_IN,na.rm=TRUE),TPH0=mean(TPH,na.rm=TRUE),n=.N)]
cat("=== (1) YEAR-0 SEED MATCH (identical stands) ===\n")
cat(sprintf("  GREG arm year-0: BA=%.2f QMD=%.2f TPH=%.1f n=%d\n",y0g$BA0,y0g$QMD0,y0g$TPH0,y0g$n))
if(!is.null(seed_ref)){
  ## match on common STAND_CN at year 0
  rr<-seed_ref[PROJ_YEAR==0]
  common<-intersect(g[PROJ_YEAR==0]$STAND_CN, rr$STAND_CN)
  gg<-g[PROJ_YEAR==0][STAND_CN %in% common]; rrc<-rr[STAND_CN %in% common]
  setkey(gg,STAND_CN); setkey(rrc,STAND_CN); rrc<-rrc[gg$STAND_CN]
  cat(sprintf("  REF (conus_b2) year-0 on %d common stands: BA=%.2f QMD=%.2f TPH=%.1f\n",
    length(common),mean(rrc$BA_FT2AC,na.rm=TRUE),mean(rrc$QMD_IN,na.rm=TRUE),mean(rrc$TPH,na.rm=TRUE)))
  cat(sprintf("  GREG  year-0 on same %d stands:        BA=%.2f QMD=%.2f TPH=%.1f\n",
    length(common),mean(gg$BA_FT2AC,na.rm=TRUE),mean(gg$QMD_IN,na.rm=TRUE),mean(gg$TPH,na.rm=TRUE)))
  ba_d<-abs(mean(gg$BA_FT2AC,na.rm=TRUE)-mean(rrc$BA_FT2AC,na.rm=TRUE))/mean(rrc$BA_FT2AC,na.rm=TRUE)*100
  tph_d<-abs(mean(gg$TPH,na.rm=TRUE)-mean(rrc$TPH,na.rm=TRUE))/mean(rrc$TPH,na.rm=TRUE)*100
  cat(sprintf("  year-0 |dBA|=%.2f%%  |dTPH|=%.2f%%  -> %s\n",ba_d,tph_d,if(ba_d<2 && tph_d<2)"PASS (identical seed)" else "CHECK"))
} else cat("  (no v3 reference metrics found for seed comparison)\n")

## ---- (2) DG + mortality spot-check vs Greg's OWN equations ----
cat("\n=== (2) SPECIES SPOT-CHECK vs Greg's own est_dg / survival (his params) ===\n")
dg<-as.data.table(readRDS(file.path(RDS,"dg_parms.RDS")))
mo<-as.data.table(readRDS(file.path(RDS,"mort_parm_base_rate_cr_cch.RDS"))); mo<-mo[order(SPCD,nll)][!duplicated(SPCD)]
est_dg<-function(n,dbh0,bal,cr,ht,elev,emt,B0,B1,B2,B3,B4,B5,B6){cdbh<-dbh0
  for(i in 1:n){dgh<-exp(B0+B1*log((cdbh+1)^2/(cr*ht+1)^B3)+B2*bal^B4/log(cdbh+2.7)+B5*elev+B6*emt); cdbh<-cdbh+dgh}; cdbh}
surv1<-function(b,cr,cch)1-exp(-exp(b[1]+b[2]*(cr+0.01)^b[3]+b[4]*cch^b[5]))
spp<-c(316,12,97)  # red maple, balsam fir, red spruce (NE-dominant)
nm<-c("316"="red maple","12"="balsam fir","97"="red spruce")
spot<-rbindlist(lapply(spp,function(s){
  pd<-dg[spcd==s]; pm<-mo[SPCD==s]
  dgyr<-if(nrow(pd)){d5<-est_dg(5,10,80,0.5,50,940,-29,pd$B0,pd$B1,pd$B2,pd$B3,pd$B4,pd$B5,pd$B6);(d5-10)/5} else NA
  mann<-if(nrow(pm)) 1-surv1(c(pm$b0,pm$b1,pm$b2,pm$b3,pm$b4),0.3,0.6) else NA
  data.table(SPCD=s,Species=nm[as.character(s)],DG_in_yr=round(dgyr,4),
             has_DG=nrow(pd)>0,annMort_cr.3_cch.6=round(mann,4),has_MORT=nrow(pm)>0)}))
print(spot)
cat("  (DG at dbh0=10in bal=80 cr=.5 ht=50ft elev=940 emt=-29; mort at cr=0.3 cch=0.6 -- Greg eqs verbatim)\n")

## ---- (3) EMT/TD present for NE ----
cat("\n=== (3) EMT/TD COVERAGE (NE) ===\n")
look<-as.data.table(readRDS(file.path(P,"greg_emt_td_lookup.rds")))
si<-fread("/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant/standinit_NE.csv",colClasses=list(character="STAND_CN"),select="STAND_CN")
si[,STAND_CN:=sub("\\..*$","",STAND_CN)]; u<-unique(si$STAND_CN); m<-look[match(u,STAND_CN)]
cat(sprintf("  NE stands=%d  EMT/TD coverage=%.2f%%  EMT range %.1f..%.1f  TD range %.1f..%.1f (raster: ClimateNA 1991-2020 EMT.tif; TD=MWMT-MCMT)\n",
  length(u),100*mean(is.finite(m$EMT)),min(m$EMT,na.rm=TRUE),max(m$EMT,na.rm=TRUE),min(m$TD,na.rm=TRUE),max(m$TD,na.rm=TRUE)))

## ---- (4) trajectory boundedness ----
cat("\n=== (4) TRAJECTORY (stand-mean) ===\n")
traj<-g[,.(BA=round(mean(BA_FT2AC,na.rm=TRUE),1),QMD=round(mean(QMD_IN,na.rm=TRUE),2),TPH=round(mean(TPH,na.rm=TRUE),0)),by=PROJ_YEAR][order(PROJ_YEAR)]
print(traj[PROJ_YEAR %in% c(0,20,40,60,80,100)])
peak<-traj[which.max(BA)]; bounded<-all(is.finite(traj$BA)) && max(traj$BA)<400 && min(traj$TPH)>0
cat(sprintf("  peak BA=%.1f at yr %d; bounded & finite=%s\n",peak$BA,peak$PROJ_YEAR,bounded))

## ---- (5) fallback fraction ----
cat("\n=== (5) PER-STAND FALLBACK FRACTION (tree-count weighted) ===\n")
if(!is.null(fb)){
  cat(sprintf("  DG fallback: mean=%.1f%% (median %.1f%%)  HG: mean=%.1f%%  MORT: mean=%.1f%%  n=%d stands\n",
    100*mean(fb$fb_dg),100*median(fb$fb_dg),100*mean(fb$fb_hg),100*mean(fb$fb_mo),nrow(fb)))
  cat(sprintf("  Greg coverage: DG=%.1f%%  HG=%.1f%%  MORT=%.1f%%\n",100*(1-mean(fb$fb_dg)),100*(1-mean(fb$fb_hg)),100*(1-mean(fb$fb_mo))))
}

## ---- write summary ----
summ<-data.table(
  arm="greg",variant="NE",
  n_stands=y0g$n,
  BA0=round(y0g$BA0,2),QMD0=round(y0g$QMD0,2),TPH0=round(y0g$TPH0,1),
  BA50=traj[PROJ_YEAR==50]$BA,QMD50=traj[PROJ_YEAR==50]$QMD,TPH50=traj[PROJ_YEAR==50]$TPH,
  BA100=traj[PROJ_YEAR==100]$BA,QMD100=traj[PROJ_YEAR==100]$QMD,TPH100=traj[PROJ_YEAR==100]$TPH,
  peakBA=peak$BA,peakBA_yr=peak$PROJ_YEAR,
  EMT_cov_pct=round(100*mean(is.finite(m$EMT)),2),
  fb_DG_pct=if(!is.null(fb))round(100*mean(fb$fb_dg),2) else NA,
  fb_HG_pct=if(!is.null(fb))round(100*mean(fb$fb_hg),2) else NA,
  fb_MORT_pct=if(!is.null(fb))round(100*mean(fb$fb_mo),2) else NA,
  bounded=bounded)
fwrite(summ,file.path(P,"validation_summary_NE_greg.csv"))
fwrite(spot,file.path(P,"validation_spotcheck_NE_greg.csv"))
cat("\nWrote validation_summary_NE_greg.csv and validation_spotcheck_NE_greg.csv\n")

## ---- thumbnail ----
png(file.path(P,"thumb_BA_QMD_TPH_NE_greg.png"),width=1100,height=360,res=110)
par(mfrow=c(1,3),mar=c(4,4,2,1))
plot(traj$PROJ_YEAR,traj$BA,type="o",pch=19,col="#2c7fb8",xlab="Projection year",ylab="BA (ft2/ac)",main="NE GREG: BA")
plot(traj$PROJ_YEAR,traj$QMD,type="o",pch=19,col="#31a354",xlab="Projection year",ylab="QMD (in)",main="NE GREG: QMD")
plot(traj$PROJ_YEAR,traj$TPH,type="o",pch=19,col="#de2d26",xlab="Projection year",ylab="TPH",main="NE GREG: TPH")
dev.off()
cat("Wrote thumb_BA_QMD_TPH_NE_greg.png\n")
cat("\nDONE.\n")
