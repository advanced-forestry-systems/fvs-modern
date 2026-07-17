#!/usr/bin/env Rscript
# build_singlecohort_traj.R -- SINGLE-COHORT (even-aged) site-class trajectories.
# Follows the SAME stands forward by PROJ_YEAR (so trajectories do not loop from cohort mixing).
# Selects a young, even-aged seed set: starting AGE in [AGE_LO,AGE_HI], then groups by
# SITE_INDEX quartile and reports mean metrics vs PROJ_YEAR. REAL HT (HT_M_DOM) used when present.
# Usage: Rscript build_singlecohort_traj.R <metrics.csv> <variant> <label> <outdir> [age_lo age_hi]
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE)
MET<-a[1]; VAR<-a[2]; LAB<-a[3]; OUT<-a[4]
AGE_LO<-if(length(a)>=5) as.numeric(a[5]) else 0
AGE_HI<-if(length(a)>=6) as.numeric(a[6]) else 40
SIDIR<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
FT2AC_to_M2HA<-0.2296; IN_to_CM<-2.54; FF<-0.40
dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
si<-fread(file.path(SIDIR,paste0("standinit_",VAR,".csv")),colClasses=list(character="STAND_CN"),select=c("STAND_CN","AGE","SITE_INDEX"))
si[,STAND_CN:=sub("\\\\..*$","",STAND_CN)]
si[,AGE:=suppressWarnings(as.numeric(AGE))]; si[,SITE_INDEX:=suppressWarnings(as.numeric(SITE_INDEX))]
si<-unique(si,by="STAND_CN")
qb<-quantile(si[SITE_INDEX>0]$SITE_INDEX,c(0,.25,.5,.75,1),na.rm=TRUE)
si[,siteclass:=cut(SITE_INDEX,qb,include.lowest=TRUE,labels=1:4)]
# single-cohort seed: even-aged young stands
seed<-si[is.finite(AGE)&AGE>=AGE_LO&AGE<=AGE_HI&!is.na(siteclass)]
cat(sprintf("single-cohort seed: %d stands with start AGE in [%g,%g]\n",nrow(seed),AGE_LO,AGE_HI))
m<-fread(MET,colClasses=list(character="STAND_CN"))
m[,STAND_CN:=sub("\\\\..*$","",STAND_CN)]
m<-merge(m,seed[,.(STAND_CN,siteclass)],by="STAND_CN")
m[,BAPH:=BA_FT2AC*FT2AC_to_M2HA]; m[,QMD:=QMD_IN*IN_to_CM]; m[,TPH:=as.numeric(TPH)]
# REAL height: prefer HT_M_DOM (top height) -> Eichhorn; fall back HT_M_MEAN
if("HT_M_DOM" %in% names(m)){ m[,HT:=as.numeric(HT_M_DOM)] } else m[,HT:=NA_real_]
if("HT_M_MEAN" %in% names(m)){ miss<-!is.finite(m$HT); m[miss,HT:=as.numeric(HT_M_MEAN)] }
m<-m[is.finite(BAPH)&is.finite(QMD)&is.finite(TPH)&TPH>0&QMD>0]
m[,age:=as.numeric(PROJ_YEAR)]    # PROJ_YEAR = cohort age since seeding (even-aged origin)
m<-m[is.finite(HT)]
m[,VOL:=BAPH*HT*FF]; m[,agebin:=round(age/5)*5]
agg<-m[,.(HT=mean(HT),TPH=mean(TPH),BAPH=mean(BAPH),QMD=mean(QMD),VOL=mean(VOL),n=.N),by=.(siteclass,agebin)]
setnames(agg,c("siteclass","agebin"),c("site","age")); agg<-agg[n>=5][order(site,age)]; agg[,origin:="all"]
of<-file.path(OUT,paste0(LAB,"_sc_traj.csv"))
fwrite(agg[,.(origin,site,age,HT,TPH,BAPH,QMD,VOL,n)],of)
cat(sprintf("wrote %s  (rows=%d; HT=%s)\n",of,nrow(agg), if("HT_M_DOM"%in%names(m))"real HT_M_DOM" else "none"))
