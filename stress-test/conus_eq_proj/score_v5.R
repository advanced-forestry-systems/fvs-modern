#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE); MET<-a[1]; REGION<-toupper(a[2]); TAG<-a[3]
SIDIR<-'/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_seed_v4'
IN2CM<-2.54; AGE_LO<-10; AGE_HI<-40
SI<-fread(file.path(SIDIR,paste0('standinit_',REGION,'.csv')),colClasses=list(character='STAND_CN'),select=c('STAND_CN','AGE','SITE_INDEX'))
SI[,STAND_CN:=sub('\\..*$','',STAND_CN)]; SI[,AGE:=as.numeric(AGE)]; SI[,SITE_INDEX:=as.numeric(SITE_INDEX)]
SI<-unique(SI,by='STAND_CN'); SI[,seed:=is.finite(AGE)&AGE>=AGE_LO&AGE<=AGE_HI]
M<-fread(MET,colClasses=list(character='STAND_CN')); M[,STAND_CN:=sub('\\..*$','',STAND_CN)]
M[,QMD:=QMD_IN*IN2CM]; M<-M[is.finite(QMD)&is.finite(TPH)&TPH>0&QMD>0]
M<-merge(M,SI[,.(STAND_CN,AGE0=AGE,SI=SITE_INDEX,seed)],by='STAND_CN',all.x=TRUE)
M<-M[seed==TRUE & is.finite(SI)&SI>0]
M[,age:=AGE0+as.integer(PROJ_YEAR)]; M[,agebin:=round(age/5)*5]
# per-stand realized Reineke (median)
setorder(M,STAND_CN,PROJ_YEAR)
sl<-M[,{lnN<-log(TPH);lnD<-log(QMD); if(.N>=4&&var(lnD)>0&&last(TPH)<first(TPH)*0.98) as.numeric(cov(lnN,lnD)/var(lnD)) else NA_real_},by=STAND_CN]$V1
# within-region quartiles -> per-site aggregated trajectory + per-site Reineke slope on aggregated means
qb<-quantile(M$SI,c(0,.25,.5,.75,1),na.rm=TRUE)
M[,site:=cut(SI,breaks=qb,include.lowest=TRUE,labels=1:4)]
agg<-M[,.(HT_DOM=mean(HT_M_DOM,na.rm=TRUE),HT_MEAN=mean(HT_M_MEAN,na.rm=TRUE),QMD=mean(QMD),TPH=mean(TPH),BA=mean(BA_FT2AC),n=.N),by=.(site,agebin)]
agg<-agg[n>=5][order(site,agebin)]
res<-agg[,{ok<-is.finite(QMD)&is.finite(TPH)&QMD>0&TPH>0
  s<-if(sum(ok)>=3) as.numeric(cov(log(TPH[ok]),log(QMD[ok]))/var(log(QMD[ok]))) else NA_real_
  htmono<-all(diff(HT_DOM)>=-1e-6,na.rm=TRUE); qmono<-all(diff(QMD)>=-1e-6,na.rm=TRUE)
  htlast<-HT_DOM[which.max(agebin)]; htmax<-max(HT_DOM,na.rm=TRUE)
  .(reineke=round(s,2),ht_dom_mono=htmono,qmd_mono=qmono,ht_drop=round(htmax-htlast,2))},by=site]
cat(sprintf('=== TAG=%s REGION=%s per-stand median Reineke=%.2f ===\n',TAG,REGION,median(sl,na.rm=TRUE)))
print(res)
