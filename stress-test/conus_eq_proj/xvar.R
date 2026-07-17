suppressPackageStartupMessages(library(data.table))
V2<-'/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq'
V3<-'/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq_v3'
slope<-function(f){ if(!file.exists(f)) return(c(NA,NA)); m<-fread(f); if(!all(c('STAND_CN','PROJ_YEAR','TPH','QMD_IN','BA_FT2AC')%in%names(m))) return(c(NA,NA))
  m<-m[is.finite(TPH)&is.finite(QMD_IN)&TPH>0&QMD_IN>0]; setorder(m,STAND_CN,PROJ_YEAR)
  st<-m[,{lnN<-log(TPH);lnD<-log(QMD_IN); s<-if(.N>=4&&var(lnD)>0&&last(TPH)<first(TPH)*0.98) cov(lnN,lnD)/var(lnD) else NA_real_; .(s=s)},by=STAND_CN]
  ba<-m[PROJ_YEAR==max(PROJ_YEAR),mean(BA_FT2AC,na.rm=TRUE)]; c(median(st$s,na.rm=TRUE),ba) }
man<-fread(file.path(V2,'../eq_manifest.tsv'),header=FALSE)
out<-list()
for(i in 1:nrow(man)){ v<-tolower(man$V2[i]); md<-man$V3[i]; cfg<-if(md=='dependent')'conus_b2' else 'conus_b1'
  fn<-sprintf('conus_eq_%s_%s_metrics.csv',v,cfg); a<-slope(file.path(V2,fn)); b<-slope(file.path(V3,fn))
  out[[length(out)+1]]<-data.table(variant=toupper(v),mode=ifelse(md=='dependent','b2','b1'),
    v2_slope=round(a[1],2),v3_slope=round(b[1],2),v2_yr100BA=round(a[2],0),v3_yr100BA=round(b[2],0)) }
r<-rbindlist(out); fwrite(r,file.path(V3,'v3_vs_v2_full_constraint_check.csv'))
cat(sprintf('cells=%d | median slope v2=%.2f v3=%.2f | median yr100 BA v2=%.0f v3=%.0f\n',
  nrow(r),median(r$v2_slope,na.rm=T),median(r$v3_slope,na.rm=T),median(r$v2_yr100BA,na.rm=T),median(r$v3_yr100BA,na.rm=T)))
