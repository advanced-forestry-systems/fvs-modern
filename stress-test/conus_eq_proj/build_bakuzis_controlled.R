#!/usr/bin/env Rscript
# Controlled Bakuzis trajectories for the headline (constrained vs unconstrained):
#  (A) MATCHED: conus_b2/b1 vs sdicon_b2/b1 on the SAME 16 variants, SAME HT method
#      (derived HT=f(QMD)) so Reineke / Eichhorn / site-ordering differences isolate the
#      stand-level SDI constraint, not the HT source or the variant mix.
#  (B) WITHIN-REGION site ordering: PN variant only, origin=PN, so site classes are a single
#      species pool (SITE_INDEX comparable) -> site ordering is interpretable.
suppressPackageStartupMessages(library(data.table))
SCR <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
EQ  <- file.path(SCR,"conus_eq_proj")
OUT <- file.path(EQ,"bakuzis_redo"); dir.create(OUT, showWarnings=FALSE)
SIDIR <- file.path(SCR,"standinit_by_variant")
FF<-0.40; FT2AC<-0.2296; IN2CM<-2.54
SDI_VARS <- toupper(c("ak","bm","ca","ci","cr","cs","ec","em","ie","nc","pn","so","tt","ut","wc","ws"))

ff_ <- readRDS(file.path(OUT,"htqmd_fit.rds")); a<-ff_$a; b<-ff_$b
ht_qmd <- function(QMD) a*QMD^b
cat(sprintf("HT=f(QMD): %.3f*QMD^%.3f (cached)\n", a, b))

# standinit AGE + SITE_INDEX
SI <- rbindlist(lapply(list.files(SIDIR,pattern="^standinit_[A-Z]+\\.csv$",full.names=TRUE),
  function(f) tryCatch(fread(f,colClasses=list(character="STAND_CN"),
    select=c("STAND_CN","AGE","SITE_INDEX")),error=function(e) NULL)),fill=TRUE)
SI<-SI[!is.na(STAND_CN)]; SI[,AGE:=as.numeric(AGE)]; SI[,SITE_INDEX:=as.numeric(SITE_INDEX)]
SI<-unique(SI,by="STAND_CN")
setkey(SI,STAND_CN)

# global quartile boundaries (for matched, pooled across the 16 vars) - reuse global SI dist
qb_g <- quantile(SI$SITE_INDEX[is.finite(SI$SITE_INDEX)&SI$SITE_INDEX>0],c(0,.25,.5,.75,1),na.rm=TRUE)

build <- function(metrics_files, label, config_filter=NULL, variant_filter=NULL,
                  origin_by_variant=FALSE, qb=qb_g){
  M<-rbindlist(lapply(metrics_files,function(f) tryCatch(fread(f,colClasses=list(character="STAND_CN")),
        error=function(e) NULL)),fill=TRUE)
  if(!nrow(M)){cat(" NO ",label,"\n");return(invisible())}
  if(!is.null(config_filter)) M<-M[CONFIG %in% config_filter]
  if(!is.null(variant_filter)) M<-M[toupper(VARIANT) %in% variant_filter]
  M[,BAPH:=BA_FT2AC*FT2AC][,QMD:=QMD_IN*IN2CM][,TPH:=as.numeric(TPH)]
  M<-M[is.finite(BAPH)&is.finite(QMD)&is.finite(TPH)&TPH>0&QMD>0]
  M<-merge(M,SI[,.(STAND_CN,AGE0=AGE,SI=SITE_INDEX)],by="STAND_CN",all.x=TRUE)
  M[,age:=AGE0+as.numeric(PROJ_YEAR)]; M<-M[is.finite(age)]
  # site class: within-variant quartiles if origin_by_variant, else global
  if(origin_by_variant){
    M<-M[is.finite(SI)&SI>0]
    M[,siteclass:=cut(SI,breaks=quantile(SI,c(0,.25,.5,.75,1),na.rm=TRUE),include.lowest=TRUE,labels=1:4),by=VARIANT]
    M[,origin:=toupper(VARIANT)]
  } else {
    M[,siteclass:=cut(SI,breaks=qb,include.lowest=TRUE,labels=1:4)]
    M[,origin:="all"]
  }
  M<-M[!is.na(siteclass)]
  M[,HT:=ht_qmd(QMD)]; M[,VOL:=BAPH*HT*FF]; M[,agebin:=round(age/5)*5]
  agg<-M[,.(HT=mean(HT),TPH=mean(TPH),BAPH=mean(BAPH),QMD=mean(QMD),VOL=mean(VOL),n=.N),
        by=.(origin,siteclass,agebin)]
  setnames(agg,c("siteclass","agebin"),c("site","age")); agg<-agg[n>=5][order(origin,site,age)]
  of<-file.path(OUT,paste0(label,"_traj.csv")); fwrite(agg[,.(origin,site,age,HT,TPH,BAPH,QMD,VOL,n)],of)
  cat(sprintf("  %-24s rows=%4d origins=%d -> %s\n",label,nrow(agg),uniqueN(agg$origin),basename(of)))
}

eqd<-file.path(EQ,"out_conus_eq"); sdd<-file.path(EQ,"out_conus_eq_sdiconstrained")
eq_b2<-list.files(eqd,pattern="_conus_b2_metrics\\.csv$",full.names=TRUE)
eq_b1<-list.files(eqd,pattern="_conus_b1_metrics\\.csv$",full.names=TRUE)
sd_b2<-list.files(sdd,pattern="_conus_b2_metrics\\.csv$",full.names=TRUE)
sd_b1<-list.files(sdd,pattern="_conus_b1_metrics\\.csv$",full.names=TRUE)

cat("\n(A) MATCHED 16-variant, derived-HT comparison:\n")
build(eq_b2,"conus_b2_m16", variant_filter=SDI_VARS)
build(eq_b1,"conus_b1_m16", variant_filter=SDI_VARS)
build(sd_b2,"sdicon_b2_m16",variant_filter=SDI_VARS)
build(sd_b1,"sdicon_b1_m16",variant_filter=SDI_VARS)

cat("\n(B) WITHIN-REGION (PN only), origin=PN, within-variant site quartiles:\n")
build(eq_b2[grepl("_pn_",eq_b2)], "conus_b2_pn", origin_by_variant=TRUE)
build(sd_b2[grepl("_pn_",sd_b2)], "sdicon_b2_pn", origin_by_variant=TRUE)
cat("\nDONE controlled.\n")
