#!/usr/bin/env Rscript
## validate_ne_v2.R -- compare conus_b2 / conus_b1 (v2 projector, 3 fixes) vs engine
## NE (default & calibrated) on IDENTICAL STAND_CNs. Engine output has AGB only;
## projector arms also carry BA/QMD/TPH. Emits validation_summary_NE_v2.csv +
## thumbnails. Comparison is restricted to the stand set the projector produced
## (which is a subset of the engine's stands -> identical-stand comparison).
suppressPackageStartupMessages({library(data.table)})
OUTD<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
ENG <-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_wo1"
SUB <-commandArgs(trailingOnly=TRUE); SUBDIR<-if(length(SUB)) SUB[1] else "."
b2<-fread(file.path(OUTD,SUBDIR,"conus_eq_ne_conus_b2_metrics.csv"),colClasses=list(character="STAND_CN"))
b1<-fread(file.path(OUTD,SUBDIR,"conus_eq_ne_conus_b1_metrics.csv"),colClasses=list(character="STAND_CN"))
ef<-list.files(ENG,pattern="^conus_ne_b[0-9]+\\.csv$",full.names=TRUE)
eng<-rbindlist(lapply(ef,function(f) fread(f,colClasses=list(character="STAND_CN"))),fill=TRUE)
eng<-eng[CONFIG %in% c("default","calibrated")]
## identical-stand set: stands present in BOTH projector arms
common<-intersect(unique(b2$STAND_CN),unique(b1$STAND_CN))
common<-intersect(common,unique(eng$STAND_CN))
cat("identical-stand comparison set:",length(common),"stands\n")
b2<-b2[STAND_CN %in% common]; b1<-b1[STAND_CN %in% common]; eng<-eng[STAND_CN %in% common]
yrs<-c(0,50,100)
sumrow<-function(dt,arm,hasmetrics){
  s<-dt[PROJ_YEAR %in% yrs, .(
     AGB_mean=mean(AGB_TONS_AC,na.rm=TRUE), AGB_sd=sd(AGB_TONS_AC,na.rm=TRUE),
     BA_mean=if(hasmetrics) mean(BA_FT2AC,na.rm=TRUE) else NA_real_,
     QMD_mean=if(hasmetrics) mean(QMD_IN,na.rm=TRUE) else NA_real_,
     TPH_mean=if(hasmetrics) mean(TPH,na.rm=TRUE) else NA_real_, n=.N), by=PROJ_YEAR]
  s[,arm:=arm]; s}
S<-rbindlist(list(
  sumrow(b2,"conus_b2",TRUE),
  sumrow(b1,"conus_b1",TRUE),
  sumrow(eng[CONFIG=="default"],"engine_default",FALSE),
  sumrow(eng[CONFIG=="calibrated"],"engine_calibrated",FALSE)),fill=TRUE)
setcolorder(S,c("arm","PROJ_YEAR","n","AGB_mean","AGB_sd","BA_mean","QMD_mean","TPH_mean"))
S<-S[order(arm,PROJ_YEAR)]
fwrite(S,file.path(OUTD,"validation_summary_NE_v2.csv"))
cat("=== NE v2 validation summary (identical stands; year 0/50/100) ===\n"); print(S,nrow=40)

## --- explicit GO/NO-GO diagnostics ---
gv<-function(dt,col,yr) dt[PROJ_YEAR==yr,mean(get(col),na.rm=TRUE)]
agb0_b2<-gv(b2,"AGB_TONS_AC",0); agb0_e<-gv(eng[CONFIG=="default"],"AGB_TONS_AC",0)
agb0_b1<-gv(b1,"AGB_TONS_AC",0)
cat(sprintf("\nYEAR-0 AGB: b2=%.2f b1=%.2f engine=%.2f | b2 dev=%.1f%% b1 dev=%.1f%%\n",
  agb0_b2,agb0_b1,agb0_e,100*(agb0_b2-agb0_e)/agb0_e,100*(agb0_b1-agb0_e)/agb0_e))
agb100_b2<-gv(b2,"AGB_TONS_AC",100); agb100_ed<-gv(eng[CONFIG=="default"],"AGB_TONS_AC",100); agb100_ec<-gv(eng[CONFIG=="calibrated"],"AGB_TONS_AC",100)
cat(sprintf("YEAR-100 AGB: b2=%.2f b1=%.2f engine_def=%.2f engine_cal=%.2f | b2/eng=%.2fx\n",
  agb100_b2,gv(b1,"AGB_TONS_AC",100),agb100_ed,agb100_ec,agb100_b2/agb100_ed))
## monotone QMD & BA plateau check (b2)
qmd<-b2[,.(q=mean(QMD_IN,na.rm=TRUE)),by=PROJ_YEAR][order(PROJ_YEAR)]
cat("QMD monotone up (b2):",all(diff(qmd$q)>=-1e-6),"\n")
ba<-b2[,.(b=mean(BA_FT2AC,na.rm=TRUE)),by=PROJ_YEAR][order(PROJ_YEAR)]
cat(sprintf("BA y50=%.1f y100=%.1f (plateau, not blowing up): %s\n",ba[PROJ_YEAR==50,b],ba[PROJ_YEAR==100,b],ba[PROJ_YEAR==100,b] < 1.5*ba[PROJ_YEAR==50,b]))
tph<-b2[,.(t=mean(TPH,na.rm=TRUE)),by=PROJ_YEAR][order(PROJ_YEAR)]
cat(sprintf("TPH y0=%.0f y50=%.0f y100=%.0f (declines then stabilizes w/ ingrowth)\n",tph[PROJ_YEAR==0,t],tph[PROJ_YEAR==50,t],tph[PROJ_YEAR==100,t]))

## --- thumbnails ---
agg<-function(dt,arm,col){ d<-dt[,.(v=mean(get(col),na.rm=TRUE)),by=PROJ_YEAR]; d[,arm:=arm]; d }
cols<-c(conus_b2="#1b9e77",conus_b1="#d95f02",engine_default="#7570b3",engine_calibrated="#e7298a")
png(file.path(OUTD,"thumb_AGB_NE_v2.png"),width=720,height=460,res=110); par(mar=c(4,4,2.5,1))
series<-list(agg(b2,"conus_b2","AGB_TONS_AC"),agg(b1,"conus_b1","AGB_TONS_AC"),
             agg(eng[CONFIG=="default"],"engine_default","AGB_TONS_AC"),
             agg(eng[CONFIG=="calibrated"],"engine_calibrated","AGB_TONS_AC"))
xr<-range(unlist(lapply(series,function(s)s$PROJ_YEAR))); yr<-range(unlist(lapply(series,function(s)s$v[is.finite(s$v)])))
plot(NA,xlim=xr,ylim=yr,xlab="Projection year",ylab="AGB (tons/ac)",main="NE AGB (identical stands, v2)")
for(s in series){a<-s$arm[1];lines(s$PROJ_YEAR,s$v,col=cols[a],lwd=2);points(s$PROJ_YEAR,s$v,col=cols[a],pch=19,cex=.5)}
legend("topleft",legend=names(cols),col=cols,lwd=2,bty="n",cex=.8); dev.off()
png(file.path(OUTD,"thumb_BA_QMD_TPH_NE_v2.png"),width=800,height=300,res=100); par(mfrow=c(1,3),mar=c(4,4,2,1))
cc2<-c(conus_b2="#1b9e77",conus_b1="#d95f02")
for(cc in c("BA_FT2AC","QMD_IN","TPH")){a2<-agg(b2,"conus_b2",cc);a1<-agg(b1,"conus_b1",cc);yr<-range(c(a2$v,a1$v),na.rm=TRUE)
  plot(a2$PROJ_YEAR,a2$v,type="l",col=cc2[1],lwd=2,ylim=yr,xlab="year",ylab=cc,main=cc);lines(a1$PROJ_YEAR,a1$v,col=cc2[2],lwd=2)
  if(cc=="BA_FT2AC") legend("bottomright",names(cc2),col=cc2,lwd=2,bty="n",cex=.8)}
dev.off()
cat("\nWrote: validation_summary_NE_v2.csv, thumb_AGB_NE_v2.png, thumb_BA_QMD_TPH_NE_v2.png\n")
