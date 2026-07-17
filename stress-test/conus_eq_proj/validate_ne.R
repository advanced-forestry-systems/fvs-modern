#!/usr/bin/env Rscript
## validate_ne.R -- compare conus_b2 / conus_b1 projection arms vs engine NE
## (default & calibrated) on AGB, plus the projector's own BA/QMD/TPH. Emits a
## compact summary CSV + thumbnail PNGs (<=800px) to scratch.
suppressPackageStartupMessages({library(data.table)})
OUTD<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
ENG<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_wo1"
b2<-fread(file.path(OUTD,"conus_eq_ne_conus_b2_metrics.csv"))
b1<-fread(file.path(OUTD,"conus_eq_ne_conus_b1_metrics.csv"))
## engine NE: concat all batch files, keep CONFIG default/calibrated
ef<-list.files(ENG,pattern="^conus_ne_b[0-9]+\\.csv$",full.names=TRUE)
eng<-rbindlist(lapply(ef,fread),fill=TRUE)
eng<-eng[CONFIG %in% c("default","calibrated")]
## --- summary table: year 0/50/100, per arm, mean+sd of AGB; BA/QMD/TPH for eq arms ---
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
fwrite(S,file.path(OUTD,"validation_summary_NE.csv"))
cat("=== NE validation summary (year 0/50/100) ===\n"); print(S,nrow=40)

## --- thumbnails (<=800px) ---
agg<-function(dt,arm,col){ d<-dt[,.(v=mean(get(col),na.rm=TRUE)),by=PROJ_YEAR]; d[,arm:=arm]; d }
make_png<-function(file,title,series,ylab){
  png(file,width=720,height=460,res=110)
  par(mar=c(4,4,2.5,1))
  cols<-c(conus_b2="#1b9e77",conus_b1="#d95f02",engine_default="#7570b3",engine_calibrated="#e7298a")
  xr<-range(unlist(lapply(series,function(s) s$PROJ_YEAR)))
  yr<-range(unlist(lapply(series,function(s) s$v[is.finite(s$v)])))
  plot(NA,xlim=xr,ylim=yr,xlab="Projection year",ylab=ylab,main=title)
  for(s in series){ a<-s$arm[1]; lines(s$PROJ_YEAR,s$v,col=cols[a],lwd=2); points(s$PROJ_YEAR,s$v,col=cols[a],pch=19,cex=.5) }
  legend("topleft",legend=names(cols),col=cols,lwd=2,bty="n",cex=.8)
  dev.off()
}
## AGB: all 4 arms
make_png(file.path(OUTD,"thumb_AGB_NE.png"),"NE AGB trajectory (stand mean)",
  list(agg(b2,"conus_b2","AGB_TONS_AC"),agg(b1,"conus_b1","AGB_TONS_AC"),
       agg(eng[CONFIG=="default"],"engine_default","AGB_TONS_AC"),
       agg(eng[CONFIG=="calibrated"],"engine_calibrated","AGB_TONS_AC")),"AGB (tons/ac)")
## BA/QMD/TPH: eq arms only (engine arms have no BA in this output)
png(file.path(OUTD,"thumb_BA_QMD_TPH_NE.png"),width=800,height=300,res=100)
par(mfrow=c(1,3),mar=c(4,4,2,1))
cols<-c(conus_b2="#1b9e77",conus_b1="#d95f02")
for(cc in c("BA_FT2AC","QMD_IN","TPH")){
  a2<-agg(b2,"conus_b2",cc); a1<-agg(b1,"conus_b1",cc)
  yr<-range(c(a2$v,a1$v),na.rm=TRUE)
  plot(a2$PROJ_YEAR,a2$v,type="l",col=cols[1],lwd=2,ylim=yr,xlab="year",ylab=cc,main=cc)
  lines(a1$PROJ_YEAR,a1$v,col=cols[2],lwd=2)
  if(cc=="BA_FT2AC") legend("bottomright",names(cols),col=cols,lwd=2,bty="n",cex=.8)
}
dev.off()
cat("\nWrote: validation_summary_NE.csv, thumb_AGB_NE.png, thumb_BA_QMD_TPH_NE.png\n")
