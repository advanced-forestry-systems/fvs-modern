suppressPackageStartupMessages(library(data.table))
d<-as.data.table(readRDS("/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"))
d<-d[fvs_variant=="NE"]
cat("TPH_UNADJ1 q:\n"); print(round(quantile(d$TPH_UNADJ1,c(.05,.5,.95),na.rm=TRUE),1))
cat("TPA1 q:\n"); print(round(quantile(d$TPA1,c(.05,.5,.95),na.rm=TRUE),1))
cat("TPA_UNADJ1 q:\n"); print(round(quantile(d$TPA_UNADJ1,c(.05,.5,.95),na.rm=TRUE),2))
## per-plot stand BA from pairs using TPA1 vs TPH_UNADJ
d[,baft:=pi/4*(DBH1^2)*(TPA1)/144]
pl<-d[,.(BA_TPA1=sum(baft,na.rm=TRUE), n=.N),by=plot_key]
cat("\nstand BA using TPA1 (ft2/ac) q over plots:\n"); print(round(quantile(pl$BA_TPA1,c(.05,.5,.95),na.rm=TRUE),1))
d[,baft2:=pi/4*(DBH1^2)*(TPH_UNADJ1/2.4710538)/144]
pl2<-d[,.(BA=sum(baft2,na.rm=TRUE)),by=plot_key]
cat("stand BA using TPH_UNADJ1/2.471 (ft2/ac) q:\n"); print(round(quantile(pl2$BA,c(.05,.5,.95),na.rm=TRUE),1))
