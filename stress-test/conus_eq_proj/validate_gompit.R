suppressMessages({library(data.table)})
D <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/ne_gompit_smoke"
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
g <- fread(file.path(D,"conus_eq_ne_conus_b2_gompit_metrics.csv"))
l <- fread(file.path(D,"conus_eq_ne_conus_b2_metrics.csv"))
g[,STAND_CN:=as.character(STAND_CN)]; l[,STAND_CN:=as.character(STAND_CN)]

g0 <- g[PROJ_YEAR==0]; l0 <- l[PROJ_YEAR==0]
common <- intersect(g0$STAND_CN, l0$STAND_CN)
g0 <- g0[STAND_CN %in% common][order(STAND_CN)]; l0 <- l0[STAND_CN %in% common][order(STAND_CN)]
max_ba <- max(abs(g0$BA_FT2AC - l0$BA_FT2AC), na.rm=TRUE)
max_qmd<- max(abs(g0$QMD_IN  - l0$QMD_IN ), na.rm=TRUE)
max_tph<- max(abs(g0$TPH     - l0$TPH    ), na.rm=TRUE)
cat(sprintf("[Y0 IDENTITY] n=%d  max|dBA|=%.6g  max|dQMD|=%.6g  max|dTPH|=%.6g\n",
            length(common), max_ba, max_qmd, max_tph))

cc <- g[is.finite(CCH_MEAN)]
cat(sprintf("[CCH SANITY] mean=%.3f median=%.3f range[%.3f,%.3f] pct>8=%.1f%%\n",
   mean(cc$CCH_MEAN), median(cc$CCH_MEAN), min(cc$CCH_MEAN), max(cc$CCH_MEAN), 100*mean(cc$CCH_MEAN>8)))

summ <- function(dt,tag) dt[PROJ_YEAR %in% c(0,25,50,75,100),
   .(arm=tag, BA=round(mean(BA_FT2AC,na.rm=TRUE),1), QMD=round(mean(QMD_IN,na.rm=TRUE),2),
     TPH=round(mean(TPH,na.rm=TRUE),0), n=.N), by=PROJ_YEAR]
S <- rbind(summ(g,"gompit"), summ(l,"logit"))[order(PROJ_YEAR,arm)]
cat("\n[TRAJECTORIES]\n"); print(S)

co <- fread("/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv")
ov <- fread("/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/mortality_compare_overall.csv")
samp <- fread("/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/cch_validation_sample.csv")
samp <- samp[is.finite(CR)&is.finite(CCH1)&CR>0&CR<=1&CCH1>=0]
samp[,SPCD:=as.integer(SPCD)]; co[,SPCD:=as.integer(SPCD)]
m <- merge(samp, co[,.(SPCD,b0,b1,b2,b3,b4)], by="SPCD")
m[, eta := b0 + b1*(CR+0.01)^b2 + b3*ifelse(CCH1>0,CCH1^b4,0)]
m[, eta := pmin(pmax(eta,-30),30)]
m[, surv_ann := exp(-exp(eta))]
cat(sprintf("\n[GOMPIT SURV REPRO] banked-coef annual survival on %d trees: mean=%.4f\n", nrow(m), mean(m$surv_ann)))
cat(sprintf("  Greg validated overall: obs_surv=%.4f pred_surv_new=%.4f pred_surv_base=%.4f (AUC_new=%.3f AUC_base=%.3f)\n",
            ov$obs_surv, ov$pred_surv_new, ov$pred_surv_base, ov$auc_new, ov$auc_base))

vs <- data.table(
  metric=c("Y0_max_abs_dBA","Y0_max_abs_dQMD","Y0_max_abs_dTPH","n_stands_common",
           "CCH_mean","CCH_median","CCH_max","CCH_pct_gt8",
           "gompit_Y100_BA","gompit_Y100_TPH","logit_Y100_BA","logit_Y100_TPH",
           "gompit_surv_ann_mean_on_sample","greg_pred_surv_new","greg_obs_surv","greg_auc_new"),
  value=c(max_ba,max_qmd,max_tph,length(common),
          mean(cc$CCH_MEAN),median(cc$CCH_MEAN),max(cc$CCH_MEAN),100*mean(cc$CCH_MEAN>8),
          g[PROJ_YEAR==100,mean(BA_FT2AC,na.rm=TRUE)], g[PROJ_YEAR==100,mean(TPH,na.rm=TRUE)],
          l[PROJ_YEAR==100,mean(BA_FT2AC,na.rm=TRUE)], l[PROJ_YEAR==100,mean(TPH,na.rm=TRUE)],
          mean(m$surv_ann), ov$pred_surv_new, ov$obs_surv, ov$auc_new))
fwrite(vs, file.path(OUT,"validation_summary_NE_gompit.csv"))
cat("\nWrote validation_summary_NE_gompit.csv\n")

png(file.path(OUT,"thumb_NE_gompit.png"), width=1100, height=380, res=110)
par(mfrow=c(1,3), mar=c(4,4,2,1))
ga <- g[, .(BA=mean(BA_FT2AC,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE)), by=PROJ_YEAR][order(PROJ_YEAR)]
la <- l[, .(BA=mean(BA_FT2AC,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE)), by=PROJ_YEAR][order(PROJ_YEAR)]
plot(ga$PROJ_YEAR,ga$BA,type="l",col="firebrick",lwd=2,ylim=range(0,ga$BA,la$BA),xlab="year",ylab="BA ft2/ac",main="BA (NE)")
lines(la$PROJ_YEAR,la$BA,col="steelblue",lwd=2); legend("topleft",c("gompit","logit"),col=c("firebrick","steelblue"),lwd=2,bty="n",cex=0.9)
plot(ga$PROJ_YEAR,ga$TPH,type="l",col="firebrick",lwd=2,ylim=range(0,ga$TPH,la$TPH),xlab="year",ylab="TPH",main="TPH (NE)")
lines(la$PROJ_YEAR,la$TPH,col="steelblue",lwd=2)
hist(cc$CCH_MEAN,breaks=30,col="grey70",border="white",xlab="CCH (gompit scale)",main="stand-mean CCH",xlim=c(0,max(8,max(cc$CCH_MEAN))))
abline(v=8,lty=2,col="red"); dev.off()
cat("Wrote thumb_NE_gompit.png\n")
