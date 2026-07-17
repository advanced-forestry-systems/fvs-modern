library(data.table)
D<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
v4<-fread(file.path(D,"out_v4_seed","conus_eq_ne_conus_b2_metrics.csv"),colClasses=list(character="STAND_CN"))
v6<-fread(file.path(D,"out_v6_seed","conus_eq_ne_conus_b2_metrics.csv"),colClasses=list(character="STAND_CN"))
# pick a stand present in both with many cycles
cn<-v4[,.N,by=STAND_CN][order(-N)][1,STAND_CN]
a<-v4[STAND_CN==cn,.(PROJ_YEAR,QMD_IN,HT4_DOM=HT_M_DOM,HT4_MEAN=HT_M_MEAN)]
b<-v6[STAND_CN==cn,.(PROJ_YEAR,HT6_DOM=HT_M_DOM,HT6_MEAN=HT_M_MEAN)]
m<-merge(a,b,by="PROJ_YEAR")
m<-m[PROJ_YEAR %in% seq(0,100,10)]
cat("STAND",cn,"\n"); print(m)
# per-stand monotonicity of HT_M_DOM over full 0..100 for both
chk<-function(dt,col){d<-dt[order(PROJ_YEAR)];h<-d[[col]];sum(diff(h)< -1e-6)}
cat(sprintf("v4 HT_DOM downsteps=%d  v6 HT_DOM downsteps=%d  (over %d cycles)\n",
  chk(v4[STAND_CN==cn],"HT_M_DOM"),chk(v6[STAND_CN==cn],"HT_M_DOM"),v4[STAND_CN==cn,.N]))
# fraction of stands with monotone HT_DOM to year 100
mono_frac<-function(dt){dt2<-dt[PROJ_YEAR<=100];dt2[order(STAND_CN,PROJ_YEAR)][,.(mono=all(diff(HT_M_DOM)>=-1e-6)),by=STAND_CN][,mean(mono,na.rm=TRUE)]}
cat(sprintf("frac stands monotone HT_DOM to yr100: v4=%.3f v6=%.3f\n",mono_frac(v4),mono_frac(v6)))
cat(sprintf("frac stands monotone HT_MEAN to yr100: v4=%.3f v6=%.3f\n",
  v4[PROJ_YEAR<=100][order(STAND_CN,PROJ_YEAR)][,.(m=all(diff(HT_M_MEAN)>=-1e-6)),by=STAND_CN][,mean(m,na.rm=TRUE)],
  v6[PROJ_YEAR<=100][order(STAND_CN,PROJ_YEAR)][,.(m=all(diff(HT_M_MEAN)>=-1e-6)),by=STAND_CN][,mean(m,na.rm=TRUE)]))
