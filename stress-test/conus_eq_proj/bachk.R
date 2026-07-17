suppressPackageStartupMessages(library(data.table))
d<-as.data.table(readRDS("/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"))
d<-d[fvs_variant=="NE"]
cat("BA1 q:\n");print(round(quantile(d$BA1,c(.05,.5,.95),na.rm=T),3))
cat("ba_x_rd q:\n");print(round(quantile(d$ba_x_rd,c(.05,.5,.95),na.rm=T),3))
cat("bal_x_rd q:\n");print(round(quantile(d$bal_x_rd,c(.05,.5,.95),na.rm=T),3))
cat("rd_sdimax q:\n");print(round(quantile(d$rd_sdimax,c(.05,.5,.95),na.rm=T),3))
## infer: ba_x_rd / BA1 ~ rd ; compare to rd_sdimax
d[,rd_implied:=ba_x_rd/BA1]
cat("implied rd (ba_x_rd/BA1) q:\n");print(round(quantile(d$rd_implied,c(.05,.5,.95),na.rm=T),3))
## is BA1 ft2/ac or m2/ha? a single tree DBH 7in TPA? check BA1 vs sum check not possible per-tree. 
## relationship: ba_x_rd = BA1 * (rd). And b2 uses BA*0.2296 (ft2ac->m2ha). 
## So if BA1 is m2/ha, ba_x_rd is m2/ha-scale. My dynamic BA is ft2/ac. ratio ft2ac/m2ha = 4.356.
cat("\nmedian BA1:",median(d$BA1,na.rm=T)," (28~ m2/ha if metric; *4.356=122 ft2/ac)\n")
