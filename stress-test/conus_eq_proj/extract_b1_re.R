suppressPackageStartupMessages({library(data.table)})
OUT<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/b1_re_means.rds"
B1FIT<-"/users/PUOM0008/crsfaaron/fvs-conus/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds"
B1META<-"/users/PUOM0008/crsfaaron/fvs-conus/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_meta.rds"
cat("loading 4.7GB b1 fit (one-time)...\n"); fit<-readRDS(B1FIT)
pull<-function(v) as.numeric(fit$summary(v,"mean")$mean)
RE<-list(z_L1=pull("z_L1"), z_L2=pull("z_L2"), z_L3=pull("z_L3"))
saveRDS(RE,OUT)
cat("z_L1 n:",length(RE$z_L1)," z_L2 n:",length(RE$z_L2)," z_L3 n:",length(RE$z_L3),"\n")
cat("wrote:",OUT,"\n")
