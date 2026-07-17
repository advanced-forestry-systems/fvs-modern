suppressPackageStartupMessages(library(data.table))
LK <- readRDS("/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/greg_emt_td_lookup.rds")
setDT(LK); LK[,STAND_CN:=as.character(STAND_CN)]
SI_DIR <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
vf <- list.files(SI_DIR, pattern="^standinit_[A-Z]+\\.csv$", full.names=TRUE)
vf <- vf[!grepl("_BLANK", vf)]
elev <- rbindlist(lapply(vf, function(f){
  s <- fread(f, colClasses=list(character="STAND_CN"), select=c("STAND_CN","ELEVFT"))
  s[,STAND_CN:=sub("\\..*$","",STAND_CN)]; s}), fill=TRUE)
elev <- unique(elev[!is.na(STAND_CN)&STAND_CN!=""], by="STAND_CN")
out <- merge(LK[,.(STAND_CN,EMT,TD)], elev, by="STAND_CN", all.x=TRUE)
setnames(out,"ELEVFT","ELEV")
out[is.na(ELEV), ELEV := 0]
setcolorder(out, c("STAND_CN","EMT","TD","ELEV"))
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/greg_climate_perstand.csv"
fwrite(out, OUT)
cat("Wrote", OUT, "n=", nrow(out), "\n")
cat("cols:", paste(names(out),collapse=","), "\n")
cat("EMT median", round(median(out$EMT,na.rm=TRUE),2), " TD median", round(median(out$TD,na.rm=TRUE),2),
    " ELEV median", round(median(out$ELEV,na.rm=TRUE),1), " ELEV>0 frac", round(mean(out$ELEV>0),3), "\n")
print(head(out,3))
