#!/usr/bin/env Rscript
## extract_emt_td.R -- per-STAND_CN EMT and TD (=MWMT-MCMT) for the faithful Greg arm.
## Greg's exact ClimateNA 1991-2020 normals (1km, LAEA). EMT = Normal_1991_2020_EMT.tif (degC),
## TD  = Normal_1991_2020_MWMT.tif - Normal_1991_2020_MCMT.tif (ClimateNA continentality).
## Stand coords are the standinit public/fuzzed LAT/LON (acceptable at 1km normal resolution).
## Output keeps ONLY STAND_CN, EMT, TD (+ELEVFT passthrough already in standinit). No coords leave.
suppressPackageStartupMessages({library(data.table); library(terra)})
RAS<-"/users/PUOM0008/crsfaaron/SiteIndex/rasters/ClimateNA/Normal_1991_2020_bioclim"
SI_DIR<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
OUT<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/greg_emt_td_lookup.rds"
emt<-rast(file.path(RAS,"Normal_1991_2020_EMT.tif"))
mwmt<-rast(file.path(RAS,"Normal_1991_2020_MWMT.tif"))
mcmt<-rast(file.path(RAS,"Normal_1991_2020_MCMT.tif"))
rcrs<-crs(emt)
vfiles<-list.files(SI_DIR,pattern="^standinit_[A-Z]+\\.csv$",full.names=TRUE)
vfiles<-vfiles[!grepl("_BLANK",vfiles)]
all<-rbindlist(lapply(vfiles,function(f){
  s<-fread(f,colClasses=list(character=c("STAND_CN")),select=c("STAND_CN","LATITUDE","LONGITUDE"))
  s[,STAND_CN:=sub("\\..*$","",STAND_CN)]; s[,src:=basename(f)]; s}),fill=TRUE)
all<-all[is.finite(LATITUDE)&is.finite(LONGITUDE)&!is.na(STAND_CN)&STAND_CN!=""]
u<-unique(all,by="STAND_CN")
cat(sprintf("Unique STAND_CN with coords: %d (from %d variant files)\n",nrow(u),length(vfiles)))
pts<-vect(data.frame(lon=u$LONGITUDE,lat=u$LATITUDE),geom=c("lon","lat"),crs="EPSG:4326")
pts<-project(pts,rcrs)
u[,EMT:=terra::extract(emt,pts)[,2]]
u[,MWMT:=terra::extract(mwmt,pts)[,2]]
u[,MCMT:=terra::extract(mcmt,pts)[,2]]
u[,TD:=MWMT-MCMT]
off<-u[is.na(EMT)|is.na(TD)]
cat(sprintf("Off-raster / NA stands: %d (%.2f%%)\n",nrow(off),100*nrow(off)/nrow(u)))
cat("EMT range:",paste(round(range(u$EMT,na.rm=TRUE),2),collapse=" .. "),"median",round(median(u$EMT,na.rm=TRUE),2),"\n")
cat("TD  range:",paste(round(range(u$TD,na.rm=TRUE),2),collapse=" .. "),"median",round(median(u$TD,na.rm=TRUE),2),"\n")
look<-u[,.(STAND_CN,EMT,MWMT,MCMT,TD)]
saveRDS(look,OUT)
cat("Wrote",OUT,"n=",nrow(look),"\n")
## per-variant coverage
cat("\nPer-variant n stands with EMT/TD (first stand_cn match):\n")
for(f in vfiles){ s<-fread(f,colClasses=list(character=c("STAND_CN")),select=c("STAND_CN")); s[,STAND_CN:=sub("\\..*$","",STAND_CN)]
  m<-look[match(unique(s$STAND_CN),STAND_CN)]; v<-sub("standinit_","",sub("\\.csv","",basename(f)))
  cat(sprintf("  %-4s n=%6d  EMT-cov=%.1f%%\n",v,length(unique(s$STAND_CN)),100*mean(is.finite(m$EMT)))) }
