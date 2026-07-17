suppressMessages(library(data.table))
B <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
v3dir <- file.path(B,"out_conus_eq_v3"); v2dir <- file.path(B,"out_conus_eq")

summ <- function(dir){
  fs <- list.files(dir, pattern="_metrics.csv$", full.names=TRUE)
  rbindlist(lapply(fs, function(f){
    dt <- tryCatch(fread(f, showProgress=FALSE), error=function(e) NULL)
    if(is.null(dt) || !all(c("STAND_CN","VARIANT","CONFIG","TPH","QMD_IN","BA_FT2AC","PROJ_YEAR") %in% names(dt))) return(NULL)
    d <- dt[is.finite(TPH)&is.finite(QMD_IN)&TPH>0&QMD_IN>0, .(STAND_CN,TPH,QMD_IN,PROJ_YEAR)]
    d[, `:=`(x=log(QMD_IN), y=log(TPH))]
    # keep stands that lose stems (final TPH < initial TPH) and have >=5 obs
    chg <- d[, .(loses = TPH[which.max(PROJ_YEAR)] < TPH[which.min(PROJ_YEAR)], n=.N), by=STAND_CN]
    keep <- chg[loses==TRUE & n>=5, STAND_CN]
    d <- d[STAND_CN %in% keep]
    # vectorized per-stand slope = cov(x,y)/var(x)
    sl <- d[, .(slope = { mx<-mean(x); my<-mean(y); sum((x-mx)*(y-my))/sum((x-mx)^2) }), by=STAND_CN]
    sl <- sl[is.finite(slope)]
    ba100 <- dt[PROJ_YEAR==100 & is.finite(BA_FT2AC), mean(BA_FT2AC)]
    data.table(variant=tolower(dt$VARIANT[1]), config=dt$CONFIG[1],
               med_slope=median(sl$slope), n_thin_stands=nrow(sl), yr100_mean_BA=ba100)
  }), fill=TRUE)
}
v3 <- summ(v3dir); v2 <- summ(v2dir)
m <- merge(v2[,.(variant,config,v2_slope=med_slope,v2_ba100=yr100_mean_BA)],
           v3[,.(variant,config,v3_slope=med_slope,v3_ba100=yr100_mean_BA,n_thin=n_thin_stands)],
           by=c("variant","config"))
m[, steepened := v3_slope < v2_slope]
m[, ba_collapse := v3_ba100 < 20]
setorder(m, variant, config)
cat("=== Cross-variant constraint check: v3 vs v2 (variants with v3 done) ===\n")
print(m, digits=3)
cat("\nReineke target slope ~ -1.605\n")
cat("\nv3 did NOT steepen:\n"); print(m[steepened==FALSE,.(variant,config,v2_slope,v3_slope)], digits=3)
cat("\nv3 yr100 BA < 20 (implausible):\n"); print(m[ba_collapse==TRUE,.(variant,config,v3_ba100)], digits=3)
fwrite(m, file.path(v3dir,"v3_vs_v2_constraint_check.csv"))
cat("\nWROTE_CSV_OK rows=",nrow(m),"\n")
