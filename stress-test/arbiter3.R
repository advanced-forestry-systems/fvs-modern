suppressPackageStartupMessages(library(data.table))
p<-readRDS("/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds")
p<-as.data.table(p)
cat("cols:\n"); print(grep("BA|CN|PLT|STAND|STATE",names(p),value=TRUE,ignore.case=TRUE))
# BA1 in m2/ha -> ft2/ac x4.356; report by variant-state if variant col present
vcol<-grep("^VARIANT$|variant",names(p),value=TRUE,ignore.case=TRUE)[1]
ba1<-grep("^BA1$|BA_1|baa1",names(p),value=TRUE,ignore.case=TRUE)[1]
cat("vcol=",vcol," ba1=",ba1,"\n")
if(!is.na(ba1)){
  p[,BA1_ft2ac:=get(ba1)*4.356]
  if(!is.na(vcol)){
    print(p[toupper(get(vcol)) %in% c("CR","PN","WC"), .(mean_BA1_ft2ac=mean(BA1_ft2ac,na.rm=TRUE), median=median(BA1_ft2ac,na.rm=TRUE), n=.N), by=.(VAR=toupper(get(vcol)))])
  } else {
    cat("no variant col; overall BA1 mean ft2ac=",mean(p$BA1_ft2ac,na.rm=TRUE),"\n")
  }
}
