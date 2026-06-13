#!/usr/bin/env Rscript
# mortality_by_size.R -- is large-tree (senescence) mortality real in FIA? Observed annual
# mortality by absolute DBH and by RELATIVE size (DBH/species max), to test the U-shape premise.
suppressPackageStartupMessages({ library(data.table) })
d <- as.data.table(readRDS("data/conus_remeasurement_pairs_metric_cond_v2.rds")); nm<-names(d)
st2 <- intersect(c("TREESTATUS2","STATUS2"), nm)[1]; st1 <- intersect(c("TREESTATUS1","STATUS1"), nm)[1]
d <- d[is.finite(YEARS)&YEARS>=1&YEARS<=20 & is.finite(DBH1)&DBH1>=2.54 & get(st1)==1 & get(st2) %in% c(1,2)]
d[, died := as.integer(get(st2)==2)]
d[, surv_ann := (1-died)^(1/1)]   # placeholder
# annual mortality per tree via exposure: annual survival p so that p^years = survived(0/1) is degenerate;
# instead aggregate by bin: annual mort = 1 - (sum_survived / n)^(1/mean_years)
tr <- intersect(c("max_dbh_cm","max_dbh"), nm)
mdcol <- if(length(tr)) tr[1] else NA
cat("rows:", nrow(d), " overall raw dead frac:", round(mean(d$died),4), " mean yrs:", round(mean(d$YEARS),1), "\n\n")
binmort <- function(dd, by, labels){
  dd[, .(n=.N, surv=mean(1-died), yrs=mean(YEARS),
         ann_mort = 1 - (mean(1-died))^(1/mean(YEARS))), by=by][order(get(by))]
}
cat("== annual mortality by ABSOLUTE DBH class (cm) ==\n")
d[, dbh_cls := cut(DBH1, c(2.5,10,20,30,40,50,70,100,300), right=FALSE)]
print(d[, .(n=.N, ann_mort=round(1-(mean(1-died))^(1/mean(YEARS)),4)), by=dbh_cls][order(dbh_cls)])
if(!is.na(mdcol)){
  d[, relsz := DBH1/get(mdcol)]
  d <- d[is.finite(relsz)&relsz>0&relsz<1.5]
  d[, rel_cls := cut(relsz, c(0,0.2,0.4,0.6,0.8,1.0,1.5), right=FALSE)]
  cat("\n== annual mortality by RELATIVE size (DBH/species max) ==\n")
  print(d[, .(n=.N, ann_mort=round(1-(mean(1-died))^(1/mean(YEARS)),4)), by=rel_cls][order(rel_cls)])
} else cat("\n(no species max_dbh column found for relative-size analysis)\n")
cat("\nDONE_MORT_SIZE\n")
