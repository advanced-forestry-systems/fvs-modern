#!/usr/bin/env Rscript
# mortality_relsize.R -- the decisive senescence test: observed annual mortality by RELATIVE
# size (DBH/species max), joining species max_dbh from traits. Senescence = mortality rising
# for trees near their species maximum size.
suppressPackageStartupMessages({ library(data.table) })
d <- as.data.table(readRDS("data/conus_remeasurement_pairs_metric_cond_v2.rds")); nm<-names(d)
tr <- as.data.table(readRDS("traits/species_traits_v2.rds"))
mdcol <- intersect(c("max_dbh_cm","max_dbh"), names(tr))[1]
cat("traits max_dbh col:", mdcol, "\n")
md <- tr[, .(SPCD, max_dbh=get(mdcol))]
st2 <- intersect(c("TREESTATUS2","STATUS2"), nm)[1]; st1 <- intersect(c("TREESTATUS1","STATUS1"), nm)[1]
d <- d[is.finite(YEARS)&YEARS>=1&YEARS<=20 & is.finite(DBH1)&DBH1>=2.54 & get(st1)==1 & get(st2) %in% c(1,2)]
d[, died := as.integer(get(st2)==2)]
d <- merge(d, md, by="SPCD", all.x=TRUE)
d <- d[is.finite(max_dbh)&max_dbh>0]
d[, relsz := DBH1/max_dbh]; d <- d[relsz>0 & relsz<1.3]
cat("rows w/ max_dbh:", nrow(d), " median relsz:", round(median(d$relsz),2), "\n\n")
d[, rel_cls := cut(relsz, c(0,0.15,0.3,0.45,0.6,0.75,0.9,1.3), right=FALSE)]
cat("== annual mortality by RELATIVE size (DBH / species max DBH) ==\n")
print(d[, .(n=.N, mean_dbh=round(mean(DBH1),1), ann_mort=round(1-(mean(1-died))^(1/mean(YEARS)),4)), by=rel_cls][order(rel_cls)])
cat("\nSenescence present if ann_mort rises in the top relative-size classes (>0.75).\n")
cat("DONE_RELSIZE\n")
