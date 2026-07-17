suppressPackageStartupMessages(library(data.table))
P6 <- "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"
d <- as.data.table(readRDS(P6))
cat("== cspiv6 columns ==\n"); print(grep("DBH|BAL|BA1|TPH|TPA|SDI|bgi|cspi|CR1|FORTYP|EPA|SPCD|fvs_var|plot|YEAR|SICOND|HT1",names(d),value=TRUE,ignore.case=TRUE))
ne <- d[fvs_variant=="NE"]
cat("\n== NE DBH1 quantiles ==\n"); print(round(quantile(ne$DBH1,c(.05,.5,.95),na.rm=TRUE),2))
cat("max DBH1:",max(ne$DBH1,na.rm=TRUE),"\n")
cat("\n== BA1 quantiles (units?) ==\n"); print(round(quantile(ne$BA1,c(.05,.5,.95),na.rm=TRUE),2))
cat("== TPH1 quantiles ==\n"); print(round(quantile(ne$TPH1,c(.05,.5,.95),na.rm=TRUE),1))
cat("== BAL_SW1 q ==\n"); print(round(quantile(ne$BAL_SW1,c(.05,.5,.95),na.rm=TRUE),2))
cat("== bgi q ==\n"); print(round(quantile(ne$bgi,c(.05,.5,.95),na.rm=TRUE),3))
cat("== cspi q ==\n"); print(round(quantile(ne$cspi,c(.05,.5,.95),na.rm=TRUE),3))
cat("\n== dg_obs check: (DBH2-DBH1)/YEARS q ==\n"); print(round(quantile((ne$DBH2-ne$DBH1)/ne$YEARS,c(.05,.5,.95),na.rm=TRUE),3))
cat("\n== plot key cols ==\n"); print(grep("plot_key|STAND_CN|PLT_CN|CN",names(d),value=TRUE,ignore.case=TRUE))
