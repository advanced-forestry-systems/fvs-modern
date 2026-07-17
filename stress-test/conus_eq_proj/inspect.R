suppressPackageStartupMessages(library(data.table))
ROOT<-"/users/PUOM0008/crsfaaron/fvs-conus"
PAIRS<-file.path(ROOT,"data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds")
d<-as.data.table(readRDS(PAIRS))
cat("PAIRS ncols:",ncol(d),"nrows:",nrow(d),"\n")
cat("colnames:\n"); print(names(d))
cat("\nNE rows:",nrow(d[fvs_variant=='NE']),"\n")
ne<-d[fvs_variant=='NE']
cat("key candidate cols sample (NE):\n")
kcols<-grep("CN|plot_key|PLT|STAND|cond",names(d),value=TRUE,ignore.case=TRUE)
print(kcols)
print(head(ne[,..kcols],3))
