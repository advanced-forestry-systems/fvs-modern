suppressPackageStartupMessages(library(data.table))
d <- as.data.table(readRDS("data/conus_remeasurement_pairs_metric_cond_v2.rds")); nm<-names(d)
cat("PLT-CN-ish cols:", paste(grep("CN|PLT|PREV", nm, value=TRUE), collapse=", "), "\n")
pick<-function(a,b) if(a %in% nm) a else b
fortyp<-pick("FORTYPCD_cond1","FORTYPCD"); stdage<-pick("STDAGE_cond1","STDAGE")
pltcn<-pick("PLT_CN_cond1","PLT_CN")
d <- d[fvs_variant=="NE" & is.finite(YEARS)&YEARS>=5&YEARS<=15 &
       is.finite(BA1)&BA1>0&is.finite(TPH1)&TPH1>0 & is.finite(BA2)&BA2>0&is.finite(TPH2)&TPH2>0&is.finite(QMD2)&QMD2>0]
# one row per plot with its PLT_CN + observed t1/t2 + cond fields
plots <- d[, .(PLT_CN=get(pltcn)[1], SICOND=SICOND[1], STDAGE=get(stdage)[1], ASPECT=ASPECT[1], SLOPE=SLOPE[1],
               ELEV=ELEV[1], STATECD=STATECD[1], COUNTYCD=if("COUNTYCD"%in%nm) COUNTYCD[1] else 0L,
               FORTYPCD=get(fortyp)[1], YEARS=round(YEARS[1]),
               BA1=BA1[1], TPH1=TPH1[1], QMD1=if("QMD1"%in%nm) QMD1[1] else NA_real_,
               BA2=BA2[1], TPH2=TPH2[1], QMD2=QMD2[1]), by=plot_key]
plots <- plots[is.finite(PLT_CN) & PLT_CN>0]
set.seed(11); plots <- plots[sample(.N, min(40,.N))]
fwrite(plots, "output/conus/stand_level/ne_bench2_plots.csv")
cat("exported NE plots:", nrow(plots), " sample PLT_CN:", paste(head(plots$PLT_CN,3),collapse=","), "\n")
cat("obs t1 meanBA", round(mean(plots$BA1),1), " obs t2 meanBA", round(mean(plots$BA2),1)," meanYrs",round(mean(plots$YEARS),1),"\n")
