suppressPackageStartupMessages({library(data.table)})
f<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq_greg/conus_eq_ne_greg_metrics.csv"
d<-fread(f,showProgress=FALSE)
cat("cols:",paste(names(d),collapse=","),"\n")
cat("nrow:",nrow(d)," nstands:",uniqueN(d$STAND_CN)," proj_years:",paste(sort(unique(d$PROJ_YEAR)),collapse=","),"\n\n")
mv<-intersect(c("AGB_TONS_AC","BA_FT2AC","QMD_IN","TPH","CCH_MEAN"),names(d))
tr<-d[,lapply(.SD,function(x) round(median(x,na.rm=TRUE),3)),by=PROJ_YEAR,.SDcols=mv][order(PROJ_YEAR)]
cat("=== NE median trajectory by PROJ_YEAR ===\n"); print(tr)
cat("\n=== CCH_MEAN populated? non-NA frac by year (first/last) ===\n")
cch<-d[,.(n=.N,nonNA=sum(is.finite(CCH_MEAN)),frac=round(mean(is.finite(CCH_MEAN)),3)),by=PROJ_YEAR][order(PROJ_YEAR)]
print(cch[c(1,.N)])
cat("\n=== per-stand QMD change yr0 -> yr100 (does QMD grow for MOST stands?) ===\n")
w<-dcast(d[PROJ_YEAR %in% c(0,100)],STAND_CN~PROJ_YEAR,value.var="QMD_IN",fun.aggregate=function(x) x[1])
setnames(w,c("STAND_CN","q0","q100"))
w<-w[is.finite(q0)&is.finite(q100)]
cat(sprintf("  stands: %d ; QMD grew: %.1f%% ; declined: %.1f%% ; median delta=%.2f in\n",
  nrow(w),100*mean(w$q100>w$q0),100*mean(w$q100<w$q0),median(w$q100-w$q0)))
cat("\n=== per-stand TPH yr0->yr100 (mortality severity) ===\n")
wt<-dcast(d[PROJ_YEAR %in% c(0,100)],STAND_CN~PROJ_YEAR,value.var="TPH",fun.aggregate=function(x) x[1])
setnames(wt,c("STAND_CN","t0","t100")); wt<-wt[is.finite(t0)&is.finite(t100)&t0>0]
cat(sprintf("  median TPH retention yr100/yr0 = %.2f ; stands losing >90%% of TPH: %.1f%%\n",
  median(wt$t100/wt$t0), 100*mean(wt$t100/wt$t0 < 0.10)))
