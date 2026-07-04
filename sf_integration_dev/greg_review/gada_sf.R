suppressPackageStartupMessages(library(data.table))
sit <- fread("/users/PUOM0008/crsfaaron/SiteIndex/NA_SITREE.csv")
sit <- sit[!is.na(HT) & !is.na(AGEDIA) & HT>0 & AGEDIA>=10 & AGEDIA<=200]
sf <- sit[SPCD %in% c(12,95,97)]   # balsam fir, black spruce, red spruce
cat("SF n =", nrow(sf), "\n")
cat("HT units check (range):", paste(round(range(sf$HT),1),collapse=" - "), "\n")
fit <- nls(HT ~ b1*(1-exp(-b2*AGEDIA))^b3, data=sf,
           start=list(b1=25,b2=0.03,b3=1.2),
           control=nls.control(maxiter=300, warnOnly=TRUE))
cf <- coef(fit)
cat(sprintf("SF GADA: b1=%.3f b2=%.5f b3=%.4f\n", cf["b1"],cf["b2"],cf["b3"]))
# SI50 per record, then quantiles of plot-level site asymptote b1_site = SI50/(1-exp(-b2*50))^b3
sf[, SI50 := HT*((1-exp(-cf["b2"]*50))/(1-exp(-cf["b2"]*AGEDIA)))^cf["b3"]]
q <- quantile(sf$SI50, c(.1,.25,.5,.75,.9), na.rm=TRUE)
cat("SI50 quantiles (m):", paste(round(q,1),collapse=" "), "\n")
b1site <- q / (1-exp(-cf["b2"]*50))^cf["b3"]
cat("implied b1 (asymptote) quantiles (m):", paste(round(b1site,1),collapse=" "), "\n")
cat("DONE_GADA\n")
