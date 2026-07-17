#!/usr/bin/env Rscript
# Parse bakuzis_check.R output txt files into a tidy per-(arm,origin,site) Reineke + per-(arm,origin) Eichhorn table.
suppressPackageStartupMessages(library(data.table))
DIR <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/bakuzis_singlecohort"
arms <- c(unconstrained="conus_b2", v3_constrained="sdicon_b2", v4_constrained="v4_b2",
          engine_derivedHT="engine_default", engine_realHT="engine_default_ht")
origins <- c("all","ne","sn","pn")
rein <- list(); eich <- list(); flagL <- list(); htm <- list()
for(an in names(arms)){ b <- arms[[an]]
  for(o in origins){
    f <- file.path(DIR, sprintf("%s_sc_%s_check.txt", b, o))
    if(!file.exists(f)) next
    L <- readLines(f, warn=FALSE)
    # Reineke
    sl <- grep("site [0-9]: slope=", L, value=TRUE)
    for(ln in sl){
      s  <- as.integer(sub(".*site ([0-9]):.*","\\1", ln))
      v  <- as.numeric(sub(".*slope=(-?[0-9.]+).*","\\1", ln))
      ok <- grepl("OK", ln)
      rein[[length(rein)+1]] <- data.table(arm=an, build=b, origin=o, site=s, slope=v, reineke_ok=ok)
    }
    # Eichhorn
    ev <- grep("mean CV=", L, value=TRUE)
    if(length(ev)){ cv <- as.numeric(sub(".*mean CV=([0-9.]+)%.*","\\1", ev[1])); eok <- grepl("OK", ev[1])
      eich[[length(eich)+1]] <- data.table(arm=an, build=b, origin=o, eichhorn_cv=cv, eichhorn_ok=eok) }
    # flags total
    fl <- grep("^FLAGS:", L, value=TRUE)
    if(length(fl)){ nf <- as.integer(sub("FLAGS:\\s*([0-9]+).*","\\1", fl[1]))
      flagL[[length(flagL)+1]] <- data.table(arm=an, build=b, origin=o, flags=nf) }
    # HT monotone fails (count of HT_up=FALSE)
    hm <- sum(grepl("HT_up=FALSE", L))
    htm[[length(htm)+1]] <- data.table(arm=an, build=b, origin=o, ht_up_fail=hm)
  }
}
R <- rbindlist(rein); E <- rbindlist(eich); F <- rbindlist(flagL); H <- rbindlist(htm)
fwrite(R, file.path(DIR,"SUMMARY_reineke.csv"))
fwrite(E, file.path(DIR,"SUMMARY_eichhorn.csv"))
fwrite(merge(F,H,by=c("arm","build","origin")), file.path(DIR,"SUMMARY_flags.csv"))

cat("\n===== REINEKE slope by arm x origin x site =====\n")
W <- dcast(R, arm+origin~site, value.var="slope")
print(W)
cat("\n===== REINEKE OK count (of 4) by arm x origin =====\n")
print(R[, .(ok=sum(reineke_ok), n=.N), by=.(arm,origin)])
cat("\n===== EICHHORN (real/derived HT) =====\n")
print(E)
cat("\n===== FLAGS + HT_up fails =====\n")
print(merge(F,H,by=c("arm","build","origin"))[order(origin,arm)])
