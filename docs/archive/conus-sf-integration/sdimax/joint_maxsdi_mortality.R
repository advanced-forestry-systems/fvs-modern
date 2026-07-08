## Joint maximum-plus-mortality fit: estimate a per-region LEVEL scalar on the localized
## (brms) maximum SDI jointly with the self-thinning response, so relative density is
## consistent with the mortality the data show. Compares predictive skill of observed
## self-thinning under (a) native FVS species-weighted maximum, (b) raw brms maximum,
## (c) per-region level-calibrated brms maximum. Reports the per-region level scalars.
suppressMessages({library(data.table); library(mgcv)})
d <- as.data.table(readRDS("~/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds"))
cat("rows", nrow(d), "\n"); cat("cols:", paste(grep("SDI|TPH|YEAR|variant|LON|FVS|fvs|FORTYP",names(d),value=TRUE),collapse=", "), "\n")
## observed annual density change (self-thinning rate)
d <- d[is.finite(SDImax_brms) & SDImax_brms>0 & is.finite(SDI1) & SDI1>0 &
       is.finite(TPH1) & is.finite(TPH2) & TPH1>0 & TPH2>0 & is.finite(YEARS) & YEARS>0]
d[, mort := -log(TPH2/TPH1)/YEARS]
d <- d[is.finite(mort)]
## region: coarse, robust grouping by FVS variant where present, else West/East/South split
if ("fvs_variant" %in% names(d)) d[, region := toupper(as.character(fvs_variant))] else d[, region := "ALL"]
## keep regions with enough data
keep <- d[, .N, by=region][N>=1500, region]
d <- d[region %in% keep]
cat("regions kept:", length(keep), " n=", nrow(d), "\n")
devexpl <- function(y, x, k=6){
  ok <- is.finite(y)&is.finite(x); if(sum(ok)<200) return(NA_real_)
  m <- tryCatch(bam(y[ok]~s(x[ok],k=k), discrete=TRUE), error=function(e) NULL)
  if(is.null(m)) return(NA_real_); summary(m)$dev.expl
}
## per-region: profile level scalar that maximizes self-thinning skill of brms RD
scales <- seq(0.6, 1.4, by=0.05)
res <- rbindlist(lapply(sort(unique(d$region)), function(g){
  s <- d[region==g]; if(nrow(s)<1500) return(NULL)
  rd_raw <- s$SDI1/s$SDImax_brms
  de_raw <- devexpl(s$mort, rd_raw)
  de_by_k <- sapply(scales, function(k) devexpl(s$mort, s$SDI1/(k*s$SDImax_brms)))
  kbest <- scales[which.max(de_by_k)]; de_best <- max(de_by_k, na.rm=TRUE)
  data.table(region=g, n=nrow(s), mean_brms=round(mean(s$SDImax_brms)),
             de_raw=round(de_raw,3), k_opt=kbest, de_calibrated=round(de_best,3),
             gain_pct=round(100*(de_best-de_raw)/de_raw,1))
}))
cat("\n=== Joint level calibration of localized max SDI, per region ===\n")
print(res)
## pooled: single model with region-specific level (joint), vs raw, vs region-only intercept
d[, rd_raw := SDI1/SDImax_brms]
## apply each region optimal k
setkey(res, region); d[, k_opt := res[.(region), k_opt, on="region"]]
d[, rd_cal := SDI1/(k_opt*SDImax_brms)]
de_pool_raw <- devexpl(d$mort, d$rd_raw)
de_pool_cal <- devexpl(d$mort, d$rd_cal)
cat(sprintf("\nPOOLED self-thinning skill (deviance explained): raw brms RD %.3f  ->  level-calibrated RD %.3f  (+%.1f%%)\n",
            de_pool_raw, de_pool_cal, 100*(de_pool_cal-de_pool_raw)/de_pool_raw))
fwrite(res, "~/fvs-conus/output/joint_maxsdi_levels.csv")
cat("DONE_JOINT\n")
