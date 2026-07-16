#!/usr/bin/env Rscript
# hg_climate_vs_none.R
# Does climate help HEIGHT GROWTH? P11 crosswalk now supplies TD/EMT, so Greg's
# est_hg climate terms can be tested on our pairs. Single-step annual form of
# est_hg, max-height fixed per species at Greg's B0, starts from Greg's B1..B8:
#   hg = B0 * B1*B2*cr^B3 * exp(-B1*ht - B4*ccfl - B8*sqrt(cch) - B5*elev
#                               [+ B6*sqrt(td) + B7*emt]) * (1-exp(-B1*ht))^(B2-1)
#   climate: with B6*sqrt(td)+B7*emt ; none: drop those two terms (refit rest).
# Fragile 8-param nls; per-species Greg starts + fixed max-height stabilize it.
# Reports convergence rate + held-out RMSE climate vs none. Author: Cowork autopilot 2026-07-03.

suppressPackageStartupMessages({ library(data.table); library(minpack.lm) })
set.seed(20260703)
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/hg_climate"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rmse <- function(a,b) sqrt(mean((a-b)^2)); bias <- function(a,b) mean(a-b)

hp <- as.data.table(readRDS("~/fvs_remodeling/rds/hg_parms.RDS"))   # spcd,B0..B8
setnames(hp, "spcd", "SPCD")
xw <- fread("/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/emt_crosswalk/plot_key_emt_td_crosswalk.csv")
d <- as.data.table(readRDS("~/fvs-conus/calibration/data/conus_remeasurement_pairs.rds"))
d <- merge(d, xw[, .(plot_key, EMT, TD)], by="plot_key", all.x=TRUE)
d <- d[is.finite(HT1)&is.finite(HT2)&is.finite(CR1)&is.finite(CCFL1)&is.finite(CCH1)&
       is.finite(ELEV)&is.finite(EMT)&is.finite(TD)&is.finite(YEARS)&YEARS>0 &
       CR1>0 & CR1<1 & HT1>=4.5 & STATUS1==1 & STATUS2==1]
d[, hg_obs := (HT2-HT1)/YEARS]
d <- d[hg_obs > -0.5 & hg_obs < 5]
setnames(d, c("CR1","HT1","CCFL1","CCH1","ELEV"), c("cr","ht","ccfl","cch","elev"))

pred_hg <- function(p, mh, cr, ht, ccfl, cch, elev, td, emt, climate) {
  e <- -p["B1"]*ht - p["B4"]*ccfl - p["B8"]*sqrt(cch) - p["B5"]*elev
  if (climate) e <- e + p["B6"]*sqrt(td) + p["B7"]*emt
  mh * p["B1"]*p["B2"]*cr^p["B3"] * exp(e) * (1 - exp(-p["B1"]*ht))^(p["B2"]-1)
}
fit_hg <- function(dat, mh, st, climate) {
  if (climate)
    frm <- hg_obs ~ mh*B1*B2*cr^B3*exp(-B1*ht - B4*ccfl - B8*sqrt(cch) - B5*elev + B6*sqrt(td) + B7*emt)*(1-exp(-B1*ht))^(B2-1)
  else
    frm <- hg_obs ~ mh*B1*B2*cr^B3*exp(-B1*ht - B4*ccfl - B8*sqrt(cch) - B5*elev)*(1-exp(-B1*ht))^(B2-1)
  keep <- if (climate) c("B1","B2","B3","B4","B5","B6","B7","B8") else c("B1","B2","B3","B4","B5","B8")
  tryCatch(nlsLM(frm, data=cbind(dat, mh=mh), start=as.list(st[keep]),
                 control=nls.lm.control(maxiter=300)), error=function(e) NULL)
}

sp_n <- d[, .N, by=SPCD][N >= 3000]
teC <- list(); teN <- list(); conv <- data.table(SPCD=integer(), climate=logical())
for (s in intersect(sp_n$SPCD, hp$SPCD)) {
  ds <- d[SPCD==s]; if (nrow(ds) > 15000) ds <- ds[sample(.N, 15000)]
  mh <- hp[SPCD==s, B0][1]; st <- unlist(hp[SPCD==s, .(B1,B2,B3,B4,B5,B6,B7,B8)])
  if (!is.finite(mh) || any(!is.finite(st))) next
  i <- sample(nrow(ds), floor(nrow(ds)*0.7)); tr <- ds[i]; te <- ds[-i]
  for (cl in c(TRUE, FALSE)) {
    m <- fit_hg(tr, mh, st, cl); if (is.null(m)) next
    conv <- rbind(conv, data.table(SPCD=s, climate=cl))
    pr <- pred_hg(coef(m), mh, te$cr, te$ht, te$ccfl, te$cch, te$elev, te$td, te$emt, cl)
    pr <- pmin(pmax(pr, 0), 5)
    tab <- data.table(obs=te$hg_obs, pred=pr)
    if (cl) teC[[as.character(s)]] <- tab else teN[[as.character(s)]] <- tab
  }
}
TC <- rbindlist(teC); TN <- rbindlist(teN)
cat(sprintf("species attempted: %d | climate converged: %d | none converged: %d\n",
    length(intersect(sp_n$SPCD, hp$SPCD)), uniqueN(conv[climate==TRUE]$SPCD), uniqueN(conv[climate==FALSE]$SPCD)))
mets <- function(T) if (is.null(T) || !nrow(T)) list(n=0L,RMSE=NA_real_,bias=NA_real_,r=NA_real_) else
  list(n=nrow(T), RMSE=rmse(T$pred,T$obs), bias=bias(T$pred,T$obs), r=suppressWarnings(cor(T$pred,T$obs)))
mc <- mets(TC); mn <- mets(TN)
sc <- data.table(
  version = c("climate (td+emt)","none"),
  n_test  = c(mc$n, mn$n),
  RMSE    = c(mc$RMSE, mn$RMSE),
  bias    = c(mc$bias, mn$bias),
  r       = c(mc$r, mn$r))
cat("\n=== HG climate vs none (held-out, ft/yr; single-step annual) ===\n"); print(sc)
fwrite(sc, file.path(OUT, "hg_climate_vs_none_scorecard.csv"))
cat("\nMean observed HG:", round(mean(TN$obs),3), "ft/yr\n")
cat("Benchmark context: this is height growth (ft/yr), not DG; compare climate vs none convergence + RMSE.\n")
