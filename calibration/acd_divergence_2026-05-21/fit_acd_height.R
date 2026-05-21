## fit_acd_height.R
##
## Robust annual height multipliers (dHt.mult) for the Acadian model, fit from
## MEASURED HT_t1/HT_t2. The measured height growth is very noisy (about 23
## percent of records show negative growth from top breakage and inconsistent
## height protocols between visits), so an SSE fit saturates. Instead we use a
## ratio of trimmed means per species, which is insensitive to symmetric noise:
##   dHt.mult = trimmedMean(observed height growth) / trimmedMean(predicted)
## Predicted height growth comes from the real Acadian dHT_fun iterated annually.
##
## Merges dHt.mult into the existing acd_annual_calibration.csv.
## Run: module load gcc/12.3.0 R/4.4.0; Rscript fit_acd_height.R
suppressWarnings(suppressMessages({ library(data.table); library(plyr); library(dplyr); library(purrr); library(tibble) }))

ACD <- Sys.getenv("ACD_PATH",
  "~/upstream_fvs_check/ForestVegetationSimulator-Interface/fvsOL/inst/extdata/AcadianGY.R")
source(ACD)
cat("sourced", AcadianVersionTag, "\n")

DG  <- Sys.getenv("ACD_DG", "/users/PUOM0008/crsfaaron/fvs-conus/data/processed/acd/diameter_growth.csv")
CSV <- Sys.getenv("CAL_CSV", "/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/acd_divergence_2026-05-21/acd_annual_calibration.csv")
CAP <- as.integer(Sys.getenv("FIT_CAP", "6000"))
TRIM <- as.numeric(Sys.getenv("TRIM", "0.10"))     # symmetric trim fraction
FT2AC_to_M2HA <- 0.2295684; CCFL_CONST <- 100

xwalk <- c("12"="BF","91"="NS","94"="WS","95"="BS","97"="RS","71"="TA","105"="JP",
  "125"="RN","130"="SC","129"="WP","241"="WC","261"="EH","315"="ST","316"="RM",
  "317"="SV","318"="SM","319"="MM","371"="YB","375"="PB","379"="GB","531"="AB",
  "541"="WA","543"="BA","544"="GA","601"="BN","741"="BP","743"="BT","746"="QA",
  "762"="BC","802"="WO","833"="RO","837"="BO","951"="BW","972"="AE")
si_to_csi <- function(si_ft){ csi<-approxfun(c(0,8,14,20),c(0,8,12,14),rule=2)(si_ft*0.3048); ifelse(is.na(csi)|csi<=0,12,csi) }

cat("reading", DG, "\n")
d <- fread(DG, select=c("PLT_CN","DIA_t1","CR_t1","HT_t1","HT_t2","BAL","SI","years_interval","SPCD","TPA"))
setnames(d, c("DIA_t1","CR_t1","HT_t1","HT_t2","years_interval","SPCD","TPA"),
            c("d1","cr","h1","h2","yi","spcd","tpa"))
d <- d[is.finite(d1)&d1>0 & is.finite(h1)&h1>0 & is.finite(h2) & yi>=1&yi<=15 & is.finite(BAL) & is.finite(tpa)&tpa>0]
d[, `:=`(alpha=xwalk[as.character(spcd)], DBHcm=d1*2.54, isSW=as.integer(spcd<300),
         CSI=si_to_csi(SI), CRp=pmin(0.95,pmax(0.05,cr/100)),
         BALm2ha=BAL*FT2AC_to_M2HA, yi=as.integer(round(yi)),
         HTm=h1*0.3048, obs=(h2-h1)*0.3048)]
d[, w := DBHcm^2 * tpa]
setorder(d, PLT_CN, -DBHcm)
d[, swf := { tot<-cumsum(w)-w; sw<-cumsum(w*isSW)-w*isSW; pf<-sum(w*isSW)/sum(w); fifelse(tot>0, sw/tot, pf) }, by=PLT_CN]
d[, `:=`(BALSW=BALm2ha*swf, BALHW=BALm2ha*(1-swf))]
df <- d[!is.na(alpha)]
cat("records with measured HT:", nrow(df), "\n")

reh <- as.data.frame(dht.fun.spp)
get_reh <- function(a){ r<-reh[reh$Spp==a,]; if(!nrow(r)) c(0,0) else as.numeric(r[1,c("dht.b0.spp","dht.b2.spp")]) }
pred_ht <- function(s,b,c){ HT<-s$HTm; start<-HT; for(yr in seq_len(max(s$yi))){ inc<-c*dHT_fun("X",pmax(HT,1.4),s$CRp,CCFL_CONST,s$CSI,b[1],b[2]); inc[!is.finite(inc)|inc<0]<-0; a<-s$yi>=yr; HT[a]<-HT[a]+inc[a] }; HT-start }
tmean <- function(x) mean(x, trim=TRIM, na.rm=TRUE)

set.seed(1)
species <- names(sort(table(df$alpha), decreasing=TRUE))
res <- data.frame()
for (a in species) {
  sub <- df[alpha==a]; if (nrow(sub)>CAP) sub <- sub[sample(.N,CAP)]
  b <- get_reh(a)
  pred1 <- pred_ht(sub, b, 1)
  ## keep central records by observed growth (symmetric trim removes noise tails)
  qlo <- quantile(sub$obs, TRIM, na.rm=TRUE); qhi <- quantile(sub$obs, 1-TRIM, na.rm=TRUE)
  keep <- sub$obs>=qlo & sub$obs<=qhi & is.finite(pred1) & pred1>0
  cfac <- tmean(sub$obs[keep]) / mean(pred1[keep])
  cfac <- min(max(cfac, 0.3), 2.5)
  res <- rbind(res, data.frame(SP=a, n_ht=nrow(sub),
    dHt.mult=round(cfac,4),
    obs_grow_m_yr=round(tmean(sub$obs[keep])/mean(sub$yi[keep]),3),
    pred_grow_m_yr=round(mean(pred1[keep])/mean(sub$yi[keep]),3),
    frac_neg=round(mean(sub$obs<0),3), stringsAsFactors=FALSE))
  cat(sprintf("  %-3s nH=%5d dHt=%.3f  obs=%.3f pred=%.3f m/yr  fracNeg=%.2f\n",
              a, nrow(sub), cfac, tmean(sub$obs[keep])/mean(sub$yi[keep]),
              mean(pred1[keep])/mean(sub$yi[keep]), mean(sub$obs<0)))
}

## merge dHt.mult into the existing calibration table
cal <- read.csv(CSV, stringsAsFactors=FALSE)
cal$dHt.mult <- NULL
cal <- merge(cal, res[,c("SP","dHt.mult")], by="SP", all.x=TRUE)
cal$dHt.mult[is.na(cal$dHt.mult)] <- 1
## restore a sensible column order
front <- c("SP","dDBH.mult","dHt.mult","mort.mult")
cal <- cal[, c(front, setdiff(names(cal), front))]
cal <- cal[order(-cal$n), ]
write.csv(cal, CSV, row.names=FALSE)
write.csv(res, sub(".csv$","_height_fit.csv",CSV), row.names=FALSE)

cat(sprintf("\nupdated %s with dHt.mult for %d species\n", CSV, nrow(res)))
cat("dHt.mult range:", round(range(res$dHt.mult),3),
    " median:", round(median(res$dHt.mult),3), "\n")
cat("median fraction negative-growth records:", round(median(res$frac_neg),3),
    " (data-quality flag)\n")
cat("DONE\n")
