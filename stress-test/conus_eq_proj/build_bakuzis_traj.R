#!/usr/bin/env Rscript
# Build per-build Bakuzis trajectory CSVs (long: site-class x age) for the six-arm comparison.
# Site class = quartile of standinit SITE_INDEX (variant-native productivity class, per-stand, complete).
# Age = standinit AGE + PROJ_YEAR. HT: from treelists where available, else from a
# fitted HT~QMD curve (noted as approximation). Metric units. VOL = BAPH*HT*ff.
suppressPackageStartupMessages(library(data.table))
SCR <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
EQ  <- file.path(SCR,"conus_eq_proj")
OUT <- file.path(EQ,"bakuzis_redo"); dir.create(OUT, showWarnings=FALSE)
SIDIR <- file.path(SCR,"standinit_by_variant")
FF <- 0.40
FT2AC_to_M2HA <- 0.2296
IN_to_CM <- 2.54

## ---- site index / age from standinit (per-stand, complete), keyed by STAND_CN ----
si_files <- list.files(SIDIR, pattern="^standinit_[A-Z]+\\.csv$", full.names=TRUE)
si_list <- lapply(si_files, function(f){
  x <- tryCatch(fread(f, colClasses=list(character="STAND_CN"),
                      select=c("STAND_CN","AGE","SITE_INDEX")), error=function(e) NULL)
  x
})
SI <- rbindlist(si_list, fill=TRUE)
SI <- SI[!is.na(STAND_CN)]
SI[, AGE := suppressWarnings(as.numeric(AGE))]
SI[, SITE_INDEX := suppressWarnings(as.numeric(SITE_INDEX))]
SI <- unique(SI, by="STAND_CN")
cat(sprintf("standinit rows: %d ; SITE_INDEX nonNA: %d ; AGE nonNA: %d\n",
            nrow(SI), sum(!is.na(SI$SITE_INDEX)), sum(!is.na(SI$AGE))))

sivals <- SI$SITE_INDEX[is.finite(SI$SITE_INDEX) & SI$SITE_INDEX>0]
qb <- quantile(sivals, probs=c(0,.25,.5,.75,1), na.rm=TRUE)
cat("SITE_INDEX quartile boundaries (ft):", paste(round(qb,1),collapse=", "), "\n")
SI[, siteclass := cut(SITE_INDEX, breaks=qb, include.lowest=TRUE, labels=1:4)]
setkey(SI, STAND_CN)

tl_meanht <- function(tf){
  if(!file.exists(tf)) return(NULL)
  t <- tryCatch(fread(tf, colClasses=list(character="STAND_CN"),
                      select=c("STAND_CN","CONFIG","PROJ_YEAR","HT_M","TPA")), error=function(e) NULL)
  if(is.null(t)||!nrow(t)) return(NULL)
  # robust numeric coercion (some files have stray non-numeric tokens in HT_M/TPA)
  suppressWarnings({ t[, HT_M := as.numeric(HT_M)]; t[, TPA := as.numeric(TPA)] })
  t <- t[is.finite(HT_M) & is.finite(TPA) & TPA>0]
  if(!nrow(t)) return(NULL)
  res <- t[, .(HT = sum(HT_M*TPA)/sum(TPA)), by=.(STAND_CN,CONFIG,PROJ_YEAR)]
  rm(t); res
}
# memory-friendly: build per-stand mean HT by streaming files one at a time
tl_meanht_many <- function(files){
  out <- vector("list", length(files))
  for(i in seq_along(files)){ out[[i]] <- tl_meanht(files[i]); if(i %% 5 == 0) gc(FALSE) }
  rbindlist(out, fill=TRUE)
}

assemble <- function(metrics_files, treelist_files=NULL, label, config_filter=NULL,
                     ht_from_qmd=NULL){
  M <- rbindlist(lapply(metrics_files, function(f)
    tryCatch(fread(f, colClasses=list(character="STAND_CN")), error=function(e) NULL)), fill=TRUE)
  if(!nrow(M)){ cat("  NO metrics for",label,"\n"); return(NULL) }
  if(!is.null(config_filter)) M <- M[CONFIG %in% config_filter]
  M[, BAPH := BA_FT2AC * FT2AC_to_M2HA]
  M[, QMD  := QMD_IN * IN_to_CM]
  M[, TPH  := as.numeric(TPH)]
  M <- M[is.finite(BAPH)&is.finite(QMD)&is.finite(TPH)&TPH>0&QMD>0]
  M <- merge(M, SI[,.(STAND_CN,AGE0=AGE,siteclass)], by="STAND_CN", all.x=TRUE)
  M[, age := AGE0 + as.numeric(PROJ_YEAR)]
  M <- M[is.finite(age) & !is.na(siteclass)]
  if(!is.null(treelist_files)){
    HT <- tl_meanht_many(treelist_files)
    if(!is.null(HT) && nrow(HT)){
      M[, PROJ_YEAR := as.integer(PROJ_YEAR)]; HT[, PROJ_YEAR := as.integer(PROJ_YEAR)]
      M <- merge(M, HT, by=c("STAND_CN","CONFIG","PROJ_YEAR"), all.x=TRUE)
    } else M[, HT := NA_real_]
    ht_src <- "treelist (TPA-wtd mean HT_M)"
  } else if(!is.null(ht_from_qmd)){
    M[, HT := ht_from_qmd(QMD)]
    ht_src <- "derived HT=f(QMD) [approx; no treelists]"
  } else { M[, HT := NA_real_]; ht_src <- "none" }
  M <- M[is.finite(HT)]
  if(!nrow(M)){ cat("  NO usable rows for",label,"\n"); return(NULL) }
  M[, VOL := BAPH*HT*FF]
  M[, agebin := round(age/5)*5]
  agg <- M[, .(HT=mean(HT), TPH=mean(TPH), BAPH=mean(BAPH), QMD=mean(QMD),
               VOL=mean(VOL), n=.N), by=.(siteclass, agebin)]
  setnames(agg, c("siteclass","agebin"), c("site","age"))
  agg <- agg[n>=5]
  agg <- agg[order(site,age)]
  agg[, origin := "all"]
  of <- file.path(OUT, paste0(label,"_traj.csv"))
  fwrite(agg[,.(origin,site,age,HT,TPH,BAPH,QMD,VOL,n)], of)
  cat(sprintf("  %-28s rows=%4d  HT_src=%s  -> %s\n", label, nrow(agg), ht_src, basename(of)))
  agg
}

cat("\n=== STEP A: HT~QMD fit (reuse cached if present) ===\n")
eq_dir <- file.path(EQ,"out_conus_eq")
fitf <- file.path(OUT,"htqmd_fit.rds")
if(file.exists(fitf)){
  ff_ <- readRDS(fitf); a <- ff_$a; b <- ff_$b
  cat(sprintf("  reused cached fit: HT = %.3f * QMD^%.3f (n=%d R2=%.3f)\n", a,b,ff_$n,ff_$r2))
} else {
  eq_b2_tl_fit <- list.files(eq_dir, pattern="_conus_b2_treelists\\.csv$", full.names=TRUE)
  eq_b2_m_fit  <- list.files(eq_dir, pattern="_conus_b2_metrics\\.csv$", full.names=TRUE)
  HTfit <- tl_meanht_many(eq_b2_tl_fit)
  Mfit  <- rbindlist(lapply(eq_b2_m_fit, function(f)
    fread(f, colClasses=list(character="STAND_CN"), select=c("STAND_CN","CONFIG","PROJ_YEAR","QMD_IN"))), fill=TRUE)
  Mfit[, QMD := QMD_IN*IN_to_CM][, PROJ_YEAR:=as.integer(PROJ_YEAR)]
  HTfit[, PROJ_YEAR:=as.integer(PROJ_YEAR)]
  FTm <- merge(Mfit, HTfit, by=c("STAND_CN","CONFIG","PROJ_YEAR"))
  FTm <- FTm[is.finite(QMD)&is.finite(HT)&QMD>0&HT>0]
  hdfit <- lm(log(HT) ~ log(QMD), FTm)
  a <- exp(coef(hdfit)[1]); b <- coef(hdfit)[2]
  cat(sprintf("  HT-QMD fit (n=%d): HT = %.3f * QMD^%.3f  (R2=%.3f)\n",
              nrow(FTm), a, b, summary(hdfit)$r.squared))
  saveRDS(list(a=a,b=b,n=nrow(FTm),r2=summary(hdfit)$r.squared), fitf)
  rm(HTfit,Mfit,FTm); gc(FALSE)
}
ht_qmd <- function(QMD) a*QMD^b

cat("\n=== STEP B: assemble each build (skip already-written) ===\n")
done <- function(lab) file.exists(file.path(OUT, paste0(lab,"_traj.csv")))
eng_b0 <- list.files(file.path(SCR,"out_conus_engine_v4"), pattern="_b0\\.csv$", full.names=TRUE)
if(!done("engine_default"))    assemble(eng_b0, NULL, "engine_default",    config_filter="default",    ht_from_qmd=ht_qmd)
if(!done("engine_calibrated")) assemble(eng_b0, NULL, "engine_calibrated", config_filter="calibrated", ht_from_qmd=ht_qmd)
if(!done("conus_b2")){
  eq_b2_metrics <- list.files(eq_dir, pattern="_conus_b2_metrics\\.csv$", full.names=TRUE)
  eq_b2_tl      <- list.files(eq_dir, pattern="_conus_b2_treelists\\.csv$", full.names=TRUE)
  assemble(eq_b2_metrics, eq_b2_tl, "conus_b2"); gc(FALSE)
}
if(!done("conus_b1")){
  eq_b1_metrics <- list.files(eq_dir, pattern="_conus_b1_metrics\\.csv$", full.names=TRUE)
  eq_b1_tl      <- list.files(eq_dir, pattern="_conus_b1_treelists\\.csv$", full.names=TRUE)
  assemble(eq_b1_metrics, eq_b1_tl, "conus_b1"); gc(FALSE)
}
gd <- file.path(EQ,"out_conus_eq_greg")
if(!done("greg")){ assemble(list.files(gd,pattern="_greg_metrics\\.csv$",full.names=TRUE),
         list.files(gd,pattern="_greg_treelists\\.csv$",full.names=TRUE), "greg"); gc(FALSE) }
gpd <- file.path(EQ,"out_conus_eq_gompit")
if(!done("gompit")){ assemble(list.files(gpd,pattern="_gompit_metrics\\.csv$",full.names=TRUE),
         list.files(gpd,pattern="_gompit_treelists\\.csv$",full.names=TRUE), "gompit"); gc(FALSE) }
sd <- file.path(EQ,"out_conus_eq_sdiconstrained")
if(!done("sdicon_b2")) assemble(list.files(sd,pattern="_conus_b2_metrics\\.csv$",full.names=TRUE),
         NULL, "sdicon_b2", ht_from_qmd=ht_qmd)
if(!done("sdicon_b1")) assemble(list.files(sd,pattern="_conus_b1_metrics\\.csv$",full.names=TRUE),
         NULL, "sdicon_b1", ht_from_qmd=ht_qmd)

cat("\nDONE. Trajectory CSVs in", OUT, "\n")
saveRDS(qb, file.path(OUT,"siteclass_boundaries.rds"))
