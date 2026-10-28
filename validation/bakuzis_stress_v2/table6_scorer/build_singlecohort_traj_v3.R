#!/usr/bin/env Rscript
# build_singlecohort_traj_v2.R
# DEFINITIVE single-cohort (even-aged) Bakuzis trajectory builder for all model builds.
#
# Fixes vs prior multi-cohort builder (build_bakuzis_traj.R):
#   The prior builder pooled stands of ALL starting ages into shared age bins
#   (age = AGE0 + PROJ_YEAR over the full AGE0 distribution), so a single age bin
#   mixed "young stand, late cycle" with "old stand, early cycle" -> trajectories
#   looped and produced spurious HT/VOL inversions and chaotic Reineke slopes.
#
# Single-cohort fix (the whole point):
#   Seed ONLY young, even-aged stands with AGE0 in a TIGHT band [AGE_LO, AGE_HI],
#   then follow those SAME stands forward by PROJ_YEAR. age = AGE0 + PROJ_YEAR is
#   then a true cohort age from a common young origin, so trajectories do not loop.
#
# HT basis (fair + documented):
#   - REAL HT (TPA-weighted mean HT_M from treelists) where treelists exist:
#       conus_b2, conus_b1, greg, gompit
#   - DERIVED HT = a*QMD^b (cached fit) where NO treelists exist:
#       engine_default, engine_calibrated, sdicon_b2, sdicon_b1
#
# Origins:
#   - "all"  : global SITE_INDEX quartiles (pooled CONUS). Used for the fully
#              comparable, primary tests (Reineke self-thinning slope).
#              NOTE: pooled CONUS site index is NOT comparable across species, so
#              site-ORDERING on "all" is NOT interpretable -> reported but flagged.
#   - region : NE, SN, PN within-variant SITE_INDEX quartiles (single species pool),
#              so site ordering IS interpretable. Used for the site-ordering test.
#
# Usage: Rscript build_singlecohort_traj_v2.R <build_label> [AGE_LO AGE_HI]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
BUILD  <- a[1]
AGE_LO <- if (length(a) >= 2) as.numeric(a[2]) else 10
AGE_HI <- if (length(a) >= 3) as.numeric(a[3]) else 40

SCR   <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
EQ    <- file.path(SCR, "conus_eq_proj")
SIDIR <- file.path(SCR, "standinit_by_variant")
OUT   <- file.path(EQ, "bakuzis_singlecohort"); dir.create(OUT, showWarnings = FALSE)
FF <- 0.40; FT2AC <- 0.2296; IN2CM <- 2.54
REGIONS <- c("NE","SN","PN")   # within-region origins for site ordering

## ---- cached HT~QMD fit for derived-HT builds (built by prior run) ----
fitf <- file.path(EQ, "bakuzis_redo", "htqmd_fit.rds")
if (!file.exists(fitf)) stop("missing cached htqmd_fit.rds")
ff_ <- readRDS(fitf); AHT <- ff_$a; BHT <- ff_$b
ht_qmd <- function(QMD) AHT * QMD^BHT
cat(sprintf("derived HT = %.3f * QMD^%.3f (cached, n=%d R2=%.3f)\n", AHT, BHT, ff_$n, ff_$r2))

## ---- standinit AGE + SITE_INDEX (per-stand), keyed by STAND_CN ----
SI <- rbindlist(lapply(list.files(SIDIR, pattern="^standinit_[A-Z]+\\.csv$", full.names=TRUE),
  function(f) tryCatch(fread(f, colClasses=list(character="STAND_CN"),
    select=c("STAND_CN","AGE","SITE_INDEX")), error=function(e) NULL)), fill=TRUE)
SI <- SI[!is.na(STAND_CN)]
SI[, STAND_CN := sub("\\..*$","",STAND_CN)]
SI[, AGE := suppressWarnings(as.numeric(AGE))]
SI[, SITE_INDEX := suppressWarnings(as.numeric(SITE_INDEX))]
SI <- unique(SI, by="STAND_CN")
# global quartile boundaries from full valid SI distribution
qb_g <- quantile(SI$SITE_INDEX[is.finite(SI$SITE_INDEX) & SI$SITE_INDEX>0],
                 c(0,.25,.5,.75,1), na.rm=TRUE)
cat("global SITE_INDEX quartile boundaries (ft):", paste(round(qb_g,1),collapse=", "), "\n")
# single-cohort seed mask
SI[, seed := is.finite(AGE) & AGE >= AGE_LO & AGE <= AGE_HI]
cat(sprintf("single-cohort seed: %d stands with AGE0 in [%g,%g]\n",
            sum(SI$seed, na.rm=TRUE), AGE_LO, AGE_HI))
setkey(SI, STAND_CN)

## ---- treelist -> TPA-weighted mean HT_M per (stand,config,proj_year) ----
## SEED_SET is set per-build (only single-cohort seed stands) so the accumulated
## HT table stays tiny and we never hold the full 40+GB treelists in memory.
SEED_SET <- character(0)
tl_meanht <- function(tf){
  if(!file.exists(tf)) return(NULL)
  t <- tryCatch(fread(tf, colClasses=list(character="STAND_CN"),
        select=c("STAND_CN","CONFIG","PROJ_YEAR","HT_M","TPA")), error=function(e) NULL)
  if(is.null(t)||!nrow(t)) return(NULL)
  t[, STAND_CN := sub("\\..*$","",STAND_CN)]
  if(length(SEED_SET)) t <- t[STAND_CN %chin% SEED_SET]   # keep only seed stands
  if(!nrow(t)) return(NULL)
  suppressWarnings({ t[, HT_M := as.numeric(HT_M)]; t[, TPA := as.numeric(TPA)] })
  t <- t[is.finite(HT_M) & is.finite(TPA) & TPA>0]
  if(!nrow(t)) return(NULL)
  res <- t[, .(HT = sum(HT_M*TPA)/sum(TPA)), by=.(STAND_CN,CONFIG,PROJ_YEAR)]
  rm(t); gc(FALSE); res
}
tl_meanht_many <- function(files){
  out <- vector("list", length(files))
  for(i in seq_along(files)){ out[[i]] <- tl_meanht(files[i]); if(i %% 5 == 0) gc(FALSE) }
  rbindlist(out, fill=TRUE)
}

## ---- build registry: metrics glob, treelist glob (NULL = derived HT), config filter ----
eqd  <- file.path(EQ,"out_conus_eq")
sdd  <- file.path(EQ,"out_conus_eq_sdiconstrained")
gd   <- file.path(EQ,"out_conus_eq_greg")
gpd  <- file.path(EQ,"out_conus_eq_gompit")
engd <- file.path(SCR,"out_conus_engine_v4")
v4d  <- file.path(EQ,"out_v4_seed")
enghtd <- file.path(SCR,"out_conus_engine_seed_ht")
lf <- function(dir, pat) list.files(dir, pattern=pat, full.names=TRUE)

REG <- list(
  engine_default    = list(metrics=lf(engd,"^conus_[a-z]+_b[0-9]+\\.csv$"), tl=NULL, cfg="default"),
  engine_calibrated = list(metrics=lf(engd,"^conus_[a-z]+_b[0-9]+\\.csv$"), tl=NULL, cfg="calibrated"),
  conus_b2          = list(metrics=lf(eqd,"_conus_b2_metrics\\.csv$"),  tl=lf(eqd,"_conus_b2_treelists\\.csv$"),  cfg=NULL),
  conus_b1          = list(metrics=lf(eqd,"_conus_b1_metrics\\.csv$"),  tl=lf(eqd,"_conus_b1_treelists\\.csv$"),  cfg=NULL),
  greg              = list(metrics=lf(gd,"_greg_metrics\\.csv$"),       tl=lf(gd,"_greg_treelists\\.csv$"),       cfg=NULL),
  gompit            = list(metrics=lf(gpd,"_gompit_metrics\\.csv$"),    tl=lf(gpd,"_gompit_treelists\\.csv$"),    cfg=NULL),
  sdicon_b2         = list(metrics=lf(sdd,"_conus_b2_metrics\\.csv$"),  tl=NULL, cfg=NULL),
  sdicon_b1         = list(metrics=lf(sdd,"_conus_b1_metrics\\.csv$"),  tl=NULL, cfg=NULL),
  v4_b2             = list(metrics=lf(v4d,"_conus_b2_metrics\\.csv$"), tl=NULL, cfg=NULL, ht_col="HT_M_MEAN"),
  engine_default_ht    = list(metrics=lf(enghtd,"^conus_[a-z]+_b[0-9]+\\.csv$"), tl=NULL, cfg="default",    ht_col="HT_M_MEAN"),
  engine_calibrated_ht = list(metrics=lf(enghtd,"^conus_[a-z]+_b[0-9]+\\.csv$"), tl=NULL, cfg="calibrated", ht_col="HT_M_MEAN")
)
if(!BUILD %in% names(REG)) stop("unknown build: ", BUILD, " | choose: ", paste(names(REG),collapse=", "))
spec <- REG[[BUILD]]
cat(sprintf("\n=== BUILD %s : %d metrics files, treelists=%s, HT=%s ===\n",
    BUILD, length(spec$metrics), if(is.null(spec$tl))"NONE" else length(spec$tl),
    if(!is.null(spec$ht_col)) paste0("REAL metrics col ",spec$ht_col) else if(is.null(spec$tl))"DERIVED f(QMD)" else "REAL treelist"))

## ---- load metrics, restrict to single-cohort seed, assign age = AGE0 + PROJ_YEAR ----
M <- rbindlist(lapply(spec$metrics, function(f)
  tryCatch(fread(f, colClasses=list(character="STAND_CN")), error=function(e) NULL)), fill=TRUE)
if(!nrow(M)) stop("no metrics for ", BUILD)
M[, STAND_CN := sub("\\..*$","",STAND_CN)]
if(!is.null(spec$cfg)) M <- M[CONFIG == spec$cfg]
M[, BAPH := BA_FT2AC*FT2AC][, QMD := QMD_IN*IN2CM][, TPH := as.numeric(TPH)]
M <- M[is.finite(BAPH)&is.finite(QMD)&is.finite(TPH)&TPH>0&QMD>0]
M <- merge(M, SI[,.(STAND_CN, AGE0=AGE, SI=SITE_INDEX, seed)], by="STAND_CN", all.x=TRUE)
M <- M[seed == TRUE]                              # SINGLE-COHORT: only young even-aged seeds
M[, PROJ_YEAR := as.integer(PROJ_YEAR)]
M[, age := AGE0 + PROJ_YEAR]                      # true cohort age from a common young origin
M <- M[is.finite(age)]
cat(sprintf("rows after single-cohort seed restriction: %d (stands=%d)\n", nrow(M), uniqueN(M$STAND_CN)))

## ---- HT: real treelist (TPA-wtd) or derived ----
if(!is.null(spec$ht_col)){
  hc <- spec$ht_col
  if(!hc %in% names(M)) stop("ht_col ", hc, " not in metrics for ", BUILD)
  M[, HT := suppressWarnings(as.numeric(get(hc)))]
  M <- M[is.finite(HT) & HT>0]
  ht_src <- paste0("real metrics-column HT (", hc, ", meters)")
} else if(!is.null(spec$tl)){
  SEED_SET <- unique(M$STAND_CN)                  # restrict treelist reads to seed stands only
  cat(sprintf("treelist read restricted to %d seed stands\n", length(SEED_SET)))
  HT <- tl_meanht_many(spec$tl)
  HT[, PROJ_YEAR := as.integer(PROJ_YEAR)]
  M <- merge(M, HT, by=c("STAND_CN","CONFIG","PROJ_YEAR"), all.x=TRUE)
  M <- M[is.finite(HT)]
  ht_src <- "real treelist TPA-wtd mean HT_M"
} else {
  M[, HT := ht_qmd(QMD)]
  ht_src <- "derived HT=a*QMD^b"
}
M[, VOL := BAPH*HT*FF]; M[, agebin := round(age/5)*5]

## ---- writer: aggregate one origin to site x agebin ----
write_traj <- function(MM, origin_lab, qb_use, fn_suffix){
  MM <- MM[is.finite(SI) & SI>0]
  if(!nrow(MM)) return(invisible())
  MM[, siteclass := cut(SI, breaks=qb_use, include.lowest=TRUE, labels=1:4)]
  MM <- MM[!is.na(siteclass)]
  agg <- MM[, .(HT=mean(HT), TPH=mean(TPH), BAPH=mean(BAPH), QMD=mean(QMD),
                VOL=mean(VOL), n=.N), by=.(siteclass, agebin)]
  setnames(agg, c("siteclass","agebin"), c("site","age"))
  agg <- agg[n>=5][order(site,age)]; agg[, origin := origin_lab]
  if(!nrow(agg)) return(invisible())
  of <- file.path(OUT, paste0(BUILD, fn_suffix, ".csv"))
  fwrite(agg[,.(origin,site,age,HT,TPH,BAPH,QMD,VOL,n)], of)
  # QMD-monotone sanity per (site)
  mono <- agg[, .(qmd_mono = all(diff(QMD) >= -1e-6)), by=site]
  cat(sprintf("  %-8s rows=%4d  QMD-monotone-by-site: %s  -> %s\n",
      origin_lab, nrow(agg), paste(mono$qmd_mono, collapse=","), basename(of)))
  invisible(agg)
}

cat("HT source:", ht_src, "\n")
## (1) "all" origin, global quartiles -> Reineke / self-thinning (primary, comparable)
write_traj(copy(M), "all", qb_g, "_sc_all")

## (2) within-region origins -> site ordering (interpretable, single species pool)
if("VARIANT" %in% names(M)){
  for(rg in REGIONS){
    Mr <- M[toupper(VARIANT) == rg]
    if(!nrow(Mr)) { cat(sprintf("  region %s: no rows\n", rg)); next }
    siv <- Mr$SI[is.finite(Mr$SI) & Mr$SI>0]
    if(length(siv) < 20){ cat(sprintf("  region %s: too few SI\n", rg)); next }
    qb_r <- quantile(siv, c(0,.25,.5,.75,1), na.rm=TRUE)
    write_traj(copy(Mr), rg, qb_r, paste0("_sc_", tolower(rg)))
  }
}
cat("DONE", BUILD, "\n")
