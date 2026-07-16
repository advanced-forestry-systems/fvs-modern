##=============================================================================
## 50hg_extract_hg_v8rd_residuals.R
## Modifier-residual extractor for HG v8rd (v7 + RD). Mirrors the v5 extractor
## but reconstructs the v8rd eta exactly (adds z_sp, a2_quad*ln_ht^2, and the
## RD interactions a_bard*BA*RD + a_blrd*BAL_raw*RD). Requires the FULL fit
## object (re-run without --compact). Hard-gates on |mean(residual)| to catch
## any eta-reconstruction error before downstream modifier fits run.
##=============================================================================
suppressPackageStartupMessages({ library(data.table); library(cmdstanr) })
`%||%` <- function(a,b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(n,d=NULL){ m<-grep(paste0("^--",n,"="),args,value=TRUE); if(!length(m)) d else sub(paste0("^--",n,"="),"",m[1]) }
FIT_FILE  <- get_arg("fit")
META_FILE <- get_arg("meta", sub("_fit\\.rds$","_meta.rds", FIT_FILE %||% ""))
PAIRS_FILE<- get_arg("pairs","calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_FILE  <- get_arg("out")
GATE      <- as.numeric(get_arg("gate","0.05"))
stopifnot(!is.null(FIT_FILE), file.exists(FIT_FILE), file.exists(META_FILE), !is.null(OUT_FILE))
cat("== 50hg_extract_hg_v8rd_residuals.R ==\n  fit:",FIT_FILE,"\n  out:",OUT_FILE,"\n\n")
dir.create(dirname(OUT_FILE), recursive=TRUE, showWarnings=FALSE)
fit<-readRDS(FIT_FILE); meta<-readRDS(META_FILE); pairs<-as.data.table(readRDS(PAIRS_FILE))
sp_levels<-meta$sp_levels; L1_levels<-meta$L1_levels; L2_levels<-meta$L2_levels
L3_levels<-meta$L3_levels; FT_levels<-meta$FT_levels
knot1<-meta$bgi_knots[1]; knot2<-meta$bgi_knots[2]
## derived columns (match driver 32d)
pairs[, hg_obs_a := (HT2-HT1)/YEARS]
pairs[, sqrt_years := sqrt(YEARS)]
pairs[, ln_dbh := log(DBH1)]
pairs[, ln_ht := log(pmax(HT1,1.5))]
pairs[, ln_cr_adj := log((CR1+0.2)/1.2)]
pairs[, bal_log := log((BAL_SW1+BAL_HW1)+5)]
if(!"SLOPE"%in%names(pairs)) pairs[,SLOPE:=0]; if(!"ASPECT"%in%names(pairs)) pairs[,ASPECT:=0]
pairs[!is.finite(SLOPE),SLOPE:=0]; pairs[!is.finite(ASPECT),ASPECT:=0]
pairs[, slope_pct := as.numeric(SLOPE)]
pairs[, cos_aspect := cos(as.numeric(ASPECT)*pi/180)]
pairs[, ba_metric := BA1*0.2296]
pairs[, bal_raw := BAL_SW1 + BAL_HW1]
pairs[, rd_additive := sdi_additive1 / SDImax_brms]
filt <- with(pairs,
  is.finite(DBH1) & DBH1>=2.54 & is.finite(HT1) & HT1>1.5 & is.finite(HT2) & HT2>1.5 &
  is.finite(CR1) & CR1>0 & CR1<=1.0 & is.finite(YEARS) & YEARS>=1 & YEARS<=20 &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE!="" & TREESTATUS1==1 & TREESTATUS2==1 &
  is.finite(BAL_SW1) & BAL_SW1>=0 & is.finite(BAL_HW1) & BAL_HW1>=0 &
  is.finite(bgi) & is.finite(BA1) & BA1>=0 & is.finite(rd_additive) & is.finite(bal_raw) &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1>0 & hg_obs_a>0.01 & hg_obs_a<5.0)
d <- pairs[filt]
cat("After filters:",nrow(d),"rows\n")
d <- d[SPCD %in% sp_levels]
d[, sp_idx := match(SPCD,sp_levels)]
d[, L1_idx := match(as.character(EPA_L1_CODE),L1_levels)]
d[, L2_idx := match(as.character(EPA_L2_CODE),L2_levels)]
d[, L3_idx := match(as.character(EPA_L3_CODE),L3_levels)]
d[, FT_idx := match(as.integer(FORTYPCD_cond1),FT_levels)]
d <- d[!is.na(sp_idx)&!is.na(L1_idx)&!is.na(L2_idx)&!is.na(L3_idx)&!is.na(FT_idx)]
cat("After level match:",nrow(d),"rows\n")
gm  <- function(v){ o<-tryCatch(fit$draws(variables=v,format="draws_matrix"),error=function(e)NULL); if(is.null(o))NA_real_ else mean(as.numeric(o)) }
vm  <- function(v,n){ o<-tryCatch(fit$draws(variables=v,format="draws_matrix"),error=function(e)NULL); if(is.null(o))rep(0,n) else colMeans(as.matrix(o)) }
a0<-gm("a0");a1<-gm("a1");a2<-gm("a2");a2q<-gm("a2_quad");a3<-gm("a3");a4<-gm("a4")
a5<-gm("a5");a6<-gm("a6");a7<-gm("a7");a8<-gm("a8");a9a<-gm("a9a");a9b<-gm("a9b");a10<-gm("a10")
a_bard<-gm("a_bard");a_blrd<-gm("a_blrd")
trait_effect<-vm("trait_effect",length(sp_levels)); species_site_slope<-vm("species_site_slope",length(sp_levels))
z_sp<-vm("z_sp",length(sp_levels)); z_L1<-vm("z_L1",length(L1_levels)); z_L2<-vm("z_L2",length(L2_levels))
z_L3<-vm("z_L3",length(L3_levels)); z_FT<-vm("z_FT",length(FT_levels)); z_L1_bgi<-vm("z_L1_bgi",length(L1_levels))
stopifnot(is.finite(a2q), is.finite(a_bard), is.finite(a_blrd))   # confirm v8rd params present
bgi<-d$bgi; bgi_b2<-pmax(bgi-knot1,0); bgi_b3<-pmax(bgi-knot2,0); ln_ht_sq<-d$ln_ht^2
ba_x_rd<-d$ba_metric*d$rd_additive; bal_x_rd<-d$bal_raw*d$rd_additive
b_site<-a4 + z_L1_bgi[d$L1_idx] + species_site_slope[d$sp_idx]
eta <- a0 + trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
  z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
  a1*d$ln_dbh + a2*d$ln_ht + a2q*ln_ht_sq + a3*d$ln_cr_adj +
  b_site*bgi + a9a*bgi_b2 + a9b*bgi_b3 +
  a5*d$bal_log + a6*d$ba_metric + a7*d$slope_pct + a8*d$cos_aspect +
  a10*(bgi*d$bal_log) + a_bard*ba_x_rd + a_blrd*bal_x_rd
obs_raw<-d$hg_obs_a; residual<-log(obs_raw)-eta; weight<-d$sqrt_years
rm<-mean(residual,na.rm=TRUE); rs<-sd(residual,na.rm=TRUE)
cat(sprintf("Residual: n=%d mean=%.4f sd=%.4f p01=%.3f p99=%.3f\n",
    sum(is.finite(residual)),rm,rs,quantile(residual,0.01,na.rm=TRUE),quantile(residual,0.99,na.rm=TRUE)))
## HARD VERIFICATION GATE
if(!is.finite(rm) || abs(rm) > GATE || !is.finite(rs) || rs > 2){
  stop(sprintf("ETA RECONSTRUCTION GATE FAILED: |mean|=%.4f (max %.3f), sd=%.4f. Aborting; downstream modifier fit will NOT run.", abs(rm), GATE, rs))
}
cat(sprintf("GATE PASSED (|mean| %.4f <= %.3f).\n", abs(rm), GATE))
keep<-intersect(c("SPCD","sp_idx","EPA_L1_CODE","L1_idx","YEARS","is_plantation",
  "had_fire_t1","had_insect_t1","had_disease_t1","had_wind_t1","had_harvest_t1",
  "had_cutting_t1","had_site_prep_t1","years_since_dstrb","years_since_trt",
  "dstrb_decay_5yr","dstrb_decay_10yr","dstrb_decay_20yr",
  "trt_decay_5yr","trt_decay_10yr","trt_decay_20yr"),names(d))
out<-d[,..keep]; out[,eta_base:=eta]; out[,obs_raw:=obs_raw]; out[,residual:=residual]; out[,weight:=weight]
for(c in grep("^(is_plantation|had_)",names(out),value=TRUE)){v<-out[[c]];v[is.na(v)]<-0L;out[[c]]<-as.integer(v)}
for(c in grep("_decay_",names(out),value=TRUE)){v<-out[[c]];v[is.na(v)]<-0.0;out[[c]]<-as.numeric(v)}
saveRDS(list(model="hg_unified_v8rd",family="log",fit_path=FIT_FILE,data=out,
  sp_levels=sp_levels,L1_levels=L1_levels,n_rows=nrow(out),resid_sd=rs), OUT_FILE)
cat("\nSaved:",OUT_FILE,"\nsigma_resid:",round(rs,3),"\n")
