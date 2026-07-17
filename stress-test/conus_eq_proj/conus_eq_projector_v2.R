#!/usr/bin/env Rscript
## conus_eq_projector_v2.R -- 100yr (20x5yr) stand projection using fvs-conus fitted equations.
## v2 adds the THREE FIXES on top of the validated v1 engine:
##   FIX 1  IDENTICAL-STAND SEEDING. Seed from the SAME standinit_<VARIANT>.csv +
##          <STATE>_FVS_TREEINIT_PLOT.csv that the engine arms used (keyed STAND_CN),
##          exactly as run_conus_task_wo1.py. Fitted-equation site covariates
##          (bgi, cspi/cspi_v6, EPA_L1/L2/L3, FORTYPCD, SDImax_brms, climate_si) are
##          joined onto each stand by nearest pairs-plot in LAT/LON within state
##          (county-median then state-median fallback). STAND_CN emitted is the
##          engine's standinit STAND_CN so the comparison is on identical stands.
##   FIX 2  INGROWTH/RECRUITMENT. Per-cycle recruits added from the banked, converged
##          per-variant empirical ingrowth lookup (output/comparisons_overstory/
##          intermediate/ingrowth_lookup.rds: med_ann_TPA recruits/ac/yr,
##          med_ann_BA ft2/ac/yr). Recruit DBH derived from BA/TPA so BA & TPA
##          increments are internally consistent; recruits enter the live tree list.
##   FIX 3  SDIMAX SELF-THINNING CAP. Density-dependent mortality multiplier ramps up
##          as relative density (SDI/SDImax) approaches/exceeds 1, so BA plateaus near
##          but not above the physical SDImax ceiling.
## species_mode toggle: --mode=dependent -> conus_b2 (v8 species-aware DG);
##                       --mode=free      -> conus_b1 (speciesfree DG, W*gamma only).
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){ m<-grep(paste0("^--",n,"="),args,value=TRUE); if(!length(m)) return(d); sub(paste0("^--",n,"="),"",m[1]) }
ROOT<-"/users/PUOM0008/crsfaaron/fvs-conus"; SCR<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
OUTD<-ga("outdir",file.path(SCR,"conus_eq_proj")); MODE<-ga("mode","dependent"); VARIANT<-toupper(ga("variant","NE"))
NSTAND<-as.integer(ga("nstands","0")); SEED<-as.integer(ga("seed","7")); NCYC<-as.integer(ga("ncycles","20")); CYCLEN<-as.integer(ga("cyclelen","5"))
STANDINIT_DIR<-ga("standinit_dir",file.path(SCR,"standinit_by_variant"))
TREEINIT_DIR<-ga("treeinit_dir","/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h")
PAIRS<-ga("pairs",file.path(ROOT,"data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"))
TRAITS<-ga("traits",file.path(ROOT,"traits/species_traits.rds"))
IGLOOK<-ga("ingrowth",file.path(ROOT,"output/comparisons_overstory/intermediate/ingrowth_lookup.rds"))
V8DIR<-file.path(ROOT,"output/conus/dg_kue/v8"); V8SUM<-file.path(V8DIR,"dg_kuehne_v8_100k_prod_summary.csv"); V8META<-file.path(V8DIR,"dg_kuehne_v8_100k_prod_meta.rds")
RE_RDS<-file.path(SCR,"smoke_conus_eq/v8_re_means.rds")
B1FIT<-file.path(ROOT,"output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds"); B1META<-file.path(ROOT,"output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_meta.rds")
CC<-file.path(ROOT,"output/conus")
MORT_SUM<-file.path(CC,"mort_logit_simple_cspi_traits1_fixed_summary.csv"); MORT_SP<-file.path(CC,"mort_logit_simple_cspi_traits1_species_intercepts.csv")
HTD_SUM<-file.path(CC,"htdbh_wykoff_lognormal_cspi_traits1_fixed_summary.csv"); HTD_SP<-file.path(CC,"htdbh_wykoff_lognormal_cspi_traits1_species_intercepts.csv")
CR_SUM<-file.path(CC,"cr_recession_cspi_traits1_fixed_summary.csv"); CR_SP<-file.path(CC,"cr_recession_cspi_traits1_species_intercepts.csv")
MORTMODE<-tolower(ga("mort","logit")); stopifnot(MORTMODE %in% c("logit","gompit"))
stopifnot(MODE %in% c("dependent","free")); CONFIG<-if(MODE=="dependent")"conus_b2" else "conus_b1"
if(MORTMODE=="gompit") CONFIG<-paste0(CONFIG,"_gompit")
GOMP_COEF<-ga("gompit_coef","/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv")
CCH_MODULE<-ga("cch_module",file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE)[1])),"cch_module.R"))
if(!file.exists(CCH_MODULE)) CCH_MODULE<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/cch_module.R"
dir.create(OUTD,recursive=TRUE,showWarnings=FALSE)
FIPS<-c("1"="AL","2"="AK","4"="AZ","5"="AR","6"="CA","8"="CO","9"="CT","10"="DE","12"="FL","13"="GA","16"="ID","17"="IL","18"="IN","19"="IA","20"="KS","21"="KY","22"="LA","23"="ME","24"="MD","25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT","31"="NE","32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND","39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD","47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV","55"="WI","56"="WY")
cat("==== conus_eq_projector_v2 (3 fixes) ====\n"); cat(sprintf("  variant=%s mode=%s (CONFIG=%s) nstands=%s %dx%dyr\n",VARIANT,MODE,CONFIG,if(NSTAND>0)NSTAND else "ALL",NCYC,CYCLEN))
cat("  WIRED: DG(b2/b1), mortality(logit), ht-dbh(Wykoff), CR-recession, dynamic competition.\n")
cat("  FIXES: (1) identical-stand seeding standinit+treeinit, (2) ingrowth lookup, (3) SDImax self-thinning cap.\n\n")
for(f in c(PAIRS,TRAITS,V8SUM,V8META,RE_RDS,MORT_SUM,MORT_SP,HTD_SUM,HTD_SP,CR_SUM,CR_SP,IGLOOK)) if(!file.exists(f)) stop("MISSING: ",f)
if(MODE=="free") for(f in c(B1FIT,B1META)) if(!file.exists(f)) stop("MISSING b1: ",f)
if(MORTMODE=="gompit"){
  if(!file.exists(GOMP_COEF)) stop("MISSING gompit coefficients: ",GOMP_COEF)
  if(!file.exists(CCH_MODULE)) stop("MISSING cch_module.R: ",CCH_MODULE)
  source(CCH_MODULE)
  gcoef<-fread(GOMP_COEF)
  gomp_b0<-setNames(gcoef$b0,as.character(gcoef$SPCD)); gomp_b1<-setNames(gcoef$b1,as.character(gcoef$SPCD))
  gomp_b2<-setNames(gcoef$b2,as.character(gcoef$SPCD)); gomp_b3<-setNames(gcoef$b3,as.character(gcoef$SPCD))
  gomp_b4<-setNames(gcoef$b4,as.character(gcoef$SPCD))
  ## genus/hardwood-softwood fallback for species without a fit: use the
  ## conifer-median or hardwood-median parameter vector (Greg fits 133 species).
  sw_med<-function(col,sw){ v<-gcoef[[col]][ (gcoef$SPCD<300)==sw ]; median(v,na.rm=TRUE) }
  GFB<-list(sw=c(b0=sw_med("b0",TRUE),b1=sw_med("b1",TRUE),b2=sw_med("b2",TRUE),b3=sw_med("b3",TRUE),b4=sw_med("b4",TRUE)),
            hw=c(b0=sw_med("b0",FALSE),b1=sw_med("b1",FALSE),b2=sw_med("b2",FALSE),b3=sw_med("b3",FALSE),b4=sw_med("b4",FALSE)))
  cat(sprintf("  GOMPIT mortality: %d species banked; fallback SW b0=%.2f HW b0=%.2f; CCH via ORGANON port (affine A=%.3f B=%.4f)\n",
      nrow(gcoef),GFB$sw["b0"],GFB$hw["b0"],CCH_A,CCH_B))
  ## per-tree gompit survival over a T-year cycle: P_surv=exp(-exp(eta)*T)
  gomp_surv<-function(SPCD,cr,cch,Tyr){
    k<-as.character(SPCD); n<-length(SPCD)
    b0<-gomp_b0[k]; b1<-gomp_b1[k]; b2<-gomp_b2[k]; b3<-gomp_b3[k]; b4<-gomp_b4[k]
    miss<-is.na(b0); if(any(miss)){ sw<-SPCD[miss]<300
      pick<-function(p) ifelse(sw,GFB$sw[p],GFB$hw[p])
      b0[miss]<-pick("b0"); b1[miss]<-pick("b1"); b2[miss]<-pick("b2"); b3[miss]<-pick("b3"); b4[miss]<-pick("b4") }
    cr<-pmin(pmax(cr,1e-4),1); cch<-pmax(cch,0)
    cr_term<-(cr+0.01)^b2; cch_term<-ifelse(cch>0,cch^b4,0)
    eta<-b0+b1*cr_term+b3*cch_term; eta<-pmin(pmax(eta,-30),30)
    H<-exp(eta); surv<-exp(-H*Tyr); surv[!is.finite(surv)]<-1
    pmin(pmax(surv,0),1) }
}

SI_FILE<-file.path(STANDINIT_DIR,paste0("standinit_",VARIANT,".csv")); if(!file.exists(SI_FILE)) stop("MISSING standinit: ",SI_FILE)

read_fx<-function(p){s<-fread(p); setNames(s$mean,s$variable)}
fx8<-read_fx(V8SUM); b8<-function(n) as.numeric(fx8[[n]]); RE<-readRDS(RE_RDS); m8<-readRDS(V8META)
sp_levels<-m8$sp_levels; L1_lev<-m8$L1_levels; L2_lev<-m8$L2_levels; L3_lev<-m8$L3_levels; FT_lev<-m8$FT_levels; tcols<-m8$trait_cols
k1<-m8$bgi_knots[1]; k2<-m8$bgi_knots[2]; sigma8<-b8("sigma")
fxm<-read_fx(MORT_SUM); mc<-function(n) as.numeric(fxm[[n]]); mort_sp<-fread(MORT_SP); mort_zsp<-setNames(mort_sp$mean,as.character(mort_sp$SPCD))
fxh<-read_fx(HTD_SUM); hc<-function(n) as.numeric(fxh[[n]]); htd_sp_dt<-fread(HTD_SP); htd_zsp<-setNames(htd_sp_dt$mean,as.character(htd_sp_dt$SPCD))
fxc<-read_fx(CR_SUM); rc<-function(n) as.numeric(fxc[[n]]); cr_sp_dt<-fread(CR_SP); cr_zsp<-setNames(cr_sp_dt$mean,as.character(cr_sp_dt$SPCD))
traits<-as.data.table(readRDS(TRAITS))
traits_sub<-traits[match(sp_levels,SPCD),c("SPCD",tcols),with=FALSE]; W<-as.matrix(traits_sub[,tcols,with=FALSE])
sw_by_sp<-traits_sub$softwood; sw_by_sp[is.na(sw_by_sp)]<-0
for(j in seq_len(ncol(W))){na<-is.na(W[,j]); if(any(na)) W[na,j]<-median(W[!na,j],na.rm=TRUE); W[,j]<-(W[,j]-mean(W[,j]))/sd(W[,j])}
hc_gamma<-{gn<-grep("^gamma\\[",names(fxh),value=TRUE); if(length(gn)) as.numeric(fxh[gn]) else rep(0,ncol(W))}
if(length(hc_gamma)!=ncol(W)) hc_gamma<-rep(0,ncol(W))

## ---- FIX 2: per-variant ingrowth rates (banked empirical lookup) ----
ig_lookup<-readRDS(IGLOOK); ig<-ig_lookup[[VARIANT]]; if(is.null(ig)) ig<-ig_lookup[["OVERALL"]]
IG_TPA<-if(!is.null(ig)) as.numeric(ig$med_ann_TPA) else 0   # recruits/ac/yr
IG_BA <-if(!is.null(ig)) as.numeric(ig$med_ann_BA)  else 0   # ft2/ac/yr
cat(sprintf("FIX2 ingrowth (%s): med_ann_TPA=%.3f rec/ac/yr  med_ann_BA=%.3f ft2/ac/yr\n",VARIANT,IG_TPA,IG_BA))
## recruit mean DBH (in) from BA/TPA identity: ba_tree=BA/TPA(ft2); dbh=sqrt(ba/(pi/4)/(1/144))
ig_dbh_in<-if(IG_TPA>1e-6 && IG_BA>0){ ba_tree<-IG_BA/IG_TPA; sqrt(ba_tree/(pi/4)*144) } else 1.5
ig_dbh_in<-min(max(ig_dbh_in,1.0),3.0)   # recruits are small trees (1-3 in)
cat(sprintf("       recruit DBH=%.2f in (small-tree entry); recruits added each cycle.\n",ig_dbh_in))

b1_pack<-NULL
load_b1<-function(){
  B1SUM<-"/users/PUOM0008/crsfaaron/fvs-conus/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_summary.csv"
  B1RE<-file.path(OUTD,"b1_re_means.rds")
  m<-readRDS(B1META); pm<-m$prep_meta
  s<-fread(B1SUM); gv<-setNames(s$mean,s$variable)
  if(file.exists(B1RE)){ RErb<-readRDS(B1RE)
  } else { cat("  b1_re_means.rds not found; loading 4.7GB fit one-time...\n")
    fit<-readRDS(B1FIT); pull<-function(v) as.numeric(fit$summary(v,"mean")$mean)
    RErb<-list(z_L1=pull("z_L1"),z_L2=pull("z_L2"),z_L3=pull("z_L3")); rm(fit); gc(verbose=FALSE)
    saveRDS(RErb,B1RE) }
  list(pm=pm,tcols=m$trait_cols,g=function(n) as.numeric(gv[[n]]),
    z_L1=RErb$z_L1,z_L2=RErb$z_L2,z_L3=RErb$z_L3,
    gamma=as.numeric(gv[grep("^gamma\\[",names(gv))]))}
if(MODE=="free") b1_pack<-load_b1()
sw_mean<-mean(sw_by_sp); sw_tree_map<-setNames(traits$softwood,as.character(traits$SPCD))
## modal SPCD for recruits = most common softwood-balanced species in NE? use stand's own dominant later.

## ============================ FIX 1: SEEDING ============================
cat("STEP 1: load standinit + treeinit (engine identical stands) + join covariates ...\n")
si<-fread(SI_FILE,colClasses=list(character=c("STAND_CN","STAND_ID")))
si[,STAND_CN:=sub("\\..*$","",STAND_CN)]
si[,INV_YEAR:=suppressWarnings(as.integer(INV_YEAR))]; si[is.na(INV_YEAR),INV_YEAR:=2010L]
si[,STATE:=suppressWarnings(as.integer(STATE))]
si<-si[is.finite(STATE)&!is.na(STAND_CN)&STAND_CN!=""]
## build covariate plot table from pairs (one row per FIA plot, any variant -> site covars)
d<-as.data.table(readRDS(PAIRS))
pl<-d[,.(bgi=first(bgi),cspi=first(cspi),cspi_v6=first(cspi_v6),SDImax_brms=first(SDImax_brms),
  climate_si=first(climate_si),
  EPA_L1_CODE=as.character(first(EPA_L1_CODE)),EPA_L2_CODE=as.character(first(EPA_L2_CODE)),
  EPA_L3_CODE=as.character(first(EPA_L3_CODE)),FORTYPCD_cond1=as.character(first(FORTYPCD_cond1)),
  LAT=first(LAT),LON=first(LON)),by=.(STATECD,COUNTYCD,PLOT)]
modef<-function(x){x<-as.character(x);x<-x[!is.na(x)&x!=""]; if(!length(x))return(NA_character_); names(sort(table(x),decreasing=TRUE))[1]}
## county-median + state-median fallbacks
cty<-pl[,.(bgi=median(bgi,na.rm=TRUE),cspi=median(cspi,na.rm=TRUE),cspi_v6=median(cspi_v6,na.rm=TRUE),
  SDImax_brms=median(SDImax_brms,na.rm=TRUE),climate_si=median(climate_si,na.rm=TRUE),
  EPA_L1_CODE=modef(EPA_L1_CODE),EPA_L2_CODE=modef(EPA_L2_CODE),EPA_L3_CODE=modef(EPA_L3_CODE),
  FORTYPCD_cond1=modef(FORTYPCD_cond1)),by=.(STATECD,COUNTYCD)]
sta<-pl[,.(bgi=median(bgi,na.rm=TRUE),cspi=median(cspi,na.rm=TRUE),cspi_v6=median(cspi_v6,na.rm=TRUE),
  SDImax_brms=median(SDImax_brms,na.rm=TRUE),climate_si=median(climate_si,na.rm=TRUE),
  EPA_L1_CODE=modef(EPA_L1_CODE),EPA_L2_CODE=modef(EPA_L2_CODE),EPA_L3_CODE=modef(EPA_L3_CODE),
  FORTYPCD_cond1=modef(FORTYPCD_cond1)),by=.(STATECD)]
si[,sid_county:=suppressWarnings(as.integer(substr(STAND_ID,5,7)))]
numcov<-c("bgi","cspi","cspi_v6","SDImax_brms","climate_si")
chrcov<-c("EPA_L1_CODE","EPA_L2_CODE","EPA_L3_CODE","FORTYPCD_cond1")
covnm<-c(numcov,chrcov)
## nearest-plot fill in LAT/LON within state (best per-stand site values & EPA codes)
## initialise with correct column types so set() below does not coerce
for(cc in numcov) si[,(cc):=NA_real_]
for(cc in chrcov) si[,(cc):=NA_character_]
plf<-pl[is.finite(LAT)&is.finite(LON)]
for(st in unique(si$STATE)){
  idx<-which(si$STATE==st & is.finite(si$LATITUDE) & is.finite(si$LONGITUDE)); if(!length(idx)) next
  cand<-plf[STATECD==st]; if(!nrow(cand)){
    cand<-plf[abs(LAT-mean(si$LATITUDE[idx]))<3 & abs(LON-mean(si$LONGITUDE[idx]))<5] }
  if(!nrow(cand)) next
  la<-si$LATITUDE[idx]; lo<-si$LONGITUDE[idx]
  nn<-sapply(seq_along(idx),function(i){dx<-cand$LAT-la[i];dy<-cand$LON-lo[i]; which.min(dx*dx+dy*dy)})
  for(cc in covnm) set(si,i=idx,j=cc,value=cand[[cc]][nn])
}
## county-median fallback for any still-NA
si[,kc:=paste(STATE,sid_county,sep="_")]; cty[,kc:=paste(STATECD,COUNTYCD,sep="_")]
mc_<-cty[match(si$kc,cty$kc)]
for(cc in covnm){ nas<-is.na(si[[cc]]); if(any(nas)) set(si,i=which(nas),j=cc,value=mc_[[cc]][which(nas)]) }
## state-median fallback
ms_<-sta[match(si$STATE,sta$STATECD)]
for(cc in covnm){ nas<-is.na(si[[cc]]); if(any(nas)) set(si,i=which(nas),j=cc,value=ms_[[cc]][which(nas)]) }
cat(sprintf("  standinit stands: %d  covariate coverage bgi=%.1f%% EPA_L1=%.1f%% SDImax=%.1f%%\n",
  nrow(si),100*mean(is.finite(si$bgi)),100*mean(!is.na(si$EPA_L1_CODE)),100*mean(is.finite(si$SDImax_brms))))

## stand subset selection (idempotent: deterministic by NSTAND>0)
si<-si[is.finite(bgi)&!is.na(EPA_L1_CODE)&is.finite(SDImax_brms)]
if(NSTAND>0 && NSTAND<nrow(si)){ set.seed(SEED); si<-si[sort(sample(.N,NSTAND))] }
cat("  stands to project:",nrow(si),"\n")

## treeinit cache per state
tcache<-new.env()
get_state_trees<-function(stcode){
  st<-FIPS[as.character(stcode)]; if(is.na(st)) return(NULL)
  if(!is.null(tcache[[st]])) return(tcache[[st]])
  tf<-file.path(TREEINIT_DIR,paste0(st,"_FVS_TREEINIT_PLOT.csv")); if(!file.exists(tf)){tcache[[st]]<-list(NULL);return(NULL)}
  tt<-fread(tf,colClasses=list(character=c("STAND_CN")),select=c("STAND_CN","TREE_COUNT","SPECIES","DIAMETER","HT","CRRATIO"))
  tt[,STAND_CN:=sub("\\..*$","",STAND_CN)]
  spl<-split(tt,by="STAND_CN"); tcache[[st]]<-spl; spl }

## ============================ EQUATION KERNELS ============================
eta_dg_b2<-function(tl){ln_dbh<-log(tl$dbh_in); ln_cr<-log((tl$CR+0.2)/1.2); BALsw_m<-tl$BAL_SW*0.2296; BALhw_m<-tl$BAL_HW*0.2296; ln_balsw<-log(BALsw_m+0.01)
  cb1<-tl$bgi; cb2<-pmax(tl$bgi-k1,0); cb3<-pmax(tl$bgi-k2,0); spi<-tl$sp_idx; L1i<-tl$L1_idx; L2i<-tl$L2_idx; L3i<-tl$L3_idx; FTi<-tl$FT_idx; sw_c<-tl$softwood-tl$sw_mean
  b_site<-b8("b6")+RE$z_L1_bgi[L1i]+RE$species_site_slope[spi]
  eta<-b8("b0")+RE$trait_effect[spi]+RE$z_sp[spi]+RE$z_L1[L1i]+RE$z_L2[L2i]+RE$z_L3[L3i]+RE$z_FT[FTi]+b8("b1")*ln_dbh+b8("b2")*tl$dbh_in+b8("b3")*ln_cr+b8("b4")*ln_balsw+b8("b5")*BALhw_m+b_site*cb1+b8("b9a")*cb2+b8("b9b")*cb3+b8("b7")*(tl$BA*0.2296*tl$rd_add)+b8("b8")*(BALsw_m*tl$rd_add)+b8("b11")*tl$sdi_cx+b8("b12")*(tl$bgi*tl$rd_add)+b8("b13")*(tl$bgi*ln_dbh)+b8("b14")*(tl$bgi*sw_c)+b8("b15")*(tl$bgi*ln_cr)
  pmin(pmax(eta,-30),20)}
eta_dg_b1<-function(tl,P){g<-P$g; pm<-P$pm; ln_dbh<-log(tl$dbh_in); ln_cr<-log((tl$CR+0.2)/1.2); ln_balsw<-log(tl$BAL_SW*0.2296+0.01); ln_csi<-log(tl$cspi+pm$cspi_shift)
  ba_x_rd<-(tl$BA*0.2296)*tl$rd_add; bal_x_rd<-(tl$BAL*0.2296)*tl$rd_add; trait_effect<-as.numeric(tl$Wrow %*% P$gamma)
  eta<-g("b0")+trait_effect+P$z_L1[tl$b1_L1]+P$z_L2[tl$b1_L2]+P$z_L3[tl$b1_L3]+g("b1")*ln_dbh+g("b2")*tl$dbh_in+g("b3")*ln_cr+g("b4")*ln_balsw+g("b5")*(tl$BAL_HW*0.2296)+g("b6")*ln_csi+g("b7")*ba_x_rd+g("b8")*bal_x_rd
  pmin(pmax(eta,-30),20)}
eta_mort<-function(tl){dbh<-tl$dbh_in; dbh2<-dbh^2; bal_over_ba<-ifelse(tl$BA>0,tl$BAL/tl$BA,0); sqrt_ba_rd<-sqrt(pmax((tl$BA*0.2296)*tl$rd_add,0)); ln_csi<-log(tl$cspi+tl$mort_cspi_shift); zsp<-tl$mort_zsp
  mc("m0")+zsp+mc("m1")*dbh+mc("m2")*dbh2+mc("m3")*bal_over_ba+mc("m4")*tl$CR+mc("m5")*sqrt_ba_rd+mc("m6")*ln_csi}
ht_from_dbh<-function(tl){BAm<-tl$BA*0.2296; BALm<-tl$BAL*0.2296; bal<-BALm; sqrt_ba<-sqrt(pmax(BAm,0)); ln_csi<-log(tl$cspi+tl$htd_cspi_shift); ba_x_rd<-BAm*tl$rd_add; bal_x_rd<-BALm*tl$rd_add
  te<-as.numeric(tl$Wrow %*% hc_gamma); te[is.na(te)]<-0; zsp<-tl$htd_zsp
  eta<-hc("b0")+te+zsp+hc("a_bal")*bal+hc("a_ba")*sqrt_ba+hc("a_cspi")*ln_csi+hc("a_bard")*ba_x_rd+hc("a_blrd")*bal_x_rd+hc("b1")/(tl$dbh_in+1.0)
  1.37+exp(eta+0.5*hc("sigma")^2)}
cr_update<-function(tl,HT2){cr<-tl$CR; ln_cr<-log(pmax(cr,1e-4)); sqrt_ba<-sqrt(pmax(tl$BA*0.2296,0)); ln_bal_ba<-log(tl$BAL/pmax(tl$BA,1e-6)+1); cr_over_rd<-cr/pmax(tl$rd_add,1e-4); ln_csi<-log(tl$cspi+tl$cr_cspi_shift); zsp<-tl$cr_zsp
  eta<-rc("r0")+zsp+rc("r1")*ln_cr+rc("r2")*cr+rc("r3")*sqrt_ba+rc("r4")*ln_bal_ba+rc("r5")*cr_over_rd+rc("r6")*ln_csi
  r<-1/(1+exp(eta)); HCB1<-(1-cr)*tl$HT; maxrec<-pmax(HT2-HCB1,0); HCB2<-HCB1+r*maxrec; CR2<-1-HCB2/pmax(HT2,1e-3); pmin(pmax(CR2,0.01),0.95)}
mort_cspi_shift<-if(!is.null(b1_pack)) b1_pack$pm$cspi_shift else 1.0; htd_cspi_shift<-mort_cspi_shift; cr_cspi_shift<-mort_cspi_shift

## build a tree list from a treeinit data.table for one stand + its joined covariates
mk_tl<-function(trows,cov){
  d2<-trows[is.finite(DIAMETER)&DIAMETER>=1.0]
  if(!nrow(d2)) return(NULL)
  SPCD<-as.integer(d2$SPECIES); dbh_in<-as.numeric(d2$DIAMETER)
  cr<-as.numeric(d2$CRRATIO); cr[!is.finite(cr)|cr<=0]<-NA   # CRRATIO is 0-100 (percent) in FVS treeinit
  cr<-cr/100; cr[is.na(cr)]<-0.5; cr<-pmin(pmax(cr,0.05),0.95)
  ht<-as.numeric(d2$HT)*0.3048; ht[!is.finite(ht)|ht<=1.37]<-NA   # ft -> m
  TPA<-as.numeric(d2$TREE_COUNT); TPA[!is.finite(TPA)|TPA<=0]<-1.0
  tl<-list(SPCD=SPCD,sp_idx=match(SPCD,sp_levels),
    L1_idx=match(as.character(cov$EPA_L1_CODE),L1_lev),L2_idx=match(as.character(cov$EPA_L2_CODE),L2_lev),
    L3_idx=match(as.character(cov$EPA_L3_CODE),L3_lev),FT_idx=match(as.integer(cov$FORTYPCD_cond1),FT_lev),
    dbh_in=dbh_in,CR=cr,HT=ht,TPA=TPA,
    softwood={v<-sw_tree_map[as.character(SPCD)];v[is.na(v)]<-0;as.integer(v)},
    bgi=cov$bgi,cspi=cov$cspi,SDImax=cov$SDImax_brms,sw_mean=sw_mean,
    mort_zsp={v<-mort_zsp[as.character(SPCD)];v[is.na(v)]<-0;v},
    htd_zsp={v<-htd_zsp[as.character(SPCD)];v[is.na(v)]<-0;v},
    cr_zsp={v<-cr_zsp[as.character(SPCD)];v[is.na(v)]<-0;v},
    mort_cspi_shift=mort_cspi_shift,htd_cspi_shift=htd_cspi_shift,cr_cspi_shift=cr_cspi_shift)
  ## L1/L2/L3/FT are scalars from cov; expand later via index. EPA idx fallback to 1 if unmatched.
  for(nm in c("L1_idx","L2_idx","L3_idx","FT_idx")) if(is.na(tl[[nm]])) tl[[nm]]<-1L
  tl$Wrow<-W[tl$sp_idx,,drop=FALSE]
  if(MODE=="free"){tl$b1_L1<-match(as.character(cov$EPA_L1_CODE),b1_pack$pm$L1);tl$b1_L2<-match(as.character(cov$EPA_L2_CODE),b1_pack$pm$L2);tl$b1_L3<-match(as.character(cov$EPA_L3_CODE),b1_pack$pm$L3)
    tl$b1_L1<-if(is.na(tl$b1_L1))1L else tl$b1_L1; tl$b1_L2<-if(is.na(tl$b1_L2))1L else tl$b1_L2; tl$b1_L3<-if(is.na(tl$b1_L3))1L else tl$b1_L3}
  ## KEEP all trees (engine projects all). Trees whose SPCD is not in the DG sp_levels
  ## (no growth params) are remapped to the stand's modal in-set species for the GROWTH
  ## kernels (sp_idx/Wrow/z-effects/softwood), but their ORIGINAL SPCD is preserved for
  ## NSBE biomass (NSBE covers 465 species, a superset). This avoids the year-0 biomass
  ## undershoot from dropping ~4% of trees.
  unk<-is.na(tl$sp_idx)
  if(any(unk)){
    inset<-tl$sp_idx[!unk]
    fb_idx<-if(length(inset)) as.integer(names(sort(table(inset),decreasing=TRUE))[1]) else match(sp_levels[1],sp_levels)
    fb_spcd<-sp_levels[fb_idx]
    tl$sp_idx[unk]<-fb_idx
    tl$Wrow[unk,]<-matrix(W[fb_idx,],nrow=sum(unk),ncol=ncol(W),byrow=TRUE)
    tl$softwood[unk]<-{v<-sw_tree_map[as.character(fb_spcd)];if(is.na(v))0L else as.integer(v)}
    tl$mort_zsp[unk]<-{v<-mort_zsp[as.character(fb_spcd)];if(is.na(v))0 else v}
    tl$htd_zsp[unk]<-{v<-htd_zsp[as.character(fb_spcd)];if(is.na(v))0 else v}
    tl$cr_zsp[unk]<-{v<-cr_zsp[as.character(fb_spcd)];if(is.na(v))0 else v}
    ## SPCD (original) kept as-is for biomass
  }
  ## expand scalar indices to vector length for kernel use
  n<-length(tl$dbh_in); if(!n) return(NULL)
  for(nm in c("L1_idx","L2_idx","L3_idx","FT_idx","bgi","cspi","SDImax","sw_mean","mort_cspi_shift","htd_cspi_shift","cr_cspi_shift"))
    tl[[nm]]<-rep(tl[[nm]],length.out=1)  # keep scalar; kernels broadcast
  if(MODE=="free") for(nm in c("b1_L1","b1_L2","b1_L3")) tl[[nm]]<-rep(tl[[nm]],length.out=1)
  tl}

## recompute competition; SDImax & EPA indices held in tl (scalars broadcast in kernels)
recompute_comp<-function(tl){DBH<-tl$dbh_in; TPA<-tl$TPA; n<-length(DBH); ba_tree<-pi/4*DBH^2/144; BA_ac<-ba_tree*TPA; BA_tot<-sum(BA_ac); ord<-order(DBH,decreasing=TRUE)
  BAL<-numeric(n); BAL_SW<-numeric(n); BAL_HW<-numeric(n); cum<-0; cum_sw<-0; cum_hw<-0
  for(i in ord){BAL[i]<-cum; BAL_SW[i]<-cum_sw; BAL_HW[i]<-cum_hw; cum<-cum+BA_ac[i]; if(tl$softwood[i]==1) cum_sw<-cum_sw+BA_ac[i] else cum_hw<-cum_hw+BA_ac[i]}
  SDI<-sum(TPA*(DBH/10)^1.605); rd_add<-SDI/pmax(tl$SDImax,1); sdi_cx<-SDI/pmax(SDI,1.0)
  tl$BA<-BA_tot; tl$BAL<-BAL; tl$BAL_SW<-BAL_SW; tl$BAL_HW<-BAL_HW; tl$SDI<-SDI; tl$rd_add<-rd_add; tl$sdi_cx<-sdi_cx; tl}
stand_metrics<-function(tl){DBH<-tl$dbh_in; TPA<-tl$TPA; ok<-is.finite(DBH)&is.finite(TPA)&TPA>0; DBH<-DBH[ok]; TPA<-TPA[ok]; sT<-sum(TPA); list(BA=if(length(DBH)) sum(pi/4*DBH^2*TPA)/144 else NA_real_, QMD=if(is.finite(sT)&&sT>1e-6) sqrt(sum(DBH^2*TPA)/sT) else NA_real_, TPH=if(is.finite(sT)) sT*2.4710538 else NA_real_)}

## FIX 2 helper: add recruits to a tree list (recruits inherit stand dominant species / softwood)
add_recruits<-function(tl,n_add_tpa){
  if(n_add_tpa<=1e-6) return(tl)
  ## recruit species = stand dominant by BA (use existing species so DG/mort params apply)
  if(length(tl$dbh_in)){ ba_ac<-pi/4*tl$dbh_in^2/144*tl$TPA; dom<-tl$SPCD[which.max(ba_ac)] } else dom<-sp_levels[1]
  spi<-match(dom,sp_levels); if(is.na(spi)) return(tl)
  tl$SPCD<-c(tl$SPCD,dom); tl$sp_idx<-c(tl$sp_idx,spi)
  tl$dbh_in<-c(tl$dbh_in,ig_dbh_in); tl$CR<-c(tl$CR,0.6); tl$HT<-c(tl$HT,NA_real_)
  tl$TPA<-c(tl$TPA,n_add_tpa)
  swv<-sw_tree_map[as.character(dom)]; swv<-if(is.na(swv))0L else as.integer(swv)
  tl$softwood<-c(tl$softwood,swv)
  tl$mort_zsp<-c(tl$mort_zsp,{v<-mort_zsp[as.character(dom)];if(is.na(v))0 else v})
  tl$htd_zsp<-c(tl$htd_zsp,{v<-htd_zsp[as.character(dom)];if(is.na(v))0 else v})
  tl$cr_zsp<-c(tl$cr_zsp,{v<-cr_zsp[as.character(dom)];if(is.na(v))0 else v})
  if(!is.null(tl$ht_ratio)) tl$ht_ratio<-c(tl$ht_ratio,1)
  tl$Wrow<-rbind(tl$Wrow,W[spi,,drop=FALSE])
  tl}

cat("STEP 2: projecting",nrow(si),"stands x",NCYC,"cycles ...\n")
rows<-vector("list",nrow(si)); treerows<-vector("list",nrow(si)); ri<-0; ti<-0; nproj<-0
state_groups<-split(seq_len(nrow(si)),si$STATE)
for(stcode in names(state_groups)){
  spl<-get_state_trees(as.integer(stcode)); if(is.null(spl)) next
  if(is.list(spl)&&length(spl)==1&&is.null(spl[[1]])) next
  for(ix in state_groups[[stcode]]){
    cn<-si$STAND_CN[ix]; trows<-spl[[cn]]; if(is.null(trows)||!nrow(trows)) next
    cov<-si[ix]; st<-FIPS[as.character(cov$STATE)]; inv_year<-cov$INV_YEAR
    tl<-mk_tl(trows,cov); if(is.null(tl)||length(tl$dbh_in)<1) next
    nproj<-nproj+1; tl<-recompute_comp(tl)
    ## seed HT: keep MEASURED treeinit HT where present (engine reports year-0 AGB on
    ## measured HT); predict via ht-dbh for trees missing HT. Track a per-tree
    ## calibration ratio (measured/predicted) so projected-cycle HT stays anchored to
    ## the measured tree instead of jumping to the population ht-dbh curve.
    hp0<-ht_from_dbh(tl); miss<-!is.finite(tl$HT); if(any(miss)) tl$HT[miss]<-hp0[miss]
    tl$ht_ratio<-pmin(pmax(tl$HT/pmax(hp0,1e-3),0.33),3.0); tl$ht_ratio[miss]<-1
    emit<-function(cy){sm<-stand_metrics(tl); py<-cy*CYCLEN; ri<<-ri+1
      cchm<-if(MORTMODE=="gompit" && !is.null(tl$cch_last) && length(tl$cch_last)==length(tl$TPA) && sum(tl$TPA)>0) sum(tl$cch_last*tl$TPA)/sum(tl$TPA) else NA_real_
      rows[[ri]]<<-data.table(STAND_CN=cn,STATE=st,YEAR=inv_year+py,PROJ_YEAR=py,VARIANT=VARIANT,CONFIG=CONFIG,AGB_TONS_AC=NA_real_,BA_FT2AC=sm$BA,QMD_IN=sm$QMD,TPH=sm$TPH,CCH_MEAN=cchm)
      ti<<-ti+1; treerows[[ti]]<<-data.table(STAND_CN=cn,CONFIG=CONFIG,PROJ_YEAR=py,SPCD=tl$SPCD,DBH_IN=tl$dbh_in,HT_M=tl$HT,TPA=tl$TPA)}
    emit(0)
    for(cy in 1:NCYC){
      eta<-if(MODE=="dependent") eta_dg_b2(tl) else eta_dg_b1(tl,b1_pack)
      dg_a<-exp(eta+sigma8^2/2); dg_a[!is.finite(dg_a)]<-0; dg_a<-pmin(dg_a,2.0)
      tl$dbh_in<-pmin(pmax(tl$dbh_in+dg_a*CYCLEN,0.1),200)
      HT2<-ht_from_dbh(tl)*tl$ht_ratio; tl$CR<-cr_update(tl,HT2); tl$HT<-HT2
      ## ---- MORTALITY: logit (default) OR Greg gompit (--mort=gompit) ----
      rd<-tl$rd_add
      st_mult<-pmin(1+pmax(rd-0.55,0)/0.45*2.0, 8)   # FIX 3 SDImax ramp (shared by both arms)
      if(MORTMODE=="gompit"){
        ## CCH per tree this cycle via the ORGANON crown-closure port (HT m->ft),
        ## crown group = refined genus crosswalk; cch on the gompit (CCH1) scale.
        gg<-grp_organon(tl$SPCD)
        cch<-stand_cch(tl$dbh_in, tl$HT/0.3048, tl$CR, tl$TPA, gg)
        tl$cch_last<-cch
        surv_cyc<-gomp_surv(tl$SPCD, tl$CR, cch, CYCLEN)        # cycle-length survival
        ## fold the shared SDImax self-thinning ramp into the hazard for parity
        H_cyc<- -log(pmin(pmax(surv_cyc,1e-9),1))                # cycle hazard*T
        surv<-exp(-H_cyc*st_mult); surv[!is.finite(surv)]<-1
      } else {
        p_die_a<-1/(1+exp(-eta_mort(tl))); p_die_a[!is.finite(p_die_a)]<-0
        p_die_a<-pmin(pmax(p_die_a*st_mult,0),0.999)
        surv<-(1-p_die_a)^CYCLEN
      }
      tl$TPA<-tl$TPA*surv; tl$TPA[!is.finite(tl$TPA)]<-0
      ## drop dead-empty trees (TPA collapsed)
      keep<-tl$TPA>1e-4
      if(any(!keep)){ vn<-length(tl$dbh_in); for(nm in names(tl)) if(length(tl[[nm]])==vn) tl[[nm]]<-tl[[nm]][keep]; tl$Wrow<-tl$Wrow[keep,,drop=FALSE] }
      ## FIX 2: add ingrowth recruits for this cycle (per-acre rate x cycle length)
      tl<-add_recruits(tl, IG_TPA*CYCLEN)
      tl$HT[!is.finite(tl$HT)]<-ht_from_dbh(tl)[!is.finite(tl$HT)]
      tl<-recompute_comp(tl); emit(cy)}
  }
}
out<-rbindlist(rows[seq_len(ri)]); tre<-rbindlist(treerows[seq_len(ti)]); otag<-sprintf("conus_eq_%s_%s",tolower(VARIANT),CONFIG)
fwrite(out,file.path(OUTD,paste0(otag,"_metrics.csv"))); fwrite(tre,file.path(OUTD,paste0(otag,"_treelists.csv")))
cat("\n  stands projected:",nproj,"\n")
cat("Wrote:",file.path(OUTD,paste0(otag,"_metrics.csv")),"(",nrow(out),"rows )\n"); cat("Wrote:",file.path(OUTD,paste0(otag,"_treelists.csv")),"(",nrow(tre),"rows )\n")
cat("\n=== year 0/50/100 stand-mean metrics (",CONFIG,") ===\n")
print(out[PROJ_YEAR %in% c(0,50,100),.(BA=mean(BA_FT2AC,na.rm=TRUE),QMD=mean(QMD_IN,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE),n=.N),by=PROJ_YEAR])
cat("DONE.\n")
