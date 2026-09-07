#!/usr/bin/env Rscript
## conus_eq_projector_v2.R -- 100yr (20x5yr) stand projection using fvs-conus fitted equations.
## ---- v2run PATCH (2026-07-21, height-emission + v2 size-aware mortality) ----
##   (A) TOP HEIGHT EMISSION. stand_metrics() now returns HT_DOM_M (TPA-weighted mean
##       height of the largest-DBH trees up to 100 TPA ac-1, Assmann dominant height)
##       and HT_MEAN_M; the per-step metrics table carries HT_DOM_M, HT_MEAN_M and an
##       explicit TOPHT alias (=HT_DOM_M). Makes the Eichhorn relation testable per arm.
##   (B) v2 SIZE-AWARE GOMPIT MORTALITY. Default --gompit_coef now points at
##       conus_mort/experimental/greg_mortality_coefficients_sizeaware_v2.csv and gomp_surv()
##       is rewired to the exact deployed GOMPSURV inverted-survival convention those
##       coefficients were fit to (fit_sizeaware_v2.R): DBH in INCHES,
##       eta=b0+b1*(cr+0.01)^b2+b3*cch^b4+b5*log(dbh); pa=1-exp(-exp(eta)); surv=pa^T.
##       (The old exp(-exp(eta)*T) hazard math + b5-drop would misuse these coeffs.)
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
GOMP_COEF<-ga("gompit_coef","/fs/scratch/PUOM0008/crsfaaron/conus_mort/experimental/greg_mortality_coefficients_sizeaware_v2.csv")
CCH_MODULE<-ga("cch_module",file.path(dirname(sub("--file=","",grep("--file=",commandArgs(FALSE),value=TRUE)[1])),"cch_module.R"))
if(!file.exists(CCH_MODULE)) CCH_MODULE<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/cch_module.R"
dir.create(OUTD,recursive=TRUE,showWarnings=FALSE)
FIPS<-c("1"="AL","2"="AK","4"="AZ","5"="AR","6"="CA","8"="CO","9"="CT","10"="DE","12"="FL","13"="GA","16"="ID","17"="IL","18"="IN","19"="IA","20"="KS","21"="KY","22"="LA","23"="ME","24"="MD","25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT","31"="NE","32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND","39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD","47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV","55"="WI","56"="WY")
cat("==== conus_eq_projector_v2 (3 fixes) ====\n"); cat(sprintf("  variant=%s mode=%s (CONFIG=%s) nstands=%s %dx%dyr\n",VARIANT,MODE,CONFIG,if(NSTAND>0)NSTAND else "ALL",NCYC,CYCLEN))
cat("  WIRED: DG(b2/b1), mortality(logit), ht-dbh(Wykoff), CR-recession, dynamic competition.\n")
cat("  FIXES: (1) identical-stand seeding standinit+treeinit, (2) ingrowth lookup, (3) SDImax self-thinning cap.\n\n")
for(f in c(PAIRS,TRAITS,V8SUM,V8META,RE_RDS,MORT_SUM,MORT_SP,HTD_SUM,HTD_SP,CR_SUM,CR_SP,IGLOOK)) if(!file.exists(f)) stop("MISSING: ",f)
if(MODE=="free") for(f in c(B1FIT,B1META)) if(!file.exists(f)) stop("MISSING b1: ",f)
## ==== v4 PATCH 3 (2026-07-28): COUNT THE SILENT "!is.finite -> benign default"
## COERCIONS. Five sites in this projector turned a NaN into "no growth" or "no
## mortality" without leaving a trace. A stand with zero growth and zero
## mortality for 100 cycles must be impossible to produce without appearing in
## the run record, so every one of them is now counted into the audit CSV.
.coerce <- new.env(parent = emptyenv())
for(.k in c("dg_a","p_die","surv_gompit","pa_gompit","survcyc_gompit","tpa"))
  assign(.k, 0, envir = .coerce)
.ctally <- function(v, tag){ n <- sum(!is.finite(v))
  if(n) assign(tag, get(tag, envir=.coerce) + n, envir=.coerce); invisible(n) }

if(MORTMODE=="gompit"){
  if(!file.exists(GOMP_COEF)) stop("MISSING gompit coefficients: ",GOMP_COEF)
  if(!file.exists(CCH_MODULE)) stop("MISSING cch_module.R: ",CCH_MODULE)
  source(CCH_MODULE)
  gcoef<-fread(GOMP_COEF)
  gomp_b0<-setNames(gcoef$b0,as.character(gcoef$SPCD)); gomp_b1<-setNames(gcoef$b1,as.character(gcoef$SPCD))
  gomp_b2<-setNames(gcoef$b2,as.character(gcoef$SPCD)); gomp_b3<-setNames(gcoef$b3,as.character(gcoef$SPCD))
  gomp_b4<-setNames(gcoef$b4,as.character(gcoef$SPCD))
  ## v2run: size-aware term b5 (coefficient of log(DBH_inches)). Older coefficient files
  ## without a b5 column fall back to 0 so size-blind coefficients behave as before.
  gomp_b5<-if("b5" %in% names(gcoef)) setNames(gcoef$b5,as.character(gcoef$SPCD)) else setNames(rep(0,nrow(gcoef)),as.character(gcoef$SPCD))
  ## genus/hardwood-softwood fallback for species without a fit: use the
  ## conifer-median or hardwood-median parameter vector (Greg fits 133 species).
  sw_med<-function(col,sw){ if(!col %in% names(gcoef)) return(0); v<-gcoef[[col]][ (gcoef$SPCD<300)==sw ]; median(v,na.rm=TRUE) }
  GFB<-list(sw=c(b0=sw_med("b0",TRUE),b1=sw_med("b1",TRUE),b2=sw_med("b2",TRUE),b3=sw_med("b3",TRUE),b4=sw_med("b4",TRUE),b5=sw_med("b5",TRUE)),
            hw=c(b0=sw_med("b0",FALSE),b1=sw_med("b1",FALSE),b2=sw_med("b2",FALSE),b3=sw_med("b3",FALSE),b4=sw_med("b4",FALSE),b5=sw_med("b5",FALSE)))
  cat(sprintf("  GOMPIT mortality (v2run size-aware): %d species banked from %s; fallback SW b0=%.2f HW b0=%.2f b5(SW)=%.4f; CCH via ORGANON port (affine A=%.3f B=%.4f)\n",
      nrow(gcoef),basename(GOMP_COEF),GFB$sw["b0"],GFB$hw["b0"],GFB$sw["b5"],CCH_A,CCH_B))
  ## per-tree gompit survival over a T-year cycle on the DEPLOYED GOMPSURV inverted-survival
  ## convention that greg_mortality_coefficients_sizeaware_v2.csv was fit to (fit_sizeaware_v2.R).
  ## DBH in INCHES: eta=b0+b1*(cr+0.01)^b2+b3*cch^b4+b5*log(max(dbh,0.1));
  ##                pa=1-exp(-exp(eta)) (ANNUAL survival, b0 is a positive survival index);
  ##                surv=pa^Tyr over the cycle.
  gomp_surv<-function(SPCD,cr,cch,dbh_in,Tyr){
    k<-as.character(SPCD); n<-length(SPCD)
    b0<-gomp_b0[k]; b1<-gomp_b1[k]; b2<-gomp_b2[k]; b3<-gomp_b3[k]; b4<-gomp_b4[k]; b5<-gomp_b5[k]
    miss<-is.na(b0); if(any(miss)){ sw<-SPCD[miss]<300
      pick<-function(p) ifelse(sw,GFB$sw[p],GFB$hw[p])
      b0[miss]<-pick("b0"); b1[miss]<-pick("b1"); b2[miss]<-pick("b2"); b3[miss]<-pick("b3"); b4[miss]<-pick("b4"); b5[miss]<-pick("b5") }
    cr<-pmin(pmax(cr,1e-4),1); cch<-pmax(cch,0); ldbh<-log(pmax(dbh_in,0.1))
    cch_term<-ifelse(cch>0,cch^b4,0)
    eta<-b0+b1*(cr+0.01)^b2+b3*cch_term+b5*ldbh; eta<-pmin(pmax(eta,-30),30)
    pa<-1-exp(-exp(eta)); .ctally(pa,"pa_gompit"); pa[!is.finite(pa)]<-1   # v4 P3: counted
    surv<-pa^Tyr; .ctally(surv,"survcyc_gompit"); surv[!is.finite(surv)]<-1 # v4 P3: counted
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

## ==== v4 PATCH 1 (2026-07-28): cspi LOG-DOMAIN GUARD ========================
## All four fitted site terms enter as log(cspi + shift). The projector
## evaluated that UNGUARDED, so any plot with cspi + shift <= 0 returned NaN and
## any plot below the fitted minimum was extrapolated down a log that dives to
## -Inf. On the projection base (conus_remeasurement_pairs_metric_cond_v2_cspiv6)
## cspi runs -2.508614 .. 10.699906; 7,778 of 17,769 plots (43.7729 percent) sit
## at cspi <= 0, outside the fitted support of the height, mortality and crown
## kernels, and 445 plots (2.5044 percent) sit at cspi <= -1 and returned NaN.
## Those NaN plots became frozen and immortal from cycle 1 in every arm via the
## downstream !is.finite coercions, while growth ran at full rate.
##
## WHY NO REFIT IS REQUIRED. 36_fit_ht_dbh.R, 34_fit_mortality.R and
## 35_fit_crown_recession.R each filtered the sample to cspi > 0 and then set
##   shift = max(1 - min(cspi over retained rows), 0.01)
## so BY CONSTRUCTION min(cspi_fit) + shift = 1 EXACTLY. Re-derived on the same
## data by derive_cspi_shifts.R (2026-07-28): all three retain
## min_cspi = 0.000110137495616 and shift = 0.999889862504, and min_cspi + shift
## prints as 1 at 17 significant digits. Therefore pmax(cspi + shift, 1) is the
## IDENTITY on the entire fitted support and a flat extrapolation below it.
## No in-support prediction changes; see the validation gate in
## validate_patch1_noop.R.
##
## 32_fit_dg_kuehne_speciesfree.R is the exception. Its CSPI_SHIFT is a
## hard-coded 1.0, it RETAINED cspi in (-1, 0], and its own transform already
## carries log(pmax(cspi + 1, 0.01)). The projector simply failed to transcribe
## that pmax. Its floor is 0.01, NOT 1.
##
## REJECTED, deliberately not implemented: raising the shift (silently
## invalidates a_cspi, m6 and r6 and needs a full refit); a regional-median
## fallback (moves off-support low-productivity plots to the middle of the
## gradient and manufactures a pass on the very site-ordering test being
## measured); asinh or log1p re-expression (correct long run but a different
## covariate, needs a refit).
## HOW THE CONSTANTS ARE BUILT, and why NOT as decimal literals. A first cut
## hard-coded shift = 0.999889862504, a 12 digit rounding of the fitted value.
## That is 3.8e-13 BELOW the true double, so min(cspi_fit) + shift came out at
## 0.99999999999961..., just under 1, the floor bit on the support's lower
## endpoint, and the validation gate correctly failed with a 2.2e-13 discrepancy
## in predicted height. A rounded constant silently destroys the exactness the
## whole no-op argument rests on.
##
## So the shift is RECONSTRUCTED with the fit's own arithmetic instead. Each fit
## computed shift = max(1.0 - min(cspi retained), 0.01). Doing that same double
## operation here on the same minimum reproduces the fitted shift bit for bit,
## with no decimal rounding anywhere. The minimum itself is stored at 17
## significant digits, which round-trips a double exactly.
##
## The FLOOR is then min_fit + shift, evaluated exactly as the kernels evaluate
## cspi + shift, rather than a hard 1.0. Because floating point addition is
## monotone, cspi >= min_fit implies cspi + shift >= min_fit + shift = floor, so
## pmax(cspi + shift, floor) is PROVABLY the identity on the whole fitted
## support, not merely identical to 12 digits. It lands within one ulp of 1.0.
CSPI_MIN_FIT_HTDBH <- 1.1013749561550323e-04   # 36_fit_ht_dbh.R, re-derived 2026-07-28
CSPI_MIN_FIT_MORT  <- 1.1013749561550323e-04   # 34_fit_mortality.R
CSPI_MIN_FIT_CR    <- 1.1013749561550323e-04   # 35_fit_crown_recession.R
CSPI_SHIFT_HTDBH  <- max(1.0 - CSPI_MIN_FIT_HTDBH, 0.01)   # the fit's own expression
CSPI_SHIFT_MORT   <- max(1.0 - CSPI_MIN_FIT_MORT,  0.01)
CSPI_SHIFT_CR     <- max(1.0 - CSPI_MIN_FIT_CR,    0.01)
CSPI_FLOOR_HTDBH  <- CSPI_MIN_FIT_HTDBH + CSPI_SHIFT_HTDBH
CSPI_FLOOR_MORT   <- CSPI_MIN_FIT_MORT  + CSPI_SHIFT_MORT
CSPI_FLOOR_CR     <- CSPI_MIN_FIT_CR    + CSPI_SHIFT_CR
CSPI_SHIFT_DG_B1  <- 1.0                # 32_fit_dg_kuehne_speciesfree.R L81, a literal
CSPI_FLOOR_DG_B1  <- 0.01               # ibid L96 log(pmax(cspi + 1, 0.01)), a literal
## Assert the identity the no-op argument depends on, at load time, in every run.
for(.nm in c("HTDBH","MORT","CR")){
  .mn <- get(paste0("CSPI_MIN_FIT_",.nm)); .sh <- get(paste0("CSPI_SHIFT_",.nm)); .fl <- get(paste0("CSPI_FLOOR_",.nm))
  if(!identical(pmax(.mn + .sh, .fl), .mn + .sh))
    stop(sprintf("v4 PATCH 1 INVARIANT VIOLATED for %s: the floor bites at the support minimum",.nm)) }
cat(sprintf("[v4 PATCH 1] floor at the support minimum verified exact; floor-1 = %.3g (htdbh), %.3g (mort), %.3g (cr)
",
  CSPI_FLOOR_HTDBH-1, CSPI_FLOOR_MORT-1, CSPI_FLOOR_CR-1))
## Per-call-site counters so the number of floored evaluations lands in the
## audit CSV. cspi is a STAND-LEVEL SCALAR broadcast across the tree list, so a
## call count and a tree-broadcast count are both recorded and are not the same
## quantity: calls counts stand-cycles, tevals counts kernel evaluations.
.cspi_guard <- new.env(parent = emptyenv())
for(.k in c("dg_b1","mort","htdbh","cr")){
  assign(paste0("calls_",.k),0,envir=.cspi_guard); assign(paste0("floored_",.k),0,envir=.cspi_guard)
  assign(paste0("tevals_",.k),0,envir=.cspi_guard); assign(paste0("ftevals_",.k),0,envir=.cspi_guard)
  assign(paste0("minarg_",.k),Inf,envir=.cspi_guard) }
ln_cspi_guarded <- function(cspi, shift, floor_at, tag, ntree = 1L){
  x <- cspi + shift
  bad <- !is.finite(x) | x < floor_at
  nb <- sum(bad); nc <- length(x); nt <- as.numeric(ntree)
  assign(paste0("calls_",tag),  get(paste0("calls_",tag), envir=.cspi_guard) + nc, envir=.cspi_guard)
  assign(paste0("tevals_",tag), get(paste0("tevals_",tag),envir=.cspi_guard) + nc*nt, envir=.cspi_guard)
  if(nb){
    assign(paste0("floored_",tag), get(paste0("floored_",tag),envir=.cspi_guard) + nb, envir=.cspi_guard)
    assign(paste0("ftevals_",tag), get(paste0("ftevals_",tag),envir=.cspi_guard) + nb*nt, envir=.cspi_guard) }
  if(any(is.finite(x))) assign(paste0("minarg_",tag),
    min(get(paste0("minarg_",tag),envir=.cspi_guard), min(x[is.finite(x)])), envir=.cspi_guard)
  log(pmax(x, floor_at)) }

## ==== v4 PATCH 2 (2026-07-28): ht_ratio BOUND DERIVED FROM THE FITTED SIGMA ==
## The ht-dbh fit is lognormal with sigma = hc("sigma") and the kernel returns
## the MEAN, so ln(ht_ratio) is approximately N(-sigma^2/2, sigma). The old
## bounds [0.33, 3.0] sit at -2.653 and +3.017 sigma: asymmetric, and outside
## the model own residual envelope. Derive the bound from sigma instead, so it
## tracks any refit automatically.
HTR_K           <- 2.0
HTR_SIGMA       <- hc("sigma")
HTR_LO          <- exp(-HTR_SIGMA^2/2 - HTR_K*HTR_SIGMA)
HTR_HI          <- exp(-HTR_SIGMA^2/2 + HTR_K*HTR_SIGMA)
## Second and worse defect: the ratio was computed at YEAR-0 diameter and
## competition state and then frozen for 100 years while the thing it multiplies
## grows. Evidence it is independent of the b1 density runaway: PN in
## conus_b2_gompit has n_eta_clamped_hi = 0 and max_eta_seen 3.81283 (kernel
## height 1.37 + exp(3.81283 + 0.07577) = 50.2126 m) yet max_ht_pre_cap_m =
## 143.0826, an implied ht_ratio of 2.8495, with the 130 m cap firing 100 times
## and ZERO eta clamping. Separately LS and SN in conus_b1 both report exactly
## 378.157301441426 = 3.0 * (1.37 + exp(4.75 + 0.07577)), both clamps saturated
## at once. So decay the calibration toward the population curve.
HTR_HALFLIFE_YR <- 30
cat(sprintf("[v4 PATCH 1] cspi shifts: dg_b1=%.17g (floor %.4g)  htdbh=%.17g  mort=%.17g  cr=%.17g\n",
  CSPI_SHIFT_DG_B1,CSPI_FLOOR_DG_B1,CSPI_SHIFT_HTDBH,CSPI_SHIFT_MORT,CSPI_SHIFT_CR))
cat(sprintf("[v4 PATCH 2] ht_ratio bound from fitted sigma=%.10g, K=%.2f -> [%.6f, %.6f]; decay half life %g yr\n",
  HTR_SIGMA,HTR_K,HTR_LO,HTR_HI,HTR_HALFLIFE_YR))

## ============================ EQUATION KERNELS ============================
eta_dg_b2<-function(tl){ln_dbh<-log(tl$dbh_in); ln_cr<-log((tl$CR+0.2)/1.2); BALsw_m<-tl$BAL_SW*0.2296; BALhw_m<-tl$BAL_HW*0.2296; ln_balsw<-log(BALsw_m+0.01)
  cb1<-tl$bgi; cb2<-pmax(tl$bgi-k1,0); cb3<-pmax(tl$bgi-k2,0); spi<-tl$sp_idx; L1i<-tl$L1_idx; L2i<-tl$L2_idx; L3i<-tl$L3_idx; FTi<-tl$FT_idx; sw_c<-tl$softwood-tl$sw_mean
  b_site<-b8("b6")+RE$z_L1_bgi[L1i]+RE$species_site_slope[spi]
  eta<-b8("b0")+RE$trait_effect[spi]+RE$z_sp[spi]+RE$z_L1[L1i]+RE$z_L2[L2i]+RE$z_L3[L3i]+RE$z_FT[FTi]+b8("b1")*ln_dbh+b8("b2")*tl$dbh_in+b8("b3")*ln_cr+b8("b4")*ln_balsw+b8("b5")*BALhw_m+b_site*cb1+b8("b9a")*cb2+b8("b9b")*cb3+b8("b7")*(tl$BA*0.2296*tl$rd_add)+b8("b8")*(BALsw_m*tl$rd_add)+b8("b11")*tl$sdi_cx+b8("b12")*(tl$bgi*tl$rd_add)+b8("b13")*(tl$bgi*ln_dbh)+b8("b14")*(tl$bgi*sw_c)+b8("b15")*(tl$bgi*ln_cr)
  pmin(pmax(eta,-30),20)}
eta_dg_b1<-function(tl,P){g<-P$g; pm<-P$pm; ln_dbh<-log(tl$dbh_in); ln_cr<-log((tl$CR+0.2)/1.2); ln_balsw<-log(tl$BAL_SW*0.2296+0.01)
  ln_csi<-ln_cspi_guarded(tl$cspi,pm$cspi_shift,CSPI_FLOOR_DG_B1,"dg_b1",length(tl$dbh_in))  # v4 P1: floor 0.01, the fit own pmax
  ba_x_rd<-(tl$BA*0.2296)*tl$rd_add; bal_x_rd<-(tl$BAL*0.2296)*tl$rd_add; trait_effect<-as.numeric(tl$Wrow %*% P$gamma)
  eta<-g("b0")+trait_effect+P$z_L1[tl$b1_L1]+P$z_L2[tl$b1_L2]+P$z_L3[tl$b1_L3]+g("b1")*ln_dbh+g("b2")*tl$dbh_in+g("b3")*ln_cr+g("b4")*ln_balsw+g("b5")*(tl$BAL_HW*0.2296)+g("b6")*ln_csi+g("b7")*ba_x_rd+g("b8")*bal_x_rd
  pmin(pmax(eta,-30),20)}
eta_mort<-function(tl){dbh<-tl$dbh_in; dbh2<-dbh^2; bal_over_ba<-ifelse(tl$BA>0,tl$BAL/tl$BA,0); sqrt_ba_rd<-sqrt(pmax((tl$BA*0.2296)*tl$rd_add,0)); ln_csi<-ln_cspi_guarded(tl$cspi,tl$mort_cspi_shift,CSPI_FLOOR_MORT,"mort",length(tl$dbh_in)); zsp<-tl$mort_zsp  # v4 P1
  mc("m0")+zsp+mc("m1")*dbh+mc("m2")*dbh2+mc("m3")*bal_over_ba+mc("m4")*tl$CR+mc("m5")*sqrt_ba_rd+mc("m6")*ln_csi}
## ==== v3 PATCH A (2026-07-25): BOUNDED HEIGHT KERNEL ========================
## ht_from_dbh() was the ONLY kernel in this script with an unclamped linear
## predictor. Every sibling is bounded: eta_dg_b2 and eta_dg_b1 to [-30, 20]
## (lines 216 and 220), the gompit eta to [-30, 30] (line 94), the logit p_die
## to [0, 0.999] (line 351). eta here carries a_bard*BAm*rd_add and
## a_blrd*BALm*rd_add with rd_add = SDI/max(SDImax, 1), which has no upper
## bound. Because a_blrd (+0.023478) exceeds abs(a_bard) (0.020929), the net
## rd_add coefficient rd_add*BAm*(a_bard + a_blrd*BALm/BAm) turns POSITIVE once
## BALm/BAm > 0.8914, so SUPPRESSED trees get an exponentially inflated height
## while dominants collapse toward breast height. Observed in
## out_conus_eq/conus_eq_cs_conus_b1_metrics.csv, stand 3184642010661:
## HT_MEAN_M 17.4 m (yr 45) -> 1.41e4 (yr 55) -> 7.44e126 m (yr 100) while
## HT_DOM_M pins at 2.5527 m. Run-wide: 366 of 3,898,818 b1_native stand-year
## rows carry HT_MEAN_M > 130 m (CS 329, CR 14, PN 12, IE 10, WC 1) and 18 rows
## carry HT_DOM_M > 130 m. No gompit row and no b2_native row diverges.
##
## BOUND, justified from the fit and the data. HT_kernel = 1.37 +
## exp(eta + 0.5*sigma^2) with sigma = 0.389281 (htdbh_wykoff_lognormal fixed
## summary), so 0.5*sigma^2 = 0.075770. HT_ETA_HI = 4.75 gives a kernel ceiling
## of 1.37 + exp(4.825770) = 126.0 m, just above the tallest tree ever measured
## (Sequoia sempervirens (D.Don) Endl., Hyperion, 115.9 m) and well above the
## tallest clean value anywhere in the v2 run (92.3 m, CS b1_gompit HT_DOM_M).
## HT_ETA_LO = -30 matches the sibling convention and returns HT -> 1.37 m
## (breast height). HT_MAX_M = 130 m caps the height AFTER the per-tree
## calibration ratio tl$ht_ratio (itself clamped to [0.33, 3.0]) is applied, so
## no path in the projector can emit a non-physical height.
##
## AUDIT. Every truncation is counted and written to <otag>_ht_clamp_log.csv
## and echoed to the job log, so silent truncation cannot recur. A non-zero
## clamp count is a POSITIVE finding: it localises the stands whose DBH and BA
## are still running away, which is a mortality problem this patch does not fix.
HT_ETA_LO <- -30; HT_ETA_HI <- 4.75; HT_MAX_M <- 130
## ==== v4 PATCH 5 (2026-08-04): INGROWTH HALT AUDIT ==========================
## add_recruits used to fail SILENTLY (see the block at add_recruits below).  A
## silent failure that is not counted cannot be audited, so every path that now
## declines to recruit is tallied here and written next to the height-clamp
## audit.  n_calls counts stand-cycles offered a recruit; fallback_used counts
## stand-cycles where the BA-dominant species had no coefficients and a modelled
## species was substituted; no_modelled_species counts stand-cycles skipped
## because NO species present was modelled; empty_tl counts stand-cycles whose
## live tree list was empty.
.ig_halt <- new.env(parent = emptyenv())
.ig_halt$n_calls <- 0L; .ig_halt$fallback_used <- 0L
.ig_halt$no_modelled_species <- 0L; .ig_halt$empty_tl <- 0L
.ig_halt$stands_fallback <- character(0); .ig_halt$stands_skipped <- character(0)
.ig_halt$examples <- character(0)
.ht_clamp <- new.env(parent = emptyenv())
.ht_clamp$n <- 0; .ht_clamp$eta_hi <- 0; .ht_clamp$eta_lo <- 0; .ht_clamp$eta_nf <- 0
.ht_clamp$ht_cap <- 0; .ht_clamp$eta_max_seen <- -Inf; .ht_clamp$ht_max_seen <- -Inf
## ==== v4 PATCH 3: DECOMPOSE THE NON-FINITE eta COUNT ========================
## The v3 tally recorded THAT eta was non-finite, not WHICH TERM did it, and did
## not separate NaN from +Inf from -Inf. Those are three different diseases:
## NaN from a log of a negative argument, +Inf from a covariate runaway
## (typically rd_add via SDI/SDImax), -Inf from a collapse. Tally each term
## separately plus an unattributed bucket. NOTE the rep_len: cspi and SDImax are
## STAND-LEVEL SCALARS broadcast across the tree list, so a naive sum over them
## would not be comparable to the tree count.
.ht_clamp$eta_nan <- 0; .ht_clamp$eta_pinf <- 0; .ht_clamp$eta_ninf <- 0
.ht_clamp$nf_lncsi <- 0; .ht_clamp$nf_rdadd <- 0; .ht_clamp$nf_trait <- 0; .ht_clamp$nf_unattr <- 0
.ht_clamp$ratio_only <- 0
.ht_tally <- function(eta, ln_csi = NULL, rd_add = NULL, te = NULL){ tryCatch({
    n <- length(eta); fin <- is.finite(eta); nf <- !fin
    .ht_clamp$n      <- .ht_clamp$n      + n
    .ht_clamp$eta_nf <- .ht_clamp$eta_nf + sum(nf)
    .ht_clamp$eta_hi <- .ht_clamp$eta_hi + sum(fin & eta > HT_ETA_HI)
    .ht_clamp$eta_lo <- .ht_clamp$eta_lo + sum(fin & eta < HT_ETA_LO)
    if (any(nf)) {
      ## NaN, +Inf and -Inf partition the non-finite set exactly.
      .ht_clamp$eta_nan  <- .ht_clamp$eta_nan  + sum(is.nan(eta))
      .ht_clamp$eta_pinf <- .ht_clamp$eta_pinf + sum(is.infinite(eta) & eta > 0)
      .ht_clamp$eta_ninf <- .ht_clamp$eta_ninf + sum(is.infinite(eta) & eta < 0)
      b_csi <- if (!is.null(ln_csi)) !is.finite(rep_len(ln_csi, n)) else rep(FALSE, n)
      b_rd  <- if (!is.null(rd_add)) !is.finite(rep_len(rd_add, n)) else rep(FALSE, n)
      b_te  <- if (!is.null(te))     !is.finite(rep_len(te,     n)) else rep(FALSE, n)
      .ht_clamp$nf_lncsi  <- .ht_clamp$nf_lncsi  + sum(nf & b_csi)
      .ht_clamp$nf_rdadd  <- .ht_clamp$nf_rdadd  + sum(nf & b_rd)
      .ht_clamp$nf_trait  <- .ht_clamp$nf_trait  + sum(nf & b_te)
      .ht_clamp$nf_unattr <- .ht_clamp$nf_unattr + sum(nf & !b_csi & !b_rd & !b_te)
    }
    if (any(fin)) .ht_clamp$eta_max_seen <- max(.ht_clamp$eta_max_seen, max(eta[fin]))
  }, error = function(e) invisible(NULL)) }
ht_from_dbh<-function(tl){BAm<-tl$BA*0.2296; BALm<-tl$BAL*0.2296; bal<-BALm; sqrt_ba<-sqrt(pmax(BAm,0)); ln_csi<-ln_cspi_guarded(tl$cspi,tl$htd_cspi_shift,CSPI_FLOOR_HTDBH,"htdbh",length(tl$dbh_in)); ba_x_rd<-BAm*tl$rd_add; bal_x_rd<-BALm*tl$rd_add  # v4 P1
  te<-as.numeric(tl$Wrow %*% hc_gamma); te[is.na(te)]<-0; zsp<-tl$htd_zsp
  eta<-hc("b0")+te+zsp+hc("a_bal")*bal+hc("a_ba")*sqrt_ba+hc("a_cspi")*ln_csi+hc("a_bard")*ba_x_rd+hc("a_blrd")*bal_x_rd+hc("b1")/(tl$dbh_in+1.0)
  .ht_tally(eta, ln_csi, tl$rd_add, te)             # v4 PATCH 3: attribute the non-finite eta to its term
  eta[is.infinite(eta) & eta > 0]<-HT_ETA_HI; eta[is.infinite(eta) & eta < 0]<-HT_ETA_LO
  eta<-pmin(pmax(eta,HT_ETA_LO),HT_ETA_HI)          # v3 PATCH A: bound as the siblings are bounded
  1.37+exp(eta+0.5*hc("sigma")^2)}
## v3 PATCH A: physical cap applied AFTER the ht_ratio calibration multiplier.
## NA heights are left as NA; the projector already tolerates and backfills them.
cap_ht<-function(h){ tryCatch({
    big <- is.finite(h) & h > HT_MAX_M; inf <- is.infinite(h) & h > 0
    if (any(is.finite(h))) .ht_clamp$ht_max_seen <- max(.ht_clamp$ht_max_seen, max(h[is.finite(h)]))
    .ht_clamp$ht_cap <- .ht_clamp$ht_cap + sum(big) + sum(inf)
  }, error = function(e) invisible(NULL))
  h[is.infinite(h) & h > 0]<-HT_MAX_M; h[is.infinite(h) & h < 0]<-1.37
  ok<-!is.na(h); h[ok]<-pmin(pmax(h[ok],1.37),HT_MAX_M); h}
cr_update<-function(tl,HT2){cr<-tl$CR; ln_cr<-log(pmax(cr,1e-4)); sqrt_ba<-sqrt(pmax(tl$BA*0.2296,0)); ln_bal_ba<-log(tl$BAL/pmax(tl$BA,1e-6)+1); cr_over_rd<-cr/pmax(tl$rd_add,1e-4); ln_csi<-ln_cspi_guarded(tl$cspi,tl$cr_cspi_shift,CSPI_FLOOR_CR,"cr",length(tl$dbh_in)); zsp<-tl$cr_zsp  # v4 P1
  eta<-rc("r0")+zsp+rc("r1")*ln_cr+rc("r2")*cr+rc("r3")*sqrt_ba+rc("r4")*ln_bal_ba+rc("r5")*cr_over_rd+rc("r6")*ln_csi
  r<-1/(1+exp(eta)); HCB1<-(1-cr)*tl$HT; maxrec<-pmax(HT2-HCB1,0); HCB2<-HCB1+r*maxrec; CR2<-1-HCB2/pmax(HT2,1e-3); pmin(pmax(CR2,0.01),0.95)}
## v4 PATCH 1 sub-item. This line previously gave the height, mortality and
## crown kernels the DIAMETER GROWTH model shift, and only when b1_pack was
## loaded, so arms that load b1_pack and arms that do not were using DIFFERENT
## singularity locations for the same three kernels. The numerical difference is
## immaterial (1.0 vs 0.999889862504) but the inconsistency is not. Each kernel
## now carries its OWN fitted shift, identically in every arm.
mort_cspi_shift<-CSPI_SHIFT_MORT; htd_cspi_shift<-CSPI_SHIFT_HTDBH; cr_cspi_shift<-CSPI_SHIFT_CR
if(!is.null(b1_pack)){
  .dgs<-as.numeric(b1_pack$pm$cspi_shift)
  if(!isTRUE(all.equal(.dgs,CSPI_SHIFT_DG_B1)))
    cat(sprintf("WARN [v4 PATCH 1] b1 meta cspi_shift=%.12g differs from CSPI_SHIFT_DG_B1=%.12g; using the meta value in eta_dg_b1.\n",.dgs,CSPI_SHIFT_DG_B1))
  cat(sprintf("[v4 PATCH 1] dg_b1 shift in use (from b1 meta): %.12g\n",.dgs)) }

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
## v2run: also return dominant (top) height HT_DOM_M and TPA-weighted mean height HT_MEAN_M (metres).
## HT_DOM_M = TPA-weighted mean HT of the largest-DBH trees whose cumulative TPA <= 100 trees ac-1 (Assmann).
stand_metrics<-function(tl){DBH<-tl$dbh_in; TPA<-tl$TPA; HT<-if(!is.null(tl$HT)) tl$HT else rep(NA_real_,length(DBH)); ok<-is.finite(DBH)&is.finite(TPA)&TPA>0; DBH<-DBH[ok]; TPA<-TPA[ok]; HT<-HT[ok]; sT<-sum(TPA); htmean<-{w<-is.finite(HT)&TPA>0; if(any(w)) sum(HT[w]*TPA[w])/sum(TPA[w]) else NA_real_}; htdom<-{w<-is.finite(HT)&is.finite(DBH); if(any(w)){o<-order(DBH[w],decreasing=TRUE); hh<-HT[w][o]; tt<-TPA[w][o]; cc<-cumsum(tt); k<-which(cc<=100); if(!length(k)) k<-1L; sum(hh[k]*tt[k])/sum(tt[k])} else NA_real_}; list(BA=if(length(DBH)) sum(pi/4*DBH^2*TPA)/144 else NA_real_, QMD=if(is.finite(sT)&&sT>1e-6) sqrt(sum(DBH^2*TPA)/sT) else NA_real_, TPH=if(is.finite(sT)) sT*2.4710538 else NA_real_, HT_DOM_M=htdom, HT_MEAN_M=htmean)}

## FIX 2 helper: add recruits to a tree list (recruits inherit stand dominant species / softwood)
add_recruits<-function(tl,n_add_tpa,sid=NA_character_){
  if(n_add_tpa<=1e-6) return(tl)
  .ig_halt$n_calls <- .ig_halt$n_calls + 1L
  ## ==== v4 PATCH 5 (2026-08-04) DEFECT 1 ==================================
  ## WAS:
  ##   if(length(tl$dbh_in)){ ba_ac<-...; dom<-tl$SPCD[which.max(ba_ac)] } else dom<-sp_levels[1]
  ##   spi<-match(dom,sp_levels); if(is.na(spi)) return(tl)
  ## Two defects in two lines.  First, which.max(ba_ac) selects the single
  ## largest tree RECORD, not the basal-area dominant SPECIES, so one big stem
  ## of a minor species outvoted the species that actually holds the stand.
  ## Second, mk_tl (see the unk block above) deliberately PRESERVES the original
  ## SPCD of trees whose species is absent from sp_levels; it remaps only
  ## sp_idx / Wrow / the z-effects for the growth kernels.  So tl$SPCD
  ## legitimately carries codes such as 500 that have no DG, mortality, ht-dbh
  ## or CR coefficients, match() returned NA, and the function returned the tree
  ## list UNCHANGED.  Ingrowth stopped for that stand-cycle with no warning, no
  ## counter and no log line.  Measured incidence over the v4 2x2 (555,734
  ## stands x 4 cells): 10.9 to 15.8 percent of stands halted in at least one
  ## cycle and 3.2 to 4.7 percent halted in every cycle, the rate differing
  ## significantly BETWEEN cells because the size hierarchy, and therefore the
  ## dominant record, evolves differently under the two mortality links.
  ##
  ## The repair, consistent with how the projector handles species elsewhere:
  ##   (a) rank species by SUMMED basal area, not by a single record;
  ##   (b) walk that ranking down to the first species that IS in sp_levels,
  ##       which is exactly the modal-in-set substitution mk_tl already performs
  ##       for the growth kernels, so a recruit always carries a species every
  ##       kernel can evaluate;
  ##   (c) warn, and count, when the fallback is used;
  ##   (d) if NO species present is modelled, skip ingrowth for this stand-cycle
  ##       but LOG it, never return silently.
  if(!length(tl$dbh_in)){
    .ig_halt$empty_tl <- .ig_halt$empty_tl + 1L
    dom<-sp_levels[1]; spi<-1L
  } else {
    ba_ac<-pi/4*tl$dbh_in^2/144*tl$TPA; ba_ac[!is.finite(ba_ac)]<-0
    ba_by_sp<-sort(tapply(ba_ac,as.character(tl$SPCD),sum),decreasing=TRUE)
    cand<-suppressWarnings(as.integer(names(ba_by_sp)))
    hit<-match(cand,sp_levels); ok<-which(!is.na(hit))
    if(!length(ok)){
      ## no species in this stand has coefficients: skip, but LOUDLY
      .ig_halt$no_modelled_species <- .ig_halt$no_modelled_species + 1L
      .ig_halt$stands_skipped <- unique(c(.ig_halt$stands_skipped,as.character(sid)))
      if(length(.ig_halt$examples) < 25L)
        .ig_halt$examples <- c(.ig_halt$examples,
          sprintf("SKIP stand=%s no modelled species among {%s}",
                  as.character(sid),paste(unique(cand),collapse="/")))
      if(.ig_halt$no_modelled_species <= 5L)
        warning(sprintf("add_recruits: stand %s has NO species with ingrowth/DG coefficients (SPCD present: %s); ingrowth SKIPPED this cycle",
                        as.character(sid),paste(unique(cand),collapse="/")),call.=FALSE)
      return(tl)
    }
    if(ok[1] != 1L){
      ## BA-dominant species is unmodelled: substitute the most abundant species
      ## that IS modelled.  This is the path that used to halt silently.
      .ig_halt$fallback_used <- .ig_halt$fallback_used + 1L
      .ig_halt$stands_fallback <- unique(c(.ig_halt$stands_fallback,as.character(sid)))
      if(length(.ig_halt$examples) < 25L)
        .ig_halt$examples <- c(.ig_halt$examples,
          sprintf("FALLBACK stand=%s dom=%s -> %s",as.character(sid),cand[1],cand[ok[1]]))
      if(.ig_halt$fallback_used <= 5L)
        warning(sprintf("add_recruits: BA-dominant SPCD %s has no ingrowth/DG coefficients; falling back to the most abundant modelled species SPCD %s (stand %s)",
                        cand[1],cand[ok[1]],as.character(sid)),call.=FALSE)
    }
    dom<-cand[ok[1]]; spi<-hit[ok[1]]
  }
  tl$SPCD<-c(tl$SPCD,dom); tl$sp_idx<-c(tl$sp_idx,spi)
  tl$dbh_in<-c(tl$dbh_in,ig_dbh_in); tl$CR<-c(tl$CR,0.6); tl$HT<-c(tl$HT,NA_real_)
  tl$TPA<-c(tl$TPA,n_add_tpa)
  swv<-sw_tree_map[as.character(dom)]; swv<-if(is.na(swv))0L else as.integer(swv)
  tl$softwood<-c(tl$softwood,swv)
  tl$mort_zsp<-c(tl$mort_zsp,{v<-mort_zsp[as.character(dom)];if(is.na(v))0 else v})
  tl$htd_zsp<-c(tl$htd_zsp,{v<-htd_zsp[as.character(dom)];if(is.na(v))0 else v})
  tl$cr_zsp<-c(tl$cr_zsp,{v<-cr_zsp[as.character(dom)];if(is.na(v))0 else v})
  if(!is.null(tl$ht_ratio)) tl$ht_ratio<-c(tl$ht_ratio,1)
  ## ==== v4 PATCH 5 (2026-08-04) DEFECT 2 ==================================
  ## add_recruits extended twelve tree-level vectors but NOT cch_last, so after
  ## any successful recruitment length(tl$cch_last) != length(tl$TPA) and the
  ## guard in emit() fell through to NA_real_.  CCH_MEAN was therefore written
  ## ONLY in cycles where recruitment had FAILED, making the column an inverted
  ## indicator of the very defect repaired above.  CCH itself is computed
  ## correctly (cch_module.R stand_cch, called below and consumed by gomp_surv),
  ## so this was a metrics-writer defect, not a modelling defect.
  ##
  ## Value appended for the recruit: the stand CCH already computed at the
  ## height of the SHORTEST live tree.  CCH is crown competition above a tree's
  ## tip, so it is monotone non-increasing in height (CCH_B = 0.0036 is positive
  ## in cch_module.R); a recruit entering at ig_dbh_in is the shortest stem in
  ## the stand, so the shortest tree's value is the nearest available point on
  ## the same profile.  It needs no extrapolation and no second O(n x 40) pass
  ## over the crown loop.  Indexing by minimum height rather than taking
  ## max(cch_last) keeps this correct even if the sign of CCH_B ever changes.
  ## This value never enters the REPORTED CCH_MEAN: the metric is now the scalar
  ## stored pre-recruit below, and tl$cch_last is overwritten wholesale at the
  ## top of the next cycle, so the append carries no modelling consequence.  It
  ## exists to restore the length invariant that emit() used to depend on.
  if(!is.null(tl$cch_last)){
    cch_ref<-{ h<-tl$HT[seq_along(tl$cch_last)]
               k<-which(is.finite(h))
               if(length(k)) tl$cch_last[k[which.min(h[k])]]
               else if(length(tl$cch_last)) suppressWarnings(max(tl$cch_last,na.rm=TRUE)) else 0 }
    if(!is.finite(cch_ref)) cch_ref<-0
    tl$cch_last<-c(tl$cch_last,cch_ref)
  }
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
    ## v4 PATCH 2: bound from the fitted sigma, and guard the RAW ratio. On a
    ## cspi <= -1 plot hp0 was NaN, so the raw ratio was NaN and
    ## pmin(pmax(NaN,...)) stayed NaN; the v3 repair only covered trees whose
    ## height was MISSING, so a measured tree on such a plot kept a NaN ratio.
    htr_raw<-tl$HT/pmax(hp0,1e-3); htr_raw[!is.finite(htr_raw)]<-1
    tl$ht_ratio<-pmin(pmax(htr_raw,HTR_LO),HTR_HI); tl$ht_ratio[miss]<-1
    ## v4 PATCH 5 (2026-08-04) DEFECT 2: stand-level CCH mean for the metrics
    ## writer.  Held OUTSIDE tl deliberately: the dead-tree filter below subsets
    ## every element of tl whose length equals the tree count, which would
    ## silently mangle a length-1 scalar in a one-tree stand.
    .cch_mean_cycle <- NA_real_
    if(MORTMODE=="gompit"){
      ## Populate PROJ_YEAR 0 as well, so CCH_MEAN is complete on every emitted
      ## row rather than only on cycles 1..NCYC.  HT has just been seeded above,
      ## so the crown geometry is evaluable here.
      cch0<-stand_cch(tl$dbh_in, tl$HT/0.3048, tl$CR, tl$TPA, grp_organon(tl$SPCD))
      tl[["cch_last"]]<-cch0
      .cch_mean_cycle <- {w<-tl$TPA; k<-is.finite(cch0)&is.finite(w)&w>0
                          if(any(k)) sum(cch0[k]*w[k])/sum(w[k]) else NA_real_}
    }
    emit<-function(cy){sm<-stand_metrics(tl); py<-cy*CYCLEN; ri<<-ri+1
      ## v4 PATCH 5 (2026-08-04) DEFECT 2: WAS a TPA-weighted mean recomputed
      ## here from tl$cch_last, guarded by length(tl$cch_last)==length(tl$TPA).
      ## That guard was FALSE on every cycle in which a recruit was added,
      ## because add_recruits grew TPA without growing cch_last, so CCH_MEAN was
      ## populated only when ingrowth had failed.  The mean is now formed where
      ## CCH is computed, pre-filter and pre-recruit, over exactly the trees and
      ## the TPA vector that CCH was evaluated on, and emit() only reads it.
      cchm<-if(MORTMODE=="gompit" && is.finite(.cch_mean_cycle)) .cch_mean_cycle else NA_real_
      rows[[ri]]<<-data.table(STAND_CN=cn,STATE=st,YEAR=inv_year+py,PROJ_YEAR=py,VARIANT=VARIANT,CONFIG=CONFIG,AGB_TONS_AC=NA_real_,BA_FT2AC=sm$BA,QMD_IN=sm$QMD,TPH=sm$TPH,HT_DOM_M=sm$HT_DOM_M,HT_MEAN_M=sm$HT_MEAN_M,TOPHT=sm$HT_DOM_M,CCH_MEAN=cchm)
      ti<<-ti+1; treerows[[ti]]<<-data.table(STAND_CN=cn,CONFIG=CONFIG,PROJ_YEAR=py,SPCD=tl$SPCD,DBH_IN=tl$dbh_in,HT_M=tl$HT,TPA=tl$TPA)}
    emit(0)
    for(cy in 1:NCYC){
      eta<-if(MODE=="dependent") eta_dg_b2(tl) else eta_dg_b1(tl,b1_pack)
      dg_a<-exp(eta+sigma8^2/2); .ctally(dg_a,"dg_a"); dg_a[!is.finite(dg_a)]<-0; dg_a<-pmin(dg_a,2.0)  # v4 P3: counted
      tl$dbh_in<-pmin(pmax(tl$dbh_in+dg_a*CYCLEN,0.1),200)
      ## ---- v4 PATCH 2: DECAY the frozen year-0 calibration ratio toward the
      ## population ht-dbh curve with a 30 yr half life, so a ratio fitted at
      ## year-0 diameter and competition state cannot multiply a 100 yr height.
      htr_w<-0.5^((cy*CYCLEN)/HTR_HALFLIFE_YR)
      htr_eff<-1+(tl$ht_ratio-1)*htr_w
      ht_kern<-ht_from_dbh(tl); ht_raw<-ht_kern*htr_eff
      ## Separate ht_ratio saturation from eta saturation: count the heights
      ## breaching the cap whose UNCALIBRATED kernel height alone was under it.
      .ht_clamp$ratio_only<-.ht_clamp$ratio_only+
        sum(is.finite(ht_kern) & ht_kern<=HT_MAX_M & (!is.finite(ht_raw) | ht_raw>HT_MAX_M))
      HT2<-cap_ht(ht_raw)                        # v3 PATCH A: bounded emission
      tl$CR<-cr_update(tl,HT2); tl$HT<-HT2
      ## ---- MORTALITY: logit (default) OR Greg gompit (--mort=gompit) ----
      rd<-tl$rd_add
      st_mult<-pmin(1+pmax(rd-0.55,0)/0.45*2.0, 8)   # FIX 3 SDImax ramp (shared by both arms)
      if(MORTMODE=="gompit"){
        ## CCH per tree this cycle via the ORGANON crown-closure port (HT m->ft),
        ## crown group = refined genus crosswalk; cch on the gompit (CCH1) scale.
        gg<-grp_organon(tl$SPCD)
        cch<-stand_cch(tl$dbh_in, tl$HT/0.3048, tl$CR, tl$TPA, gg)
        tl$cch_last<-cch
        ## v4 PATCH 5 (2026-08-04) DEFECT 2: form the reported CCH mean HERE,
        ## pre-filter and pre-recruit, over the same trees and the same TPA
        ## vector stand_cch was just evaluated on.  tl$TPA has not yet been
        ## multiplied by surv at this point, so cch and the weights refer to an
        ## identical stand state.
        .cch_mean_cycle <- {w<-tl$TPA; k<-is.finite(cch)&is.finite(w)&w>0
                            if(any(k)) sum(cch[k]*w[k])/sum(w[k]) else NA_real_}
        surv_cyc<-gomp_surv(tl$SPCD, tl$CR, cch, tl$dbh_in, CYCLEN)   # v2run: DBH-in-inches size term; cycle-length survival
        ## fold the shared SDImax self-thinning ramp into the hazard for parity
        H_cyc<- -log(pmin(pmax(surv_cyc,1e-9),1))                # cycle hazard*T
        surv<-exp(-H_cyc*st_mult); .ctally(surv,"surv_gompit"); surv[!is.finite(surv)]<-1  # v4 P3: counted
      } else {
        p_die_a<-1/(1+exp(-eta_mort(tl))); .ctally(p_die_a,"p_die"); p_die_a[!is.finite(p_die_a)]<-0  # v4 P3: counted
        p_die_a<-pmin(pmax(p_die_a*st_mult,0),0.999)
        surv<-(1-p_die_a)^CYCLEN
      }
      tl$TPA<-tl$TPA*surv; .ctally(tl$TPA,"tpa"); tl$TPA[!is.finite(tl$TPA)]<-0  # v4 P3: counted
      ## drop dead-empty trees (TPA collapsed)
      keep<-tl$TPA>1e-4
      if(any(!keep)){ vn<-length(tl$dbh_in); for(nm in names(tl)) if(length(tl[[nm]])==vn) tl[[nm]]<-tl[[nm]][keep]; tl$Wrow<-tl$Wrow[keep,,drop=FALSE] }
      ## FIX 2: add ingrowth recruits for this cycle (per-acre rate x cycle length)
      tl<-add_recruits(tl, IG_TPA*CYCLEN, cn)
      tl$HT[!is.finite(tl$HT)]<-ht_from_dbh(tl)[!is.finite(tl$HT)]
      tl<-recompute_comp(tl); emit(cy)}
  }
}
out<-rbindlist(rows[seq_len(ri)]); tre<-rbindlist(treerows[seq_len(ti)]); otag<-sprintf("conus_eq_%s_%s",tolower(VARIANT),CONFIG)
fwrite(out,file.path(OUTD,paste0(otag,"_metrics.csv"))); fwrite(tre,file.path(OUTD,paste0(otag,"_treelists.csv")))
## ---- v3 PATCH A: height-clamp audit (silent truncation must be impossible) --
## ==== v4 PATCH 4 (2026-07-28): SPLIT frac_eta_at_upper_bound ================
## v3 computed (eta_hi + eta_nf)/n in a field whose NAME asserts it measures
## clamping. Verified consequence: it overstated the true clamp fraction by up
## to 4,994x and INVERTED the region ranking. As reported, LS looked worst at
## 23.36 percent and SN mid-pack at 1.00 percent; corrected, SN is worst at
## 0.120 percent, LS is fourth at 0.036 percent, and seven regions carrying
## reported fractions of 0.54 to 3.93 percent never clamped at all.
## frac_eta_at_upper_bound keeps its POSITION and its name but its numerator is
## now eta_hi ONLY. Non-finite and lower-bound get their own fields.
##
## SCHEMA CONTRACT. bakuzis_v3_chain.sbatch built the clamp roll-up
## POSITIONALLY, taking the header from head -1 of whichever file ls returned
## first and then blindly concatenating every other body with tail -n +2, so
## mixed vintages misaligned SILENTLY. Every new column is therefore APPENDED AT
## THE END, never inserted in the middle, and the shell roll-up is replaced by
## rollup_ht_clamp.R, which is schema aware (rbindlist fill=TRUE,
## use.names=TRUE) and recomputes BOTH fractions from the raw counts, so
## existing v3 logs are repaired retrospectively without re-running anything.
HT_CLAMP_SCHEMA_VERSION <- 2L
.cg <- function(p,t) get(paste0(p,"_",t), envir=.cspi_guard)
ht_audit<-tryCatch(data.table(otag=otag, VARIANT=VARIANT, CONFIG=CONFIG,
    run_date=as.character(Sys.Date()), ht_eta_lo=HT_ETA_LO, ht_eta_hi=HT_ETA_HI, ht_max_m=HT_MAX_M,
    n_tree_kernel_evals=.ht_clamp$n, n_eta_clamped_hi=.ht_clamp$eta_hi,
    n_eta_clamped_lo=.ht_clamp$eta_lo, n_eta_nonfinite=.ht_clamp$eta_nf, n_ht_capped=.ht_clamp$ht_cap,
    ## v4 PATCH 4: numerator is eta_hi ONLY (v3 wrongly added eta_nf here)
    frac_eta_at_upper_bound=if(.ht_clamp$n>0) .ht_clamp$eta_hi/.ht_clamp$n else NA_real_,
    max_eta_seen=.ht_clamp$eta_max_seen, max_ht_pre_cap_m=.ht_clamp$ht_max_seen,
    ## ---------------- v4 columns, APPENDED AT THE END ONLY ----------------
    schema_version=HT_CLAMP_SCHEMA_VERSION,
    frac_eta_nonfinite=if(.ht_clamp$n>0) .ht_clamp$eta_nf/.ht_clamp$n else NA_real_,
    frac_eta_at_lower_bound=if(.ht_clamp$n>0) .ht_clamp$eta_lo/.ht_clamp$n else NA_real_,
    n_eta_nan=.ht_clamp$eta_nan, n_eta_posinf=.ht_clamp$eta_pinf, n_eta_neginf=.ht_clamp$eta_ninf,
    n_eta_nf_lncsi=.ht_clamp$nf_lncsi, n_eta_nf_rdadd=.ht_clamp$nf_rdadd,
    n_eta_nf_trait=.ht_clamp$nf_trait, n_eta_nf_unattr=.ht_clamp$nf_unattr,
    htr_sigma=HTR_SIGMA, htr_k=HTR_K, htr_lo=HTR_LO, htr_hi=HTR_HI,
    htr_halflife_yr=HTR_HALFLIFE_YR, n_ht_capped_by_ratio_only=.ht_clamp$ratio_only,
    cspi_shift_dg_b1=CSPI_SHIFT_DG_B1, cspi_floor_dg_b1=CSPI_FLOOR_DG_B1,
    cspi_shift_htdbh=CSPI_SHIFT_HTDBH, cspi_shift_mort=CSPI_SHIFT_MORT,
    cspi_shift_cr=CSPI_SHIFT_CR,
    cspi_floor_htdbh=CSPI_FLOOR_HTDBH, cspi_floor_mort=CSPI_FLOOR_MORT, cspi_floor_cr=CSPI_FLOOR_CR,
    cspi_min_fit_htdbh=CSPI_MIN_FIT_HTDBH, cspi_min_fit_mort=CSPI_MIN_FIT_MORT, cspi_min_fit_cr=CSPI_MIN_FIT_CR,
    n_cspi_calls_dg_b1=.cg("calls","dg_b1"), n_cspi_floored_dg_b1=.cg("floored","dg_b1"),
    n_cspi_tevals_dg_b1=.cg("tevals","dg_b1"), n_cspi_floored_tevals_dg_b1=.cg("ftevals","dg_b1"),
    min_cspi_arg_dg_b1=.cg("minarg","dg_b1"),
    n_cspi_calls_mort=.cg("calls","mort"), n_cspi_floored_mort=.cg("floored","mort"),
    n_cspi_tevals_mort=.cg("tevals","mort"), n_cspi_floored_tevals_mort=.cg("ftevals","mort"),
    min_cspi_arg_mort=.cg("minarg","mort"),
    n_cspi_calls_htdbh=.cg("calls","htdbh"), n_cspi_floored_htdbh=.cg("floored","htdbh"),
    n_cspi_tevals_htdbh=.cg("tevals","htdbh"), n_cspi_floored_tevals_htdbh=.cg("ftevals","htdbh"),
    min_cspi_arg_htdbh=.cg("minarg","htdbh"),
    n_cspi_calls_cr=.cg("calls","cr"), n_cspi_floored_cr=.cg("floored","cr"),
    n_cspi_tevals_cr=.cg("tevals","cr"), n_cspi_floored_tevals_cr=.cg("ftevals","cr"),
    min_cspi_arg_cr=.cg("minarg","cr"),
    n_coerce_dg_a=get("dg_a",envir=.coerce), n_coerce_p_die=get("p_die",envir=.coerce),
    n_coerce_surv_gompit=get("surv_gompit",envir=.coerce), n_coerce_pa_gompit=get("pa_gompit",envir=.coerce),
    n_coerce_survcyc_gompit=get("survcyc_gompit",envir=.coerce), n_coerce_tpa=get("tpa",envir=.coerce)),
  error=function(e){ cat("WARN [PATCH A] ht_audit assembly failed:",conditionMessage(e),"\n"); NULL })
if(!is.null(ht_audit)){
  tryCatch(fwrite(ht_audit,file.path(OUTD,paste0(otag,"_ht_clamp_log.csv"))),
    error=function(e) cat("WARN [PATCH A] could not write ht_clamp_log:",conditionMessage(e),"\n"))
  cat(sprintf("\n[PATCH A] height clamp: eta>%g on %d, eta<%g on %d, eta non-finite on %d of %d tree-kernel evaluations; HT capped at %g m on %d trees; max eta seen %.6g; max pre-cap HT %.6g m\n",
      HT_ETA_HI,.ht_clamp$eta_hi,HT_ETA_LO,.ht_clamp$eta_lo,.ht_clamp$eta_nf,.ht_clamp$n,
      HT_MAX_M,.ht_clamp$ht_cap,.ht_clamp$eta_max_seen,.ht_clamp$ht_max_seen))
  cat(sprintf("[v4 PATCH 4] frac_eta_at_upper_bound=%.6g (eta_hi ONLY)  frac_eta_nonfinite=%.6g  frac_eta_at_lower_bound=%.6g  schema_version=%d\n",
      if(.ht_clamp$n>0) .ht_clamp$eta_hi/.ht_clamp$n else NA_real_,
      if(.ht_clamp$n>0) .ht_clamp$eta_nf/.ht_clamp$n else NA_real_,
      if(.ht_clamp$n>0) .ht_clamp$eta_lo/.ht_clamp$n else NA_real_, HT_CLAMP_SCHEMA_VERSION))
  cat(sprintf("[v4 PATCH 3] non-finite eta decomposition: NaN=%d +Inf=%d -Inf=%d | ln_csi=%d rd_add=%d trait=%d unattributed=%d\n",
      .ht_clamp$eta_nan,.ht_clamp$eta_pinf,.ht_clamp$eta_ninf,
      .ht_clamp$nf_lncsi,.ht_clamp$nf_rdadd,.ht_clamp$nf_trait,.ht_clamp$nf_unattr))
  cat(sprintf("[v4 PATCH 3] benign-default coercions: dg_a=%g p_die=%g surv_gompit=%g pa_gompit=%g survcyc_gompit=%g tpa=%g\n",
      get("dg_a",envir=.coerce),get("p_die",envir=.coerce),get("surv_gompit",envir=.coerce),
      get("pa_gompit",envir=.coerce),get("survcyc_gompit",envir=.coerce),get("tpa",envir=.coerce)))
  cat(sprintf("[v4 PATCH 2] ht cap reached with the UNCALIBRATED kernel height already under the cap on %d trees (ht_ratio saturation, not eta saturation)\n",
      .ht_clamp$ratio_only))
  cat(sprintf("[v4 PATCH 1] cspi floored evaluations: dg_b1=%g/%g mort=%g/%g htdbh=%g/%g cr=%g/%g (calls); min arg seen dg_b1=%.6g mort=%.6g htdbh=%.6g cr=%.6g\n",
      .cg("floored","dg_b1"),.cg("calls","dg_b1"),.cg("floored","mort"),.cg("calls","mort"),
      .cg("floored","htdbh"),.cg("calls","htdbh"),.cg("floored","cr"),.cg("calls","cr"),
      .cg("minarg","dg_b1"),.cg("minarg","mort"),.cg("minarg","htdbh"),.cg("minarg","cr")))
  if(.ht_clamp$eta_hi+.ht_clamp$eta_nf>0) cat("[PATCH A] NOTE: the upper bound was REACHED. Those stands are still growing without bound in DBH and BA; this patch bounds the HEIGHT EMISSION only, it does not fix the underlying increment-mortality imbalance. Inspect BA_FT2AC in the affected stands.\n")
}
## ==== v4 PATCH 5 (2026-08-04): INGROWTH AUDIT OUT ===========================
## The whole point of the repair is that a declined recruitment can no longer be
## invisible, so it is written to a file, not only to stdout.
ig_audit<-tryCatch(data.table(
    VARIANT=VARIANT, CONFIG=CONFIG, MORTMODE=MORTMODE,
    IG_TPA=IG_TPA, IG_BA=IG_BA, ig_dbh_in=ig_dbh_in,
    n_recruit_calls=.ig_halt$n_calls,
    n_fallback_used=.ig_halt$fallback_used,
    n_no_modelled_species=.ig_halt$no_modelled_species,
    n_empty_treelist=.ig_halt$empty_tl,
    n_stands_fallback=length(.ig_halt$stands_fallback),
    n_stands_skipped=length(.ig_halt$stands_skipped),
    frac_calls_fallback=if(.ig_halt$n_calls>0) .ig_halt$fallback_used/.ig_halt$n_calls else NA_real_,
    frac_calls_skipped=if(.ig_halt$n_calls>0) .ig_halt$no_modelled_species/.ig_halt$n_calls else NA_real_,
    examples=paste(.ig_halt$examples,collapse=" | ")),
  error=function(e){cat("WARN [PATCH 5] ingrowth audit assembly failed:",conditionMessage(e),"\n"); NULL})
if(!is.null(ig_audit)){
  tryCatch(fwrite(ig_audit,file.path(OUTD,paste0(otag,"_ingrowth_audit.csv"))),
    error=function(e) cat("WARN [PATCH 5] could not write ingrowth_audit:",conditionMessage(e),"\n"))
  cat(sprintf("\n[v4 PATCH 5] ingrowth: %d recruit calls; BA-dominant species unmodelled and SUBSTITUTED on %d (%.4g of calls, %d stands); NO modelled species so SKIPPED on %d (%.4g of calls, %d stands); empty tree list on %d\n",
    .ig_halt$n_calls,.ig_halt$fallback_used,
    if(.ig_halt$n_calls>0) .ig_halt$fallback_used/.ig_halt$n_calls else NA_real_,
    length(.ig_halt$stands_fallback),.ig_halt$no_modelled_species,
    if(.ig_halt$n_calls>0) .ig_halt$no_modelled_species/.ig_halt$n_calls else NA_real_,
    length(.ig_halt$stands_skipped),.ig_halt$empty_tl))
  if(length(.ig_halt$examples)) cat("[v4 PATCH 5] examples: ",paste(head(.ig_halt$examples,10),collapse=" | "),"\n",sep="")
  if(.ig_halt$fallback_used==0 && .ig_halt$no_modelled_species==0)
    cat("[v4 PATCH 5] no stand-cycle declined a recruit; pre-patch this run would have been identical.\n")
}
cat("\n  stands projected:",nproj,"\n")
cat("Wrote:",file.path(OUTD,paste0(otag,"_metrics.csv")),"(",nrow(out),"rows )\n"); cat("Wrote:",file.path(OUTD,paste0(otag,"_treelists.csv")),"(",nrow(tre),"rows )\n")
cat("\n=== year 0/50/100 stand-mean metrics (",CONFIG,") ===\n")
print(out[PROJ_YEAR %in% c(0,50,100),.(BA=mean(BA_FT2AC,na.rm=TRUE),QMD=mean(QMD_IN,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE),HT_DOM_M=mean(HT_DOM_M,na.rm=TRUE),TOPHT=mean(TOPHT,na.rm=TRUE),n=.N),by=PROJ_YEAR])
cat("DONE.\n")
