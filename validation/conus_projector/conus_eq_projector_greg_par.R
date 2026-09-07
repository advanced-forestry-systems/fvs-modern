#!/usr/bin/env Rscript
## conus_eq_projector_greg_par.R -- PARALLEL (mclapply) copy of faithful Greg arm. IDENTICAL eqs/params/forms; only per-stand loop parallelized.
## ----------------------------------------------------------------------------------------
## Uses Greg Johnson's ACTUAL fitted parameters and ACTUAL equation FORMS:
##   DG  : ~/fvs_remodeling/rds/dg_parms.RDS  (84 spp; spcd,B0..B6) + est_dg annual loop
##         dg = exp( B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7)
##                    + B5*elev + B6*EMT )                       [code form is authoritative]
##   HG  : ~/fvs_remodeling/rds/hg_parms.RDS  (96 spp; spcd,B0=max_height,B1..B8) + est_hg loop
##         dht = max_height*b1*b2*cr^b3 * exp(-b1*ht -b4*ccfl -b8*cch^0.5 -b5*elev
##                    +b6*TD^0.5 +b7*EMT) * (1-exp(-b1*ht))^(b2-1)
##   MORT: ~/fvs_remodeling/rds/mort_parm_base_rate_cr_cch.RDS (92 spp; SPCD,b0..b4) + gompit
##         P_surv = 1 - exp(-exp( b0 + b1*(cr+0.01)^b2 + b3*cch^b4 ))   (annual loop)
##   CROWN CHANGE: Greg's repo has NO fitted crown-change prediction equation (the
##         Crown_Change_Equations_for_CONUS.qmd "Proposed Model" section restates the
##         mortality gompit and fits mort_parm_base_rate_cr_cch2.RDS; only an exploratory
##         data subset exists, no estimated dcl model). To keep the arm RUNNABLE without
##         fabricating a Greg crown equation, crown ratio is updated with the SAME
##         fvs-conus CR-recession kernel the other eq-projection arms use, applied
##         identically. This is the ONE non-Greg component and is flagged in the log/report.
##
## EMT/TD: per-STAND_CN lookup (greg_emt_td_lookup.rds) extracted from Greg's ClimateNA
##         1991-2020 normals (EMT.tif; TD=MWMT-MCMT) at standinit public/fuzzed LAT/LON.
## ELEV  : standinit ELEVFT (feet), as Greg's est_dg/est_hg expect (elev in feet).
##
## Seeding/competition/ingrowth/output mirror conus_eq_projector_v2.R EXACTLY (identical
## standinit + treeinit STAND_CN passthrough, recompute_comp each cycle, ingrowth lookup,
## SDImax self-thinning ramp folded into hazard). cch via the existing cch_module.R
## (ORGANON crown-closure port) -- same as the v2 --mort=gompit branch.
##
## Species without Greg DG/HG/MORT params -> fvs-conus species-free (trait-driven) fallback:
##   DG/HG  : remap to stand's modal in-Greg-set species for the growth kernels (parameters
##            borrowed) -- LOGGED as fallback. (Faithful Greg has no params for them.)
##   MORT   : softwood/hardwood median Greg parameter vector.
## Per-stand fallback fraction (tree-count weighted) is logged and emitted.
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){ m<-grep(paste0("^--",n,"="),args,value=TRUE); if(!length(m)) return(d); sub(paste0("^--",n,"="),"",m[1]) }
ROOT<-"/users/PUOM0008/crsfaaron/fvs-conus"; SCR<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
PROJ<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
OUTD<-ga("outdir",file.path(PROJ,"out_conus_eq_greg")); VARIANT<-toupper(ga("variant","NE"))
NSTAND<-as.integer(ga("nstands","0")); SEED<-as.integer(ga("seed","7")); NCYC<-as.integer(ga("ncycles","20")); CYCLEN<-as.integer(ga("cyclelen","5"))
STANDINIT_DIR<-ga("standinit_dir",file.path(SCR,"standinit_by_variant"))
TREEINIT_DIR<-ga("treeinit_dir","/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h")
PAIRS<-ga("pairs",file.path(ROOT,"data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"))
TRAITS<-ga("traits",file.path(ROOT,"traits/species_traits.rds"))
IGLOOK<-ga("ingrowth",file.path(ROOT,"output/comparisons_overstory/intermediate/ingrowth_lookup.rds"))
GREG_RDS<-ga("greg_rds","/users/PUOM0008/crsfaaron/fvs_remodeling/rds")
EMT_LOOK<-ga("emt_td",file.path(PROJ,"greg_emt_td_lookup.rds"))
CCH_MODULE<-ga("cch_module",file.path(PROJ,"cch_module.R"))
## fvs-conus CR-recession (the one non-Greg piece, used to update crown ratio)
CC<-file.path(ROOT,"output/conus")
CR_SUM<-file.path(CC,"cr_recession_cspi_traits1_fixed_summary.csv"); CR_SP<-file.path(CC,"cr_recession_cspi_traits1_species_intercepts.csv")
CONFIG<-"greg"
dir.create(OUTD,recursive=TRUE,showWarnings=FALSE)
FIPS<-c("1"="AL","2"="AK","4"="AZ","5"="AR","6"="CA","8"="CO","9"="CT","10"="DE","12"="FL","13"="GA","16"="ID","17"="IL","18"="IN","19"="IA","20"="KS","21"="KY","22"="LA","23"="ME","24"="MD","25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT","31"="NE","32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND","39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD","47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV","55"="WI","56"="WY")
cat("==== conus_eq_projector_GREG (faithful Greg fvs_remodeling arm) ====\n")
cat(sprintf("  variant=%s CONFIG=%s nstands=%s %dx%dyr\n",VARIANT,CONFIG,if(NSTAND>0)NSTAND else "ALL",NCYC,CYCLEN))
cat("  Greg DG(B0..B6,elev,EMT) + Greg HG(max_ht,b1..b8,ccfl,cch,elev,TD,EMT) + Greg gompit mort(b0..b4,cr,cch).\n")
cat("  CR update = fvs-conus CR-recession (Greg has no fitted crown-change eq -- FLAGGED).\n\n")
for(f in c(PAIRS,TRAITS,IGLOOK,EMT_LOOK,CCH_MODULE,CR_SUM,CR_SP)) if(!file.exists(f)) stop("MISSING: ",f)
source(CCH_MODULE)

## ----------------------------- GREG PARAMS -----------------------------
dgP<-as.data.table(readRDS(file.path(GREG_RDS,"dg_parms.RDS")))
hgP<-as.data.table(readRDS(file.path(GREG_RDS,"hg_parms.RDS")))
moP<-as.data.table(readRDS(file.path(GREG_RDS,"mort_parm_base_rate_cr_cch.RDS")))
## dedupe: DG/HG one row per spcd (converged preferred); MORT lowest-nll per SPCD
dgP<-dgP[order(spcd,-isConv)][!duplicated(spcd)]
hgP<-hgP[isConv==TRUE][order(spcd)][!duplicated(spcd)]
moP<-moP[order(SPCD,nll)][!duplicated(SPCD)]
dg_set<-dgP$spcd; hg_set<-hgP$spcd; mo_set<-moP$SPCD
DG<-as.list(setNames(lapply(c("B0","B1","B2","B3","B4","B5","B6"),function(c) setNames(dgP[[c]],as.character(dgP$spcd))),c("B0","B1","B2","B3","B4","B5","B6")))
HG<-as.list(setNames(lapply(c("B0","B1","B2","B3","B4","B5","B6","B7","B8"),function(c) setNames(hgP[[c]],as.character(hgP$spcd))),c("B0","B1","B2","B3","B4","B5","B6","B7","B8")))
MO<-as.list(setNames(lapply(c("b0","b1","b2","b3","b4"),function(c) setNames(moP[[c]],as.character(moP$SPCD))),c("b0","b1","b2","b3","b4")))
## softwood / hardwood median mortality fallback (SPCD<300 = conifer convention)
mo_sw<-moP[SPCD<300]; mo_hw<-moP[SPCD>=300]
mo_fb<-list(sw=sapply(c("b0","b1","b2","b3","b4"),function(c) median(mo_sw[[c]],na.rm=TRUE)),
            hw=sapply(c("b0","b1","b2","b3","b4"),function(c) median(mo_hw[[c]],na.rm=TRUE)))
cat(sprintf("  Greg params: DG %d spp, HG %d spp, MORT %d spp. Mort fallback SW b0=%.2f HW b0=%.2f\n",
  length(dg_set),length(hg_set),length(mo_set),mo_fb$sw["b0"],mo_fb$hw["b0"]))

## EMT/TD per STAND_CN
emt_look<-as.data.table(readRDS(EMT_LOOK)); setkey(emt_look,STAND_CN)
EMTv<-setNames(emt_look$EMT,emt_look$STAND_CN); TDv<-setNames(emt_look$TD,emt_look$STAND_CN)
EMT_MED<-median(emt_look$EMT,na.rm=TRUE); TD_MED<-median(emt_look$TD,na.rm=TRUE)

## traits + softwood map (for fallback species selection & mort sw/hw)
traits<-as.data.table(readRDS(TRAITS))
sw_tree_map<-setNames(traits$softwood,as.character(traits$SPCD))

## fvs-conus CR-recession (non-Greg; crown update only)
read_fx<-function(p){s<-fread(p); setNames(s$mean,s$variable)}
fxc<-read_fx(CR_SUM); rc<-function(n) as.numeric(fxc[[n]]); cr_sp_dt<-fread(CR_SP); cr_zsp<-setNames(cr_sp_dt$mean,as.character(cr_sp_dt$SPCD))
## cspi shift placeholder (CR-recession used ln(cspi+shift)); use 1.0 (same default as v2 when b1 absent)
CR_CSPI_SHIFT<-1.0

## ----------------------------- INGROWTH (FIX 2) -----------------------------
ig_lookup<-readRDS(IGLOOK); ig<-ig_lookup[[VARIANT]]; if(is.null(ig)) ig<-ig_lookup[["OVERALL"]]
IG_TPA<-if(!is.null(ig)) as.numeric(ig$med_ann_TPA) else 0
IG_BA <-if(!is.null(ig)) as.numeric(ig$med_ann_BA)  else 0
ig_dbh_in<-if(IG_TPA>1e-6 && IG_BA>0){ ba_tree<-IG_BA/IG_TPA; sqrt(ba_tree/(pi/4)*144) } else 1.5
ig_dbh_in<-min(max(ig_dbh_in,1.0),3.0)
cat(sprintf("  FIX2 ingrowth (%s): TPA=%.3f rec/ac/yr BA=%.3f ft2/ac/yr; recruit DBH=%.2f in\n",VARIANT,IG_TPA,IG_BA,ig_dbh_in))

## ============================ FIX 1: SEEDING ============================
SI_FILE<-file.path(STANDINIT_DIR,paste0("standinit_",VARIANT,".csv")); if(!file.exists(SI_FILE)) stop("MISSING standinit: ",SI_FILE)
cat("STEP 1: load standinit + treeinit (identical stands) + join covariates ...\n")
si<-fread(SI_FILE,colClasses=list(character=c("STAND_CN","STAND_ID")))
si[,STAND_CN:=sub("\\..*$","",STAND_CN)]
si[,INV_YEAR:=suppressWarnings(as.integer(INV_YEAR))]; si[is.na(INV_YEAR),INV_YEAR:=2010L]
si[,STATE:=suppressWarnings(as.integer(STATE))]
si[,ELEVFT:=suppressWarnings(as.numeric(ELEVFT))]
si<-si[is.finite(STATE)&!is.na(STAND_CN)&STAND_CN!=""]
## site covars from pairs (SDImax only needed for self-thinning cap; EPA not needed for Greg kernels)
d<-as.data.table(readRDS(PAIRS))
pl<-d[,.(SDImax_brms=first(SDImax_brms),LAT=first(LAT),LON=first(LON)),by=.(STATECD,COUNTYCD,PLOT)]
cty<-pl[,.(SDImax_brms=median(SDImax_brms,na.rm=TRUE)),by=.(STATECD,COUNTYCD)]
sta<-pl[,.(SDImax_brms=median(SDImax_brms,na.rm=TRUE)),by=.(STATECD)]
si[,sid_county:=suppressWarnings(as.integer(substr(STAND_ID,5,7)))]
si[,SDImax_brms:=NA_real_]
plf<-pl[is.finite(LAT)&is.finite(LON)]
for(st in unique(si$STATE)){
  idx<-which(si$STATE==st & is.finite(si$LATITUDE) & is.finite(si$LONGITUDE)); if(!length(idx)) next
  cand<-plf[STATECD==st]; if(!nrow(cand)) next
  la<-si$LATITUDE[idx]; lo<-si$LONGITUDE[idx]
  nn<-sapply(seq_along(idx),function(i){dx<-cand$LAT-la[i];dy<-cand$LON-lo[i]; which.min(dx*dx+dy*dy)})
  set(si,i=idx,j="SDImax_brms",value=cand$SDImax_brms[nn])
}
si[,kc:=paste(STATE,sid_county,sep="_")]; cty[,kc:=paste(STATECD,COUNTYCD,sep="_")]
mc_<-cty[match(si$kc,cty$kc)]; nas<-is.na(si$SDImax_brms); if(any(nas)) set(si,i=which(nas),j="SDImax_brms",value=mc_$SDImax_brms[which(nas)])
ms_<-sta[match(si$STATE,sta$STATECD)]; nas<-is.na(si$SDImax_brms); if(any(nas)) set(si,i=which(nas),j="SDImax_brms",value=ms_$SDImax_brms[which(nas)])
si[is.na(SDImax_brms),SDImax_brms:=600]
## attach EMT/TD per stand
si[,EMT:=EMTv[STAND_CN]]; si[,TD:=TDv[STAND_CN]]
si[!is.finite(EMT),EMT:=EMT_MED]; si[!is.finite(TD),TD:=TD_MED]
si[!is.finite(ELEVFT)|ELEVFT<0,ELEVFT:=median(si$ELEVFT[is.finite(si$ELEVFT)&si$ELEVFT>=0],na.rm=TRUE)]
cat(sprintf("  standinit stands: %d  EMT/TD coverage=%.1f%%  ELEV median=%.0f ft  EMT med=%.1f TD med=%.1f\n",
  nrow(si),100*mean(is.finite(EMTv[si$STAND_CN])),median(si$ELEVFT,na.rm=TRUE),median(si$EMT),median(si$TD)))

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

## ============================ GREG EQUATION KERNELS ============================
## DG annual increment (in/yr) at current state. dbh in, cr frac, ht ft, bal ft2/ac, elev ft, emt degC.
dg_annual<-function(SPCD,dbh,cr,ht,bal,elev,emt){
  k<-as.character(SPCD)
  B0<-DG$B0[k];B1<-DG$B1[k];B2<-DG$B2[k];B3<-DG$B3[k];B4<-DG$B4[k];B5<-DG$B5[k];B6<-DG$B6[k]
  z<- B0 + B1*log((dbh+1)^2/(cr*ht+1.0)^B3) + B2*bal^B4/log(dbh+2.7) + B5*elev + B6*emt
  z<-pmin(pmax(z,-30),5); g<-exp(z); g[!is.finite(g)]<-0; pmax(g,0)
}
## HG annual increment (ft/yr). max_height=B0; ht ft, cr frac, ccfl ft2/ac, cch frac (0-1), elev ft, td, emt.
hg_annual<-function(SPCD,ht,cr,ccfl,cch,elev,td,emt){
  k<-as.character(SPCD)
  mx<-HG$B0[k];b1<-HG$B1[k];b2<-HG$B2[k];b3<-HG$B3[k];b4<-HG$B4[k];b5<-HG$B5[k];b6<-HG$B6[k];b7<-HG$B7[k];b8<-HG$B8[k]
  crp<-pmax(cr,1e-4); cchp<-pmax(cch,0)
  dht<- mx*b1*b2*crp^b3*exp(-b1*ht -b4*ccfl -b8*cchp^0.5 -b5*elev +b6*sqrt(pmax(td,0)) +b7*emt)*(1.0-exp(-b1*ht))^(b2-1.0)
  dht[!is.finite(dht)]<-0; pmax(dht,0)
}
## Greg gompit annual survival probability. cr frac, cch frac.
gomp_surv_annual<-function(SPCD,cr,cch){
  k<-as.character(SPCD); n<-length(SPCD)
  b0<-MO$b0[k];b1<-MO$b1[k];b2<-MO$b2[k];b3<-MO$b3[k];b4<-MO$b4[k]
  miss<-is.na(b0)
  if(any(miss)){ sw<-SPCD[miss]<300
    b0[miss]<-ifelse(sw,mo_fb$sw["b0"],mo_fb$hw["b0"]); b1[miss]<-ifelse(sw,mo_fb$sw["b1"],mo_fb$hw["b1"])
    b2[miss]<-ifelse(sw,mo_fb$sw["b2"],mo_fb$hw["b2"]); b3[miss]<-ifelse(sw,mo_fb$sw["b3"],mo_fb$hw["b3"]); b4[miss]<-ifelse(sw,mo_fb$sw["b4"],mo_fb$hw["b4"]) }
  crp<-pmax(cr,1e-4); cchp<-pmax(cch,0)
  eta<- b0 + b1*(crp+0.01)^b2 + b3*ifelse(cchp>0,cchp^b4,0)
  eta<-pmin(pmax(eta,-30),30)
  ps<-1-exp(-exp(eta)); ps[!is.finite(ps)]<-1; ps[cr<=0]<-0
  pmin(pmax(ps,0),1)
}
## CR update (fvs-conus CR-recession; NON-GREG, flagged). HT2 in metres.
cr_update<-function(SPCD,cr,HT,HT2,BA_m,BAL,BA_ft,cspi,zsp){
  ln_cr<-log(pmax(cr,1e-4)); sqrt_ba<-sqrt(pmax(BA_m,0)); ln_bal_ba<-log(BAL/pmax(BA_ft,1e-6)+1)
  rd_add_proxy<-1.0; cr_over_rd<-cr/rd_add_proxy; ln_csi<-log(pmax(cspi,0)+CR_CSPI_SHIFT)
  eta<-rc("r0")+zsp+rc("r1")*ln_cr+rc("r2")*cr+rc("r3")*sqrt_ba+rc("r4")*ln_bal_ba+rc("r5")*cr_over_rd+rc("r6")*ln_csi
  r<-1/(1+exp(eta)); HCB1<-(1-cr)*HT; maxrec<-pmax(HT2-HCB1,0); HCB2<-HCB1+r*maxrec
  CR2<-1-HCB2/pmax(HT2,1e-3); pmin(pmax(CR2,0.01),0.95)
}

## build tree list (mirrors v2 mk_tl; keeps original SPCD; flags Greg coverage)
mk_tl<-function(trows,cov){
  d2<-trows[is.finite(DIAMETER)&DIAMETER>=1.0]; if(!nrow(d2)) return(NULL)
  SPCD<-as.integer(d2$SPECIES); dbh_in<-as.numeric(d2$DIAMETER)
  cr<-as.numeric(d2$CRRATIO); cr[!is.finite(cr)|cr<=0]<-NA; cr<-cr/100; cr[is.na(cr)]<-0.5; cr<-pmin(pmax(cr,0.05),0.95)
  ht<-as.numeric(d2$HT); ht[!is.finite(ht)|ht<=4.5]<-NA   # KEEP FEET for Greg eqs (est_dg/hg use ft)
  TPA<-as.numeric(d2$TREE_COUNT); TPA[!is.finite(TPA)|TPA<=0]<-1.0
  sw<-{v<-sw_tree_map[as.character(SPCD)];v[is.na(v)]<-0;as.integer(v)}
  ## growth species: SPCD in Greg DG set? else remap to stand modal in-set species (params borrowed)
  has_dg<-SPCD %in% dg_set; has_hg<-SPCD %in% hg_set
  GSPCD<-SPCD  # SPCD used for Greg DG/HG kernels
  inset<-SPCD[has_dg]
  fb<-if(length(inset)) as.integer(names(sort(table(inset),decreasing=TRUE))[1]) else dg_set[1]
  GSPCD[!has_dg]<-fb
  ## HG fallback: if remapped DG species also lacks HG, fall back to a HG-set modal
  has_hg_g<-GSPCD %in% hg_set
  insetH<-SPCD[SPCD %in% hg_set]; fbH<-if(length(insetH)) as.integer(names(sort(table(insetH),decreasing=TRUE))[1]) else hg_set[1]
  GSPCD_H<-GSPCD; GSPCD_H[!has_hg_g]<-fbH
  cr_zsp_v<-{v<-cr_zsp[as.character(SPCD)];v[is.na(v)]<-0;v}
  list(SPCD=SPCD,GSPCD=GSPCD,GSPCD_H=GSPCD_H,dbh_in=dbh_in,CR=cr,HT=ht,TPA=TPA,softwood=sw,
       has_dg=has_dg,has_hg=has_hg,has_mo=SPCD %in% mo_set,
       cr_zsp=cr_zsp_v,
       SDImax=cov$SDImax_brms,EMT=cov$EMT,TD=cov$TD,ELEV=cov$ELEVFT,cspi=NA_real_)
}
## competition (mirror v2): BA(ft2/ac), per-tree BAL(ft2/ac), CCFL(ft2/ac of larger trees), SDI, rd
recompute_comp<-function(tl){
  DBH<-tl$dbh_in; TPA<-tl$TPA; n<-length(DBH); ba_tree<-pi/4*DBH^2/144; BA_ac<-ba_tree*TPA; BA_tot<-sum(BA_ac)
  ord<-order(DBH,decreasing=TRUE); BAL<-numeric(n); cum<-0
  for(i in ord){ BAL[i]<-cum; cum<-cum+BA_ac[i] }
  ## CCFL = crown competition factor in larger trees; approximate with BAL (ft2/ac) as Greg's ccfl proxy
  ## (Greg's ccfl is "basal area in larger trees" style competition; we pass BAL).
  CCFL<-BAL
  SDI<-sum(TPA*(DBH/10)^1.605); rd_add<-SDI/pmax(tl$SDImax,1)
  tl$BA<-BA_tot; tl$BAL<-BAL; tl$CCFL<-CCFL; tl$SDI<-SDI; tl$rd_add<-rd_add; tl
}
stand_metrics<-function(tl){DBH<-tl$dbh_in; TPA<-tl$TPA; ok<-is.finite(DBH)&is.finite(TPA)&TPA>0; DBH<-DBH[ok]; TPA<-TPA[ok]; sT<-sum(TPA)
  list(BA=if(length(DBH)) sum(pi/4*DBH^2*TPA)/144 else NA_real_, QMD=if(is.finite(sT)&&sT>1e-6) sqrt(sum(DBH^2*TPA)/sT) else NA_real_, TPH=if(is.finite(sT)) sT*2.4710538 else NA_real_)}
add_recruits<-function(tl,n_add_tpa){
  if(n_add_tpa<=1e-6) return(tl)
  if(length(tl$dbh_in)){ ba_ac<-pi/4*tl$dbh_in^2/144*tl$TPA; dom<-tl$SPCD[which.max(ba_ac)] } else dom<-dg_set[1]
  domg<-if(dom %in% dg_set) dom else { inset<-tl$SPCD[tl$has_dg]; if(length(inset)) as.integer(names(sort(table(inset),decreasing=TRUE))[1]) else dg_set[1] }
  domh<-if(domg %in% hg_set) domg else hg_set[1]
  tl$SPCD<-c(tl$SPCD,dom); tl$GSPCD<-c(tl$GSPCD,domg); tl$GSPCD_H<-c(tl$GSPCD_H,domh)
  tl$dbh_in<-c(tl$dbh_in,ig_dbh_in); tl$CR<-c(tl$CR,0.6); tl$HT<-c(tl$HT,NA_real_); tl$TPA<-c(tl$TPA,n_add_tpa)
  swv<-sw_tree_map[as.character(dom)]; swv<-if(is.na(swv))0L else as.integer(swv)
  tl$softwood<-c(tl$softwood,swv)
  tl$has_dg<-c(tl$has_dg,dom%in%dg_set); tl$has_hg<-c(tl$has_hg,dom%in%hg_set); tl$has_mo<-c(tl$has_mo,dom%in%mo_set)
  tl$cr_zsp<-c(tl$cr_zsp,{v<-cr_zsp[as.character(dom)];if(is.na(v))0 else v})
  tl
}

## ========================== PARALLEL PROJECTION (mclapply) ==========================
## Identical equations/forms/params as the sequential driver. ONLY the stand loop is
## parallelized: each stand is fully independent (NO RNG inside the per-stand loop; the
## only set.seed is the stand-subsample above), so per-stand outputs are deterministic
## regardless of execution order. Tree cache for a state is built ONCE in the parent
## before forking, so workers share it copy-on-write (no re-reads). Results are
## collected and rbind'd in the SAME stand order as the sequential driver.
suppressPackageStartupMessages({ library(parallel) })
NCORES_CLI <- suppressWarnings(as.integer(ga("cores", NA)))
NCORES_ENV <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=NA)); if(is.na(NCORES_ENV)||NCORES_ENV<1){ ntsk<-suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS",unset=NA))); NCORES_ENV<-if(!is.na(ntsk)&&ntsk>1) ntsk else max(1L, detectCores()-1L) }; NCORES <- if(!is.na(NCORES_CLI)&&NCORES_CLI>=1L) NCORES_CLI else NCORES_ENV
cat(sprintf("STEP 2: projecting %d stands x %d cycles using mclapply mc.cores=%d ...\n", nrow(si), NCYC, NCORES))

project_one <- function(ix, spl){
  cn<-si$STAND_CN[ix]; trows<-spl[[cn]]; if(is.null(trows)||!nrow(trows)) return(NULL)
  cov<-si[ix]; st<-FIPS[as.character(cov$STATE)]; inv_year<-cov$INV_YEAR
  tl<-mk_tl(trows,cov); if(is.null(tl)||length(tl$dbh_in)<1) return(NULL)
  tl<-recompute_comp(tl)
  elev<-tl$ELEV; emt<-tl$EMT; td<-tl$TD
  tcw<-tl$TPA; tct<-sum(tcw)
  sfb_dg<-sum(tcw[!tl$has_dg])/tct; sfb_hg<-sum(tcw[!tl$has_hg])/tct; sfb_mo<-sum(tcw[!tl$has_mo])/tct
  fb_dg_n<-sum(tcw[!tl$has_dg]); fb_hg_n<-sum(tcw[!tl$has_hg]); fb_mo_n<-sum(tcw[!tl$has_mo])
  fbrow<-data.table(STAND_CN=cn,VARIANT=VARIANT,fb_dg=sfb_dg,fb_hg=sfb_hg,fb_mo=sfb_mo)
  miss<-!is.finite(tl$HT); if(any(miss)) tl$HT[miss]<-4.5*tl$dbh_in[miss]^0.5+4.5
  loc_rows<-vector("list",NCYC+1L); loc_tre<-vector("list",NCYC+1L); k<-0L
  emit<-function(cy){sm<-stand_metrics(tl); py<-cy*CYCLEN; k<<-k+1L
    cchm<-if(!is.null(tl$cch_last)&&length(tl$cch_last)==length(tl$TPA)&&sum(tl$TPA)>0) sum(tl$cch_last*tl$TPA)/sum(tl$TPA) else NA_real_
    loc_rows[[k]]<<-data.table(STAND_CN=cn,STATE=st,YEAR=inv_year+py,PROJ_YEAR=py,VARIANT=VARIANT,CONFIG=CONFIG,AGB_TONS_AC=NA_real_,BA_FT2AC=sm$BA,QMD_IN=sm$QMD,TPH=sm$TPH,CCH_MEAN=cchm)
    loc_tre[[k]]<<-data.table(STAND_CN=cn,CONFIG=CONFIG,PROJ_YEAR=py,SPCD=tl$SPCD,DBH_IN=tl$dbh_in,HT_M=tl$HT*0.3048,TPA=tl$TPA)}
  emit(0)
  for(cy in 1:NCYC){
    gg<-grp_organon(tl$SPCD)
    cch<-stand_cch(tl$dbh_in, tl$HT, tl$CR, tl$TPA, gg); tl$cch_last<-cch
    surv_cyc<-rep(1.0,length(tl$dbh_in))
    for(yr in 1:CYCLEN){
      dgan<-dg_annual(tl$GSPCD, tl$dbh_in, tl$CR, tl$HT, tl$CCFL, elev, emt)
      hgan<-hg_annual(tl$GSPCD_H, tl$HT, tl$CR, tl$CCFL, cch, elev, td, emt)
      psa <-gomp_surv_annual(tl$SPCD, tl$CR, cch)
      tl$dbh_in<-pmin(pmax(tl$dbh_in+dgan,0.1),200)
      tl$HT<-pmin(tl$HT+hgan,400)
      surv_cyc<-surv_cyc*psa
    }
    rd<-tl$rd_add; st_mult<-pmin(1+pmax(rd-1.0,0)*2.0, 2.0)
    H_cyc<- -log(pmin(pmax(surv_cyc,1e-9),1)); surv<-exp(-H_cyc*st_mult); surv[!is.finite(surv)]<-1
    tl$TPA<-tl$TPA*surv; tl$TPA[!is.finite(tl$TPA)]<-0
    BA_m<-tl$BA*0.2296; HT2m<-tl$HT*0.3048
    tl$CR<-cr_update(tl$SPCD,tl$CR,tl$HT*0.3048,HT2m,BA_m,tl$BAL,tl$BA,if(is.na(tl$cspi))1.0 else tl$cspi,tl$cr_zsp)
    keep<-tl$TPA>1e-4
    if(any(!keep)){ vn<-length(tl$dbh_in); for(nm in names(tl)) if(length(tl[[nm]])==vn) tl[[nm]]<-tl[[nm]][keep] }
    tl<-add_recruits(tl, IG_TPA*CYCLEN)
    tl$HT[!is.finite(tl$HT)]<-4.5*tl$dbh_in[!is.finite(tl$HT)]^0.5+4.5
    tl<-recompute_comp(tl); emit(cy)
  }
  list(rows=rbindlist(loc_rows[seq_len(k)]),
       tre =rbindlist(loc_tre[seq_len(k)]),
       fb  =fbrow,
       fb_dg_n=fb_dg_n, fb_hg_n=fb_hg_n, fb_mo_n=fb_mo_n, tc=tct)
}

state_groups<-split(seq_len(nrow(si)),si$STATE)
all_res<-vector("list", length(state_groups)); gi<-0L
for(stcode in names(state_groups)){
  spl<-get_state_trees(as.integer(stcode))
  gi<-gi+1L
  if(is.null(spl)) next
  if(is.list(spl)&&length(spl)==1&&is.null(spl[[1]])) next
  idxs<-state_groups[[stcode]]
  res<-mclapply(idxs, project_one, spl=spl, mc.cores=NCORES, mc.preschedule=TRUE)
  errs<-vapply(res, function(x) inherits(x,"try-error"), logical(1))
  if(any(errs)){ cat("  WORKER ERRORS in state",stcode,":",sum(errs),"-- first:\n"); print(res[[which(errs)[1]]]); stop("mclapply worker error") }
  all_res[[gi]]<-res
}
flat<-unlist(all_res, recursive=FALSE)
flat<-flat[!vapply(flat, is.null, logical(1))]
nproj<-length(flat)
out<-rbindlist(lapply(flat, `[[`, "rows"))
tre<-rbindlist(lapply(flat, `[[`, "tre"))
psf<-rbindlist(lapply(flat, `[[`, "fb"))
fb_dg_acc<-sum(vapply(flat, `[[`, numeric(1), "fb_dg_n"))
fb_hg_acc<-sum(vapply(flat, `[[`, numeric(1), "fb_hg_n"))
fb_mo_acc<-sum(vapply(flat, `[[`, numeric(1), "fb_mo_n"))
tc_acc<-sum(vapply(flat, `[[`, numeric(1), "tc"))

otag<-sprintf("conus_eq_%s_%s",tolower(VARIANT),CONFIG)
fwrite(out,file.path(OUTD,paste0(otag,"_metrics.csv"))); fwrite(tre,file.path(OUTD,paste0(otag,"_treelists.csv")))
fwrite(psf,file.path(OUTD,paste0(otag,"_fallback.csv")))
cat("\n  stands projected:",nproj,"\n")
cat(sprintf("  Greg coverage (tree-count weighted): DG=%.1f%% (fallback %.1f%%)  HG=%.1f%% (fallback %.1f%%)  MORT=%.1f%% (fallback %.1f%%)\n",
  100*(1-fb_dg_acc/tc_acc),100*fb_dg_acc/tc_acc,100*(1-fb_hg_acc/tc_acc),100*fb_hg_acc/tc_acc,100*(1-fb_mo_acc/tc_acc),100*fb_mo_acc/tc_acc))
cat("Wrote:",file.path(OUTD,paste0(otag,"_metrics.csv")),"(",nrow(out),"rows )\n")
cat("Wrote:",file.path(OUTD,paste0(otag,"_treelists.csv")),"(",nrow(tre),"rows )\n")
cat("Wrote:",file.path(OUTD,paste0(otag,"_fallback.csv")),"(",nrow(psf),"stands )\n")
cat("\n=== year 0/50/100 stand-mean metrics (",CONFIG,") ===\n")
print(out[PROJ_YEAR %in% c(0,50,100),.(BA=mean(BA_FT2AC,na.rm=TRUE),QMD=mean(QMD_IN,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE),n=.N),by=PROJ_YEAR])
cat("DONE.\n")
