#!/usr/bin/env Rscript
## conus_eq_projector.R -- 100yr (20x5yr) stand projection using fvs-conus fitted equations.
## species_mode toggle: --mode=dependent -> conus_b2 (v8 species-aware DG);
##                       --mode=free      -> conus_b1 (speciesfree DG, W*gamma only).
## Wired: DG(b2/b1), mortality(logit), ht-dbh(Wykoff static), CR-recession(Hann-Hanus),
##        dynamic competition recompute each cycle. Stubbed: ingrowth, height-increment HG.
## Reuses v8_re_means.rds + fourarm eta forms (cspiv6 pairs; DBH inches).
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){ m<-grep(paste0("^--",n,"="),args,value=TRUE); if(!length(m)) return(d); sub(paste0("^--",n,"="),"",m[1]) }
ROOT<-"/users/PUOM0008/crsfaaron/fvs-conus"; SCR<-"/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
OUTD<-ga("outdir",file.path(SCR,"conus_eq_proj")); MODE<-ga("mode","dependent"); VARIANT<-toupper(ga("variant","NE"))
NSTAND<-as.integer(ga("nstands","300")); SEED<-as.integer(ga("seed","7")); NCYC<-as.integer(ga("ncycles","20")); CYCLEN<-as.integer(ga("cyclelen","5"))
PAIRS<-ga("pairs",file.path(ROOT,"data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds"))
TRAITS<-ga("traits",file.path(ROOT,"traits/species_traits.rds"))
V8DIR<-file.path(ROOT,"output/conus/dg_kue/v8"); V8SUM<-file.path(V8DIR,"dg_kuehne_v8_100k_prod_summary.csv"); V8META<-file.path(V8DIR,"dg_kuehne_v8_100k_prod_meta.rds")
RE_RDS<-file.path(SCR,"smoke_conus_eq/v8_re_means.rds")
B1FIT<-file.path(ROOT,"output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds"); B1META<-file.path(ROOT,"output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_meta.rds")
CC<-file.path(ROOT,"output/conus")
MORT_SUM<-file.path(CC,"mort_logit_simple_cspi_traits1_fixed_summary.csv"); MORT_SP<-file.path(CC,"mort_logit_simple_cspi_traits1_species_intercepts.csv")
HTD_SUM<-file.path(CC,"htdbh_wykoff_lognormal_cspi_traits1_fixed_summary.csv"); HTD_SP<-file.path(CC,"htdbh_wykoff_lognormal_cspi_traits1_species_intercepts.csv")
CR_SUM<-file.path(CC,"cr_recession_cspi_traits1_fixed_summary.csv"); CR_SP<-file.path(CC,"cr_recession_cspi_traits1_species_intercepts.csv")
stopifnot(MODE %in% c("dependent","free")); CONFIG<-if(MODE=="dependent")"conus_b2" else "conus_b1"; dir.create(OUTD,recursive=TRUE,showWarnings=FALSE)
cat("==== conus_eq_projector ====\n"); cat(sprintf("  variant=%s mode=%s (CONFIG=%s) nstands=%d %dx%dyr\n",VARIANT,MODE,CONFIG,NSTAND,NCYC,CYCLEN))
cat("  WIRED: DG(b2/b1 toggle), mortality(logit), ht-dbh(Wykoff static), CR-recession, dynamic competition.\n")
cat("  STUBBED: ingrowth/recruitment; height-increment HG (static ht-dbh supplies HT).\n\n")
for(f in c(PAIRS,TRAITS,V8SUM,V8META,RE_RDS,MORT_SUM,MORT_SP,HTD_SUM,HTD_SP,CR_SUM,CR_SP)) if(!file.exists(f)) stop("MISSING: ",f)
if(MODE=="free") for(f in c(B1FIT,B1META)) if(!file.exists(f)) stop("MISSING b1: ",f)
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
cat("STEP 1: load pairs + seed treelists ...\n"); d<-as.data.table(readRDS(PAIRS)); d<-d[fvs_variant==VARIANT]; d[,dbh_in:=DBH1]
filt<-is.finite(d$DBH1)&d$DBH1>=1.0&is.finite(d$CR1)&d$CR1>0&d$CR1<=1.0&is.finite(d$HT1)&d$HT1>=1.37&is.finite(d$BA1)&d$BA1>0&is.finite(d$BAL1)&d$BAL1>=0&is.finite(d$BAL_SW1)&is.finite(d$BAL_HW1)&is.finite(d$bgi)&!is.na(d$EPA_L1_CODE)&d$EPA_L1_CODE!=""&!is.na(d$FORTYPCD_cond1)&d$FORTYPCD_cond1>0&is.finite(d$SDImax_brms)&d$SDImax_brms>0&is.finite(d$sdi_additive1)&d$sdi_additive1>0&is.finite(d$TPH_UNADJ1)&d$TPH_UNADJ1>0&d$SPCD %in% sp_levels
d<-d[filt]; cat("  usable trees:",nrow(d),"\n"); d[,TPA:=TPH_UNADJ1/2.4710538]   # per-acre expansion (verified: yields physical stand BA ~272 ft2/ac median; TPA1 is non-normalized)
sw_tree_map<-setNames(traits$softwood,as.character(traits$SPCD)); d[,softwood:={v<-sw_tree_map[as.character(SPCD)]; v[is.na(v)]<-0; as.integer(v)}]
cnt<-d[,.N,by=plot_key][N>=8&N<=80]; set.seed(SEED); keys<-cnt[sample(.N,min(NSTAND,.N))]$plot_key; d<-d[plot_key %in% keys]
cat("  stands selected:",length(keys)," trees:",nrow(d),"\n")
if(!"STATECD" %in% names(d)) d[,STATECD:=NA_integer_]
FIPS<-c("9"="CT","23"="ME","25"="MA","33"="NH","34"="NJ","36"="NY","42"="PA","44"="RI","50"="VT","24"="MD","10"="DE","54"="WV")
if(MODE=="free") b1_pack<-load_b1()
sw_mean<-mean(sw_by_sp)
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
mk_tl<-function(p){tl<-list(SPCD=p$SPCD,sp_idx=match(p$SPCD,sp_levels),L1_idx=match(as.character(p$EPA_L1_CODE),L1_lev),L2_idx=match(as.character(p$EPA_L2_CODE),L2_lev),L3_idx=match(as.character(p$EPA_L3_CODE),L3_lev),FT_idx=match(as.integer(p$FORTYPCD_cond1),FT_lev),dbh_in=p$dbh_in,CR=p$CR1,HT=p$HT1,TPA=p$TPA,softwood=p$softwood,bgi=p$bgi,cspi=p$cspi,SDImax=p$SDImax_brms,sw_mean=sw_mean,mort_zsp={v<-mort_zsp[as.character(p$SPCD)];v[is.na(v)]<-0;v},htd_zsp={v<-htd_zsp[as.character(p$SPCD)];v[is.na(v)]<-0;v},cr_zsp={v<-cr_zsp[as.character(p$SPCD)];v[is.na(v)]<-0;v},mort_cspi_shift=mort_cspi_shift,htd_cspi_shift=htd_cspi_shift,cr_cspi_shift=cr_cspi_shift)
  tl$Wrow<-W[tl$sp_idx,,drop=FALSE]
  if(MODE=="free"){tl$b1_L1<-match(as.character(p$EPA_L1_CODE),b1_pack$pm$L1);tl$b1_L2<-match(as.character(p$EPA_L2_CODE),b1_pack$pm$L2);tl$b1_L3<-match(as.character(p$EPA_L3_CODE),b1_pack$pm$L3);tl$b1_L1[is.na(tl$b1_L1)]<-1;tl$b1_L2[is.na(tl$b1_L2)]<-1;tl$b1_L3[is.na(tl$b1_L3)]<-1}
  ok<-!is.na(tl$sp_idx)&!is.na(tl$L1_idx)&!is.na(tl$L2_idx)&!is.na(tl$L3_idx)&!is.na(tl$FT_idx)
  for(nm in names(tl)) if(length(tl[[nm]])==length(ok)) tl[[nm]]<-tl[[nm]][ok]
  tl$Wrow<-tl$Wrow[ok,,drop=FALSE]; tl}
recompute_comp<-function(tl){DBH<-tl$dbh_in; TPA<-tl$TPA; n<-length(DBH); ba_tree<-pi/4*DBH^2/144; BA_ac<-ba_tree*TPA; BA_tot<-sum(BA_ac); ord<-order(DBH,decreasing=TRUE)
  BAL<-numeric(n); BAL_SW<-numeric(n); BAL_HW<-numeric(n); cum<-0; cum_sw<-0; cum_hw<-0
  for(i in ord){BAL[i]<-cum; BAL_SW[i]<-cum_sw; BAL_HW[i]<-cum_hw; cum<-cum+BA_ac[i]; if(tl$softwood[i]==1) cum_sw<-cum_sw+BA_ac[i] else cum_hw<-cum_hw+BA_ac[i]}
  SDI<-sum(TPA*(DBH/10)^1.605); rd_add<-SDI/pmax(tl$SDImax,1); sdi_cx<-SDI/pmax(SDI,1.0)
  tl$BA<-BA_tot; tl$BAL<-BAL; tl$BAL_SW<-BAL_SW; tl$BAL_HW<-BAL_HW; tl$SDI<-SDI; tl$rd_add<-rd_add; tl$sdi_cx<-sdi_cx; tl}
stand_metrics<-function(tl){DBH<-tl$dbh_in; TPA<-tl$TPA; ok<-is.finite(DBH)&is.finite(TPA)&TPA>0; DBH<-DBH[ok]; TPA<-TPA[ok]; sT<-sum(TPA); list(BA=if(length(DBH)) sum(pi/4*DBH^2*TPA)/144 else NA_real_, QMD=if(is.finite(sT)&&sT>1e-6) sqrt(sum(DBH^2*TPA)/sT) else NA_real_, TPH=if(is.finite(sT)) sT*2.4710538 else NA_real_)}
cat("STEP 2: projecting",length(keys),"stands x",NCYC,"cycles ...\n"); rows<-list(); treerows<-list(); ri<-0; ti<-0
for(kk in keys){p<-d[plot_key==kk]; tl<-mk_tl(p); if(length(tl$dbh_in)<3) next
  cn<-as.character(p$PLT_CN_cond1[1]); if(is.na(cn)||cn=="") cn<-as.character(kk); st<-if(!is.na(p$STATECD[1])) FIPS[as.character(p$STATECD[1])] else NA; inv_year<-2010L
  tl<-recompute_comp(tl)
  emit<-function(cy){sm<-stand_metrics(tl); py<-cy*CYCLEN; ri<<-ri+1
    rows[[ri]]<<-data.table(STAND_CN=cn,STATE=st,YEAR=inv_year+py,PROJ_YEAR=py,VARIANT=VARIANT,CONFIG=CONFIG,AGB_TONS_AC=NA_real_,BA_FT2AC=sm$BA,QMD_IN=sm$QMD,TPH=sm$TPH)
    HT<-ht_from_dbh(tl); ti<<-ti+1; treerows[[ti]]<<-data.table(STAND_CN=cn,CONFIG=CONFIG,PROJ_YEAR=py,SPCD=tl$SPCD,DBH_IN=tl$dbh_in,HT_M=HT,TPA=tl$TPA)}
  emit(0)
  for(cy in 1:NCYC){eta<-if(MODE=="dependent") eta_dg_b2(tl) else eta_dg_b1(tl,b1_pack); dg_a<-exp(eta+sigma8^2/2); dg_a[!is.finite(dg_a)]<-0; dg_a<-pmin(dg_a,2.0); tl$dbh_in<-pmin(pmax(tl$dbh_in+dg_a*CYCLEN,0.1),200)
    HT2<-ht_from_dbh(tl); tl$CR<-cr_update(tl,HT2); tl$HT<-HT2
    p_die_a<-1/(1+exp(-eta_mort(tl))); p_die_a[!is.finite(p_die_a)]<-0; surv<-(1-pmin(pmax(p_die_a,0),0.999))^CYCLEN; tl$TPA<-tl$TPA*surv; tl$TPA[!is.finite(tl$TPA)]<-0
    tl<-recompute_comp(tl); emit(cy)}}
out<-rbindlist(rows); tre<-rbindlist(treerows); otag<-sprintf("conus_eq_%s_%s",tolower(VARIANT),CONFIG)
fwrite(out,file.path(OUTD,paste0(otag,"_metrics.csv"))); fwrite(tre,file.path(OUTD,paste0(otag,"_treelists.csv")))
cat("\nWrote:",file.path(OUTD,paste0(otag,"_metrics.csv")),"(",nrow(out),"rows )\n"); cat("Wrote:",file.path(OUTD,paste0(otag,"_treelists.csv")),"(",nrow(tre),"rows )\n")
cat("\n=== year 0/50/100 stand-mean metrics (",CONFIG,") ===\n")
print(out[PROJ_YEAR %in% c(0,50,100),.(BA=mean(BA_FT2AC,na.rm=TRUE),QMD=mean(QMD_IN,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE),n=.N),by=PROJ_YEAR])
cat("DONE.\n")
