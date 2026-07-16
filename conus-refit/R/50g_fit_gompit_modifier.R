#!/usr/bin/env Rscript
# =============================================================================
# 50g_fit_gompit_modifier.R
# Extract gompit-survival eta_base from the refit surv_crz_100k outputs, build
# the disturbance-modifier bundle from pairs_v2, and fit
# modifier_gompit_speciesdep.stan (species-dependent modifier on the deployed
# gompit mortality form).
#
# eta_base = b0 + trait_effect[sp] + z_L1+z_L2+z_L3+z_FT
#          + b1*dbh + b2*dbh_sq + b3*cr_z + b3b*cr_z_sq + b4*ln_csi
#          + b5*bal_metric + b6*sqrt_ba_rd + b7*cch_z + b7b*cch_z_sq
#   trait_effect = W . gamma  (species-FREE trait part; z_sp not saved by the
#   refit so eta_base is the species-free survival base. CAVEAT: the modifier's
#   per-species RE absorbs species-specific base error along with the true
#   disturbance-by-species interaction; GLOBAL modifier alphas are clean. For the
#   full species-specific base, re-fit survival saving z_sp.)
# Standardization (cch_z, cr_z) reproduces the stan transformed-data block (mean
# and sd over the 100k survival-fit subsample, same 34c prep + seed 2026).
# Usage: Rscript 50g_fit_gompit_modifier.R --smoke
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(cmdstanr); library(posterior) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){ m<-grep(paste0("^--",n,"="),args,value=TRUE); if(!length(m)) return(d); sub(paste0("^--",n,"="),"",m[1]) }
SMOKE <- any(grepl("^--smoke$", args))
SURV  <- ga("surv","output/conus/mort/surv_unified_v2_crz/surv_crz_100k")
PAIRS <- ga("pairs","calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS<- ga("traits","calibration/traits/species_traits_v2.rds")
STAN  <- ga("stan","calibration/stan/modifier_gompit_speciesdep.stan")
OUT   <- ga("outdir","output/conus/mort/gompit_modifier"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
SUB   <- as.integer(ga("subsample", if(SMOKE) "60000" else "300000"))

## 1. base-fit coefficients (posterior means)
pd <- readRDS(paste0(SURV,"_param_draws.rds"))
pm <- function(v) mean(as.numeric(pd[,v]))
b0<-pm("b0");b1<-pm("b1");b2<-pm("b2");b3<-pm("b3");b3b<-pm("b3b");b4<-pm("b4")
b5<-pm("b5");b6<-pm("b6");b7<-pm("b7");b7b<-pm("b7b")
gamma <- sapply(sprintf("gamma[%d]",1:8), pm)
reL1<-fread(paste0(SURV,"_re_L1.csv")); reL2<-fread(paste0(SURV,"_re_L2.csv"))
reL3<-fread(paste0(SURV,"_re_L3.csv")); reFT<-fread(paste0(SURV,"_re_FT.csv"))
zL1<-setNames(reL1$mean,as.character(reL1$level)); zL2<-setNames(reL2$mean,as.character(reL2$level))
zL3<-setNames(reL3$mean,as.character(reL3$level)); zFT<-setNames(reFT$mean,as.character(reFT$level))
rsp_path<-paste0(SURV,"_re_sp.csv")
if(file.exists(rsp_path)){ resp<-fread(rsp_path); zsp<-setNames(resp$mean,as.character(resp$level)); cat(sprintf("z_sp loaded: %d species (species-dependent base)\n",nrow(resp))) } else { zsp<-setNames(numeric(0),character(0)); cat("z_sp absent: species-FREE base (caveat)\n") }

## 2. replicate 34c survival prep (filter + variety splits + levels)
dat <- as.data.table(readRDS(PAIRS)); traits <- as.data.table(readRDS(TRAITS))
cc<-intersect(c("CCH","cch","CCH1","CCH_TT"),names(dat))[1]
if(is.na(cc)) dat[,cch:=0] else { dat[,cch:=as.numeric(get(cc))]; dat[!is.finite(cch),cch:=0] }
dat[,alive:=as.integer(TREESTATUS2==1)]; dat[,T_years:=YEARS]
if("climate_si"%in%names(dat)){ med<-median(dat$climate_si,na.rm=TRUE); dat[!is.finite(climate_si),climate_si:=med]; dat[,ln_csi:=log(pmax(climate_si,0.1))] } else dat[,ln_csi:=0]
dat[!is.finite(ln_csi),ln_csi:=0]
dat[,rd_ratio:=sdi_additive1/SDImax_brms]
dat[,sqrt_ba_rd:=sqrt(pmax(BA1*0.2296,0)*pmax(rd_ratio,0))]
dat[,bal_metric:=BAL_SW1+BAL_HW1]
dat <- dat[TREESTATUS1==1 & !is.na(TREESTATUS2)&TREESTATUS2%in%c(1,2) &
  is.finite(DBH1)&DBH1>=2.54 & is.finite(CR1)&CR1>0&CR1<=1 & is.finite(cch) &
  is.finite(YEARS)&YEARS>=1&YEARS<=20 & !is.na(EPA_L1_CODE)&!is.na(EPA_L2_CODE)&!is.na(EPA_L3_CODE) &
  EPA_L1_CODE!=""&EPA_L2_CODE!=""&EPA_L3_CODE!="" & is.finite(BA1)&BA1>=0 &
  is.finite(BAL_SW1)&BAL_SW1>=0 & is.finite(BAL_HW1)&BAL_HW1>=0 &
  is.finite(rd_ratio)&rd_ratio>=0 & !is.na(FORTYPCD_cond1)&FORTYPCD_cond1>0]
if(any(traits$SPCD==2020L)){ dat[SPCD==202L & as.character(EPA_L1_CODE)=="7",SPCD:=2020L]; dat[SPCD==202L,SPCD:=2021L] }
if(any(traits$SPCD==1080L)){ dat[SPCD==108L & as.character(EPA_L1_CODE)=="7",SPCD:=1080L]; dat[SPCD==108L,SPCD:=1081L] }
spk<-dat[,.N,by=SPCD][N>=5000]; dat<-dat[SPCD%in%spk$SPCD]
sp_levels<-sort(unique(dat$SPCD))
trait_cols<-c("wood_specific_gravity","shade_tolerance_num","softwood","leaf_longevity_months","max_ht_m","max_dbh_cm","vulnerability_score","sensitivity")
tr<-traits[match(sp_levels,SPCD),c("SPCD",trait_cols),with=FALSE]; W<-as.matrix(tr[,trait_cols,with=FALSE])
for(j in seq_len(ncol(W))){ na<-is.na(W[,j]); if(any(na))W[na,j]<-median(W[!na,j],na.rm=TRUE); W[,j]<-(W[,j]-mean(W[,j]))/sd(W[,j]) }
trait_effect <- as.numeric(W %*% gamma)

set.seed(2026); sub_idx<-if(nrow(dat)>100000) sort(sample.int(nrow(dat),100000)) else seq_len(nrow(dat))
cch_mean<-mean(dat$cch[sub_idx]); cch_sd<-sd(dat$cch[sub_idx]); if(cch_sd<=0)cch_sd<-1
cr_mean<-mean(dat$CR1[sub_idx]); cr_sd<-sd(dat$CR1[sub_idx]); if(cr_sd<=0)cr_sd<-1

## 3. eta_base (species-free base)
dat[,sp_idx:=match(SPCD,sp_levels)]
gv<-function(map,key){ v<-map[as.character(key)]; v[is.na(v)]<-0; as.numeric(v) }
dat[,cch_z:=(cch-cch_mean)/cch_sd]; dat[,cr_z:=(CR1-cr_mean)/cr_sd]
dat[,eta_base := b0 + trait_effect[sp_idx] + gv(zsp,SPCD) +
   gv(zL1,EPA_L1_CODE)+gv(zL2,EPA_L2_CODE)+gv(zL3,EPA_L3_CODE)+gv(zFT,FORTYPCD_cond1) +
   b1*DBH1 + b2*DBH1^2 + b3*cr_z + b3b*cr_z^2 + b4*ln_csi +
   b5*bal_metric + b6*sqrt_ba_rd + b7*cch_z + b7b*cch_z^2]
dat[,p_surv:=exp(-exp(-pmin(pmax(eta_base,-5),20))*T_years)]
cat(sprintf("VALIDATION: mean predicted survival=%.4f  mean observed alive=%.4f  (diff %.4f)\n",
            mean(dat$p_surv,na.rm=TRUE), mean(dat$alive), mean(dat$p_surv,na.rm=TRUE)-mean(dat$alive)))

## 4. modifier bundle
L1_levels<-as.character(reL1$level); dat[,L1_idx:=match(as.character(EPA_L1_CODE),L1_levels)]
dat<-dat[!is.na(L1_idx) & !is.na(sp_idx)]
z0<-function(x){ x<-as.numeric(x); x[!is.finite(x)]<-0; x }
dd <- dat[, .(alive, eta_base, T_years, sp_idx, L1_idx,
   is_plantation=z0(is_plantation), d_fire=z0(had_fire_t1), d_insect=z0(had_insect_t1),
   d_disease=z0(had_disease_t1), d_wind=z0(had_wind_t1), d_harvest=z0(had_harvest_t1),
   dstrb_decay=z0(dstrb_decay_5yr), t_cutting=z0(had_cutting_t1),
   t_site_prep=z0(had_site_prep_t1), trt_decay=z0(trt_decay_5yr))]
if(SUB<nrow(dd)){ set.seed(7); dd<-dd[sort(sample.int(nrow(dd),SUB))] }
cat(sprintf("modifier rows=%s  disturbed frac: fire=%.3f insect=%.3f harvest=%.3f plant=%.3f\n",
  format(nrow(dd),big.mark=","), mean(dd$d_fire>0), mean(dd$d_insect>0), mean(dd$d_harvest>0), mean(dd$is_plantation>0)))

stan_data<-list(N_obs=nrow(dd), N_L1=length(L1_levels), N_sp=length(sp_levels), P_trait=ncol(W),
  alive=as.integer(dd$alive), eta_base=dd$eta_base, T_years=dd$T_years,
  is_plantation=dd$is_plantation, d_fire=dd$d_fire, d_insect=dd$d_insect, d_disease=dd$d_disease,
  d_wind=dd$d_wind, d_harvest=dd$d_harvest, dstrb_decay=dd$dstrb_decay,
  t_cutting=dd$t_cutting, t_site_prep=dd$t_site_prep, trt_decay=dd$trt_decay,
  L1_idx=dd$L1_idx, sp_idx=dd$sp_idx, W=W)

## 5. fit
mod<-cmdstan_model(STAN)
if(SMOKE){iw<-150;is_<-150;ch<-2}else{iw<-1000;is_<-1000;ch<-4}
fit<-mod$sample(data=stan_data,chains=ch,parallel_chains=ch,iter_warmup=iw,iter_sampling=is_,
  seed=42,adapt_delta=0.9,max_treedepth=10,refresh=50)
if(!SMOKE && nrow(dd)<=120000){ ll<-fit$draws("log_lik",format="draws_matrix"); lr<-loo::loo(ll); saveRDS(lr,file.path(OUT,"modifier_loo.rds")); cat(sprintf("LOO modifier elpd=%.1f (SE %.1f)\n",lr$estimates["elpd_loo","Estimate"],lr$estimates["elpd_loo","SE"])); rm(ll); gc() } else if(!SMOKE){ cat("skip loo save (N>120k, log_lik matrix too large)\n") }
vars<-c("alpha_0","alpha_fire","alpha_insect","alpha_harvest","alpha_cutting","alpha_plant",
  "sigma_sp_fire","sigma_sp_insect","sigma_sp_harvest","sigma_L1")
summ<-fit$summary(variables=vars,"mean","sd",~quantile(.x,c(.05,.95)),"rhat","ess_bulk")
print(summ); fwrite(summ, file.path(OUT,"gompit_modifier_speciesdep_summary.csv"))
saveRDS(list(summary=summ, n=nrow(dd), sp_levels=sp_levels), file.path(OUT,"gompit_modifier_speciesdep_meta.rds"))
cat("\nDONE. alpha_fire/insect/harvest = global disturbance survival effects (gompit scale); sigma_sp_* = species spread.\n")
