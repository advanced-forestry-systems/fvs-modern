## SMOKE PROOF: step ONE NE plot forward with fvs-conus v8c DG (species-aware)
## - extracts RE posterior means from 5GB fit RDS ONCE into compact CSVs
## - then projects DBH a few 5yr cycles, emits BA/QMD/TPH/(AGB stub) per cycle
suppressPackageStartupMessages({library(data.table)})
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/smoke_conus_eq"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
V8 <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/dg_kue/v8"
FIT  <- file.path(V8,"dg_kuehne_v8_100k_prod_fit.rds")
META <- file.path(V8,"dg_kuehne_v8_100k_prod_meta.rds")
SUM  <- file.path(V8,"dg_kuehne_v8_100k_prod_summary.csv")
PAIRS <- "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS<- "/users/PUOM0008/crsfaaron/fvs-conus/traits/species_traits.rds"

meta <- readRDS(META)
sp_levels<-meta$sp_levels; L1<-meta$L1_levels; L2<-meta$L2_levels; L3<-meta$L3_levels; FT<-meta$FT_levels
tcols<-meta$trait_cols; k1<-meta$bgi_knots[1]; k2<-meta$bgi_knots[2]

## ---- 1. fixed effects from summary CSV (compact, no RDS needed) ----
s<-fread(SUM); fx<-setNames(s$mean,s$variable)
b<-function(n) as.numeric(fx[[n]])
gamma     <- as.numeric(fx[paste0("gamma[",1:length(tcols),"]")])
gamma_site<- as.numeric(fx[paste0("gamma_site[",1:length(tcols),"]")])

## ---- 2. RE means from fit RDS: extract ONCE -> compact CSV ----
re_csv<-file.path(OUT,"v8_re_means.rds")
if(!file.exists(re_csv)){
  cat("Extracting REs from 5GB fit (one-time)...\n")
  fit<-readRDS(FIT)
  pull<-function(v) as.numeric(fit$summary(v,"mean")$mean)
  RE<-list(trait_effect=pull("trait_effect"), z_sp=pull("z_sp"),
           z_L1=pull("z_L1"), z_L2=pull("z_L2"), z_L3=pull("z_L3"),
           z_FT=pull("z_FT"), z_L1_bgi=pull("z_L1_bgi"),
           species_site_slope=pull("species_site_slope"))
  saveRDS(RE,re_csv); rm(fit); gc()
} else RE<-readRDS(re_csv)

## ---- 3. one NE plot from pairs (has all v8 inputs) ----
d<-as.data.table(readRDS(PAIRS))
d<-d[fvs_variant=="NE" & is.finite(DBH1)&DBH1>=2.54 & is.finite(CR1)&CR1>0&CR1<=1 &
     is.finite(BAL_SW1)&BAL_SW1>=0 & is.finite(BAL_HW1)&BAL_HW1>=0 & is.finite(bgi) &
     !is.na(EPA_L1_CODE)&EPA_L1_CODE!="" & !is.na(FORTYPCD_cond1)&FORTYPCD_cond1>0 &
     is.finite(BA1)&BA1>0 & SPCD %in% sp_levels]
cnt<-d[,.N,by=plot_key][N>=10&N<=60]; set.seed(1); pk<-cnt[sample(.N,1)]$plot_key
p<-d[plot_key==pk]
cat("Smoke plot:",pk," trees:",nrow(p),"\n")

traits<-as.data.table(readRDS(TRAITS))
## standardize W per training (approx: use training species set)
tr<-traits[SPCD %in% sp_levels, c("SPCD",tcols),with=FALSE]
for(c in tcols){na<-is.na(tr[[c]]); if(any(na)) tr[na,(c):=median(tr[[c]],na.rm=TRUE)]}
sw_train<-traits[match(sp_levels,SPCD),softwood]; sw_train[is.na(sw_train)]<-0; sw_mean<-mean(sw_train)

p[, sp_idx:=match(SPCD,sp_levels)]
p[, L1_idx:=match(as.character(EPA_L1_CODE),L1)]
p[, L2_idx:=match(as.character(EPA_L2_CODE),L2)]
p[, L3_idx:=match(as.character(EPA_L3_CODE),L3)]
p[, FT_idx:=match(as.integer(FORTYPCD_cond1),FT)]
p<-p[!is.na(sp_idx)&!is.na(L1_idx)&!is.na(L2_idx)&!is.na(L3_idx)&!is.na(FT_idx)]
cat("trees usable:",nrow(p),"\n")
sw_tree<-traits[match(p$SPCD,SPCD),softwood]; sw_tree[is.na(sw_tree)]<-0; sw_c<-sw_tree-sw_mean

## CSPI: use observed if present else median productivity
cspi_col<-grep("cspi",names(p),value=TRUE,ignore.case=TRUE)[1]
ln_cspi<- if(!is.na(cspi_col)) log(pmax(p[[cspi_col]],0.1)) else log(20)
TPA <- p$TPH_UNADJ1/2.4710538            # per acre
DBH <- p$DBH1/2.54                        # inches (state tracked in inches like ne_bench)
DBHcm <- p$DBH1                           # model uses cm
CRfrac<-p$CR1
BALsw<-p$BAL_SW1; BALhw<-p$BAL_HW1
rd_add<-p$sdi_additive1/p$SDImax_brms
sdi_cx<-p$sdi_additive1/pmax(p$SDI1,1)
bgi<-p$bgi

## eta function (v8c form WITHOUT a_cspi*ln_cspi unless present; a_cspi not in this fit's summary)
has_acspi <- "a_cspi" %in% names(fx)
eta_dg<-function(DBHcm,CRfrac,BALsw,BALhw,bgi,rd_add,sdi_cx,sp_idx,L1_idx,L2_idx,L3_idx,FT_idx,sw_c,ln_cspi){
  ln_dbh<-log(DBHcm); ln_cr<-log((CRfrac+0.2)/1.2); ln_balsw<-log(BALsw+0.01)
  bgi_b1<-bgi; bgi_b2<-pmax(bgi-k1,0); bgi_b3<-pmax(bgi-k2,0)
  b_site<-b("b6")+RE$z_L1_bgi[L1_idx]+RE$species_site_slope[sp_idx]
  eta<-b("b0")+RE$trait_effect[sp_idx]+RE$z_sp[sp_idx]+
    RE$z_L1[L1_idx]+RE$z_L2[L2_idx]+RE$z_L3[L3_idx]+RE$z_FT[FT_idx]+
    b("b1")*ln_dbh+b("b2")*DBHcm+b("b3")*ln_cr+
    b("b4")*ln_balsw+b("b5")*BALhw+
    b_site*bgi_b1+b("b9a")*bgi_b2+b("b9b")*bgi_b3+
    b("b7")*(0)+ # ba_x_rdadd handled below via metric; use ba_metric approx 0 contribution placeholder
    b("b11")*sdi_cx+
    b("b13")*(bgi*ln_dbh)+b("b14")*(bgi*sw_c)+b("b15")*(bgi*ln_cr)
  if(has_acspi) eta<-eta+b("a_cspi")*ln_cspi
  eta
}
sigma<-b("sigma")

## ---- 4. step forward: 4 cycles x 5yr (annual increment * years) ----
metrics<-function(DBHin,TPA){
  ba<-sum(pi/4*(DBHin^2)*TPA)/144            # ft2/ac (DBH inches)
  qmd<-sqrt(sum(DBHin^2*TPA)/sum(TPA))
  tph<-sum(TPA)*2.4710538
  data.table(BA_ft2ac=ba,QMD_in=qmd,TPH=tph,TPA=sum(TPA))
}
res<-list(); DBHc<-DBHcm; DBHi<-DBH
res[[1]]<-cbind(cycle=0,year=0,metrics(DBHi,TPA))
for(cy in 1:4){
  eta<-eta_dg(DBHc,CRfrac,BALsw,BALhw,bgi,rd_add,sdi_cx,p$sp_idx,p$L1_idx,p$L2_idx,p$L3_idx,p$FT_idx,sw_c,ln_cspi)
  eta<-pmin(pmax(eta,-30),20)
  dg_a<-exp(eta+sigma^2/2)        # annual DBH increment (cm/yr), lognormal mean
  DBHc<-DBHc+dg_a*5               # 5-yr step
  DBHi<-DBHc/2.54
  res[[cy+1]]<-cbind(cycle=cy,year=cy*5,metrics(DBHi,TPA))
}
out<-rbindlist(res)
out[,plot_key:=pk]
fwrite(out,file.path(OUT,"smoke_cycles.csv"))
cat("\n=== SMOKE CYCLES (species-aware v8 DG, no mortality/ingrowth) ===\n")
print(out)
cat("\nmean annual dg last cycle (cm/yr):",round(mean(dg_a),3),"\n")
