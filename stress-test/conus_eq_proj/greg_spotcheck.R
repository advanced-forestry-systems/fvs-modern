## Reproduce Greg's est_dg, est_hg, survival VERBATIM with his fitted params (sanity reference)
suppressPackageStartupMessages(library(data.table))
RDS<-"/users/PUOM0008/crsfaaron/fvs_remodeling/rds"
dg<-as.data.table(readRDS(file.path(RDS,"dg_parms.RDS")))
hg<-as.data.table(readRDS(file.path(RDS,"hg_parms.RDS")))
mo<-as.data.table(readRDS(file.path(RDS,"mort_parm_base_rate_cr_cch.RDS")))

# --- Greg est_dg (verbatim) ---
est_dg<-function(n,dbh0,bal0,bal1,cr0,cr1,ht0,ht1,elev,emt,B0,B1,B2,B3,B4,B5,B6){
  dht<-(ht1-ht0)/n; dbal<-(bal1-bal0)/n; dcr<-(cr1-cr0)/n
  cht<-ht0; cdbh<-dbh0; cbal<-bal0; ccr<-cr0; max_n<-max(n)
  for(i in 1:max_n){
    dg_hat<-exp(B0+B1*log((cdbh+1)^2/(ccr*cht+1.0)^B3)+B2*cbal^B4/(log(cdbh+2.7))+B5*elev+B6*emt)
    cdbh<-cdbh+ifelse(i<=n,dg_hat,0); cbal<-cbal+ifelse(i<=n,dbal,0); cht<-cht+ifelse(i<=n,dht,0); ccr<-ccr+ifelse(i<=n,dcr,0)}
  cdbh}
# --- Greg est_hg (verbatim) ---
est_hg<-function(periods,ht,cr,cr2,ccfl,ccfl2,cch,cch2,max_height,elev,td,emt,b1,b2,b3,b4,b5,b6,b7,b8){
  htc<-ht; crc<-cr; cchc<-cch; ccflc<-ccfl
  crg<-(cr2-cr)/periods; ccflg<-(ccfl2-ccfl)/periods; cchg<-(cch2-cch)/periods
  for(i in 1:max(periods)){
    htc<-htc+ifelse(i<=periods, max_height*b1*b2*(crc)^b3*exp(-b1*htc-b4*ccflc-b8*cchc^0.5-b5*elev+b6*td^0.5+b7*emt)*(1.0-exp(-b1*htc))^(b2-1.0),0.0)
    crc<-crc+ifelse(i<=periods,crg,0); ccflc<-ccflc+ifelse(i<=periods,ccflg,0); cchc<-cchc+ifelse(i<=periods,cchg,0)}
  htc}
# --- Greg survival cr+cch (verbatim per-step, annual loop) ---
surv<-function(b,periods,startCR,endCR,startCCH,endCCH){
  ccr<-startCR; dcr<-(endCR-startCR)/periods; ccch<-startCCH; dcch<-(endCCH-startCCH)/periods
  p<-rep(1.0,length(periods))
  for(i in 1:max(periods)){
    pstep<-ifelse(ccr>0, 1-exp(-exp(b[1]+b[2]*(ccr+0.01)^b[3]+b[4]*ccch^b[5])), 0)
    m<-i<=periods; p[m]<-p[m]*pstep[m]; ccr<-ccr+ifelse(m,dcr,0); ccch<-ccch+ifelse(m,dcch,0)}
  p}

spp<-c(131,202,316)  # loblolly, doug-fir, red maple
nyr<-5
cat("==== DG (5yr, dbh0=10in bal=80 cr=0.5 ht=50ft elev=1500 EMT=-20; covariates held constant) ====\n")
for(s in spp){ p<-dg[spcd==s]
  d5<-est_dg(nyr,10,80,80,0.5,0.5,50,50,1500,-20,p$B0,p$B1,p$B2,p$B3,p$B4,p$B5,p$B6)
  cat(sprintf("  SPCD %d %-14s 5yr dbh: 10.000 -> %.3f in  (%.3f in/yr)\n",s,p$Common_Name,d5,(d5-10)/nyr)) }
cat("\n==== HG (5yr, ht=50ft cr=0.5 ccfl=100 cch=0.5 elev=1500 TD=17 EMT=-20; held constant) ====\n")
for(s in spp){ p<-hg[spcd==s]
  if(!nrow(p)){cat(sprintf("  SPCD %d -- no HG params\n",s));next}
  h5<-est_hg(nyr,50,0.5,0.5,100,100,0.5,0.5,p$B0,1500,17,-20,p$B1,p$B2,p$B3,p$B4,p$B5,p$B6,p$B7,p$B8)
  cat(sprintf("  SPCD %d %-14s max_ht(B0)=%.0f  5yr ht: 50.0 -> %.2f ft (%.3f ft/yr)\n",s,p$Common_Name,p$B0,h5,(h5-50)/nyr)) }
cat("\n==== Mortality cr+cch (annual survival; cr=0.5 cch=0.5 constant) ====\n")
for(s in spp){ p<-mo[SPCD==s]
  if(!nrow(p)){cat(sprintf("  SPCD %d -- no mort params\n",s));next}
  s1<-surv(c(p$b0,p$b1,p$b2,p$b3,p$b4),1,0.5,0.5,0.5,0.5)
  s5<-surv(c(p$b0,p$b1,p$b2,p$b3,p$b4),5,0.5,0.5,0.5,0.5)
  cat(sprintf("  SPCD %d %-14s annual P_surv=%.5f  (annual mort=%.4f) 5yr surv=%.4f\n",s,p$CommonName,s1,1-s1,s5)) }
