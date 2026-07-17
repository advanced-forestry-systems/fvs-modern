suppressPackageStartupMessages({library(data.table)})
GREG_RDS<-"/users/PUOM0008/crsfaaron/fvs_remodeling/rds"
dgP<-as.data.table(readRDS(file.path(GREG_RDS,"dg_parms.RDS")))
dgP<-dgP[order(spcd,-isConv)][!duplicated(spcd)]
DG<-as.list(setNames(lapply(c("B0","B1","B2","B3","B4","B5","B6"),function(c) setNames(dgP[[c]],as.character(dgP$spcd))),c("B0","B1","B2","B3","B4","B5","B6")))
dg_annual<-function(SPCD,dbh,cr,ht,bal,elev,emt){
  k<-as.character(SPCD)
  B0<-DG$B0[k];B1<-DG$B1[k];B2<-DG$B2[k];B3<-DG$B3[k];B4<-DG$B4[k];B5<-DG$B5[k];B6<-DG$B6[k]
  z<- B0 + B1*log((dbh+1)^2/(cr*ht+1.0)^B3) + B2*bal^B4/log(dbh+2.7) + B5*elev + B6*emt
  z<-pmin(pmax(z,-30),5); g<-exp(z); g[!is.finite(g)]<-0; pmax(g,0)
}
cat("=== DG param coefficient ranges (84 spp) ===\n")
for(c in c("B0","B1","B2","B3","B4","B5","B6")) cat(sprintf("  %s: min=%.4g med=%.4g max=%.4g\n",c,min(dgP[[c]],na.rm=T),median(dgP[[c]],na.rm=T),max(dgP[[c]],na.rm=T)))
# common NE species: 316 red maple, 833 red oak, 12 balsam fir, 129 e white pine, 371 yellow birch, 531 am beech
sp<-c(316,833,12,129,371,531,318,241,95,97)
sp<-sp[as.character(sp) %in% names(DG$B0)]
cat("\n=== annual DG (in/yr) at dbh=9,cr=0.5,ht=50ft,bal=80,elev=1000,emt=6 ===\n")
for(s in sp){ g<-dg_annual(s,9,0.5,50,80,1000,6); cat(sprintf("  SPCD %d: %.4f in/yr  (5yr cycle=%.3f in)\n",s,g,g*5)) }
cat("\n=== DG vs dbh (SPCD=316 red maple), cr=0.5 ht=50 bal=80 elev=1000 emt=6 ===\n")
for(d in c(4,6,9,12,16,20)) cat(sprintf("  dbh=%2d -> %.4f in/yr\n",d,dg_annual(316,d,0.5,50,80,1000,6)))
cat("\n=== DG at LOW competition (bal=10) vs HIGH (bal=200), SPCD=316 dbh=9 ===\n")
cat(sprintf("  bal=10 -> %.4f in/yr ; bal=200 -> %.4f in/yr\n",dg_annual(316,9,0.5,50,10,1000,6),dg_annual(316,9,0.5,50,200,1000,6)))
# median annual DG across all species at a nominal state
gall<-sapply(dgP$spcd,function(s) dg_annual(s,9,0.5,50,80,1000,6))
cat(sprintf("\n=== across all 84 spp at nominal state: median=%.4f mean=%.4f in/yr; frac<0.02=%.2f ===\n",
  median(gall,na.rm=T),mean(gall,na.rm=T),mean(gall<0.02,na.rm=T)))
