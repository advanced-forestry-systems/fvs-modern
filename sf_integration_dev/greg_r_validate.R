OUT<-Sys.getenv("GREG_OUT", tempdir())
RDS<-"/users/PUOM0008/crsfaaron/fvs_remodeling/rds"
dg<-as.data.frame(readRDS(file.path(RDS,"dg_parms.RDS")))
hg<-as.data.frame(readRDS(file.path(RDS,"hg_parms.RDS")))
mo<-as.data.frame(readRDS(file.path(RDS,"mort_parm_base_rate_cr_cch.RDS")))
dg<-dg[order(dg$spcd,-dg$isConv),]; dg<-dg[!duplicated(dg$spcd),]
hg<-hg[hg$isConv==TRUE,]; hg<-hg[order(hg$spcd),]; hg<-hg[!duplicated(hg$spcd),]
mo<-mo[order(mo$SPCD,mo$nll),]; mo<-mo[!duplicated(mo$SPCD),]
gd<-function(s){r<-dg[dg$spcd==s,];as.numeric(r[c("B0","B1","B2","B3","B4","B5","B6")])}
gh<-function(s){r<-hg[hg$spcd==s,];as.numeric(r[c("B0","B1","B2","B3","B4","B5","B6","B7","B8")])}
gm<-function(s){r<-mo[mo$SPCD==s,];as.numeric(r[c("b0","b1","b2","b3","b4")])}
dg_annual<-function(s,dbh,cr,ht,bal,elev,emt){b<-gd(s)
  z<-b[1]+b[2]*log((dbh+1)^2/(cr*ht+1)^b[4])+b[3]*bal^b[5]/log(dbh+2.7)+b[6]*elev+b[7]*emt
  z<-min(max(z,-30),5);g<-exp(z);if(!is.finite(g)||g<0)0 else g}
hg_annual<-function(s,ht,cr,ccfl,cch,elev,td,emt){b<-gh(s);mx<-b[1]
  crp<-max(cr,1e-4);cchp<-max(cch,0)
  d<-mx*b[2]*b[3]*crp^b[4]*exp(-b[2]*ht-b[5]*ccfl-b[9]*cchp^0.5-b[6]*elev+b[7]*sqrt(max(td,0))+b[8]*emt)*(1-exp(-b[2]*ht))^(b[3]-1)
  if(!is.finite(d)||d<0)0 else d}
surv_annual<-function(s,cr,cch){b<-gm(s);crp<-max(cr,1e-4);cchp<-max(cch,0)
  eta<-b[1]+b[2]*(crp+0.01)^b[3]+if(cchp>0)b[4]*cchp^b[5] else 0
  eta<-min(max(eta,-30),30);ps<-1-exp(-exp(eta));if(cr<=0)ps<-0;min(max(ps,0),1)}
SP<-c(12,97,316,318,833);DBHS<-c(6,12,20)
CR<-0.5;BAL<-90;CCFL<-90;ELEV<-1200;CCH<-0.4;EMT<--28.8;TD<-24.8
ht_of<-function(dbh)4.5+4.0*dbh
out<-data.frame()
for(s in SP)for(dbh in DBHS){ht<-ht_of(dbh)
  d5<-dbh;for(i in 1:5)d5<-d5+dg_annual(s,d5,CR,ht,BAL,ELEV,EMT)
  out<-rbind(out,data.frame(SPCD=s,DBH=dbh,dg_ann=dg_annual(s,dbh,CR,ht,BAL,ELEV,EMT),
    hg_ann=hg_annual(s,ht,CR,CCFL,CCH,ELEV,TD,EMT),surv_ann=surv_annual(s,CR,CCH),dbh_inc5=d5-dbh))}
write.csv(out,file.path(OUT,"greg_r_grid.csv"),row.names=FALSE)
cat("R grid written\n")
