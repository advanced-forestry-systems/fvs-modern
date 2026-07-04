co <- read.csv("/fs/scratch/PUOM0008/crsfaaron/track_mortverify/greg_mortality_coefficients_cch2_regen.csv")
fo <- read.csv("/fs/scratch/PUOM0008/crsfaaron/track_mortverify/surv_fortran_current.csv")
g <- function(b0,b1,b2,b3,b4,cr,cch,fint){
  crc<-pmin(pmax(cr,1e-4),1); cchc<-pmax(cch,0)
  cterm<-ifelse(cchc>0, cchc^b4, 0)
  eta<-b0+b1*(crc+0.01)^b2+b3*cterm; eta<-pmin(pmax(eta,-30),30)
  hz<-1-exp(-exp(eta)); pmin(pmax(hz,0),1)^fint
}
m<-merge(fo,co,by="SPCD")
m$surv_r<-with(m,g(b0,b1,b2,b3,b4,cr,cch,fint))
m$absdiff<-abs(m$surv-m$surv_r)
cat(sprintf("rows=%d species=%d\n",nrow(m),length(unique(m$SPCD))))
cat(sprintf("max abs diff (Fortran f32 vs R f64): %.3e  mean: %.3e\n",max(m$absdiff),mean(m$absdiff)))
cat(sprintf("survival range: [%.5f, %.5f]\n",min(m$surv_r),max(m$surv_r)))
cat(if(max(m$absdiff)<1e-6) "PASS <1e-6\n" else if(max(m$absdiff)<1e-4) "PASS <1e-4 (single-prec float)\n" else "FAIL\n")
print(head(m[order(-m$absdiff),c("SPCD","cr","cch","surv","surv_r","absdiff")],4),row.names=FALSE)
