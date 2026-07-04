#!/usr/bin/env Rscript
# Cross-check the Fortran DG/HG evaluators against the projector's R forms at the
# same tree state, and sanity-check that the increments are biologically plausible.
co_dg <- read.csv("greg_dg_coefficients.csv"); fo_dg <- read.csv("dg_fortran.csv")
co_hg <- read.csv("greg_hg_coefficients.csv"); fo_hg <- read.csv("hg_fortran.csv")
dbh<-8; cr<-0.5; ht<-50; bal<-80; ccfl<-120; cch<-0.4; elev<-1500; td<-25; emt<- -15

dgf <- function(B0,B1,B2,B3,B4,B5,B6){
  z <- B0 + B1*log((dbh+1)^2/(cr*ht+1)^B3) + B2*bal^B4/log(dbh+2.7) + B5*elev + B6*emt
  z <- pmin(pmax(z,-30),5); pmax(exp(z),0)
}
co_dg$dg_r <- with(co_dg, dgf(B0,B1,B2,B3,B4,B5,B6))
m <- merge(fo_dg, co_dg[,c("SPCD","dg_r")], by="SPCD"); m$ad <- abs(m$dg - m$dg_r)

hgf <- function(mx,b1,b2,b3,b4,b5,b6,b7,b8){
  crp<-max(cr,1e-4); cchp<-max(cch,0); tdp<-max(td,0)
  mx*b1*b2*crp^b3*exp(-b1*ht -b4*ccfl -b8*cchp^0.5 -b5*elev +b6*sqrt(tdp) +b7*emt)*(1-exp(-b1*ht))^(b2-1)
}
co_hg$hg_r <- pmax(with(co_hg, mapply(hgf,B0,B1,B2,B3,B4,B5,B6,B7,B8)), 0)
m2 <- merge(fo_hg, co_hg[,c("SPCD","hg_r")], by="SPCD"); m2$ad <- abs(m2$hg - m2$hg_r)

cat(sprintf("DG: n=%d  max|Fortran-R|=%.2e  dg median=%.4f in/yr  range [%.4f, %.4f]\n",
            nrow(m), max(m$ad), median(m$dg_r), min(m$dg_r), max(m$dg_r)))
cat(sprintf("HG: n=%d  max|Fortran-R|=%.2e  hg median=%.4f ft/yr  range [%.4f, %.4f]\n",
            nrow(m2), max(m2$ad), median(m2$hg_r), min(m2$hg_r), max(m2$hg_r)))
match_ok <- max(m$ad) < 1e-4 && max(m2$ad) < 1e-4
bio_ok <- median(m$dg_r) > 0.02 && median(m$dg_r) < 1.0 && median(m2$hg_r) > 0.05 && median(m2$hg_r) < 3.0
cat(if (match_ok) "PASS match: Fortran reproduces the projector DG/HG forms\n" else "FAIL match\n")
cat(if (bio_ok) "PASS plausibility: DG/HG increments in expected range\n" else "CHECK plausibility: increments outside typical range (inspect)\n")
