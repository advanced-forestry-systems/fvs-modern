## =============================================================================
## cch_module.R  --  Crown Closure at tree tip (CCH) for the gompit mortality arm
##
## Faithful R port of cch_organon.py / validate_cch.R (ORGANON CAL_CCH.for crown
## profile), plus the genus -> ORGANON SWO crown-group crosswalk from
## refine_cch.py and the validated affine map onto the gompit CCH (CCH1) scale.
##
## Pipeline (matches the banked greg_mortality_coefficients.csv, which were fit
## on the panel's stored CCH1):
##   1. per tree: group = grp_organon(SPCD)  (18 ORGANON SWO groups)
##   2. ORGANON crown geometry -> 41-elem crown-area profile -> tree_cch (raw)
##   3. cch1 = CCH_A + CCH_B * cch_hat_raw       (affine, validated R2=0.713,
##      Pearson 0.844 / Spearman 0.925 vs stored CCH1 on 117k-tree sample)
##   4. cch1 clamped to >= 0 ; this is the `cch` covariate Greg's gompit expects.
##
## Inputs imperial: DBH (in), HT (ft), CR (0-1), EXPAN (TPA). The projector keeps
## HT in metres, so convert before calling stand_cch().
## =============================================================================

## ---- ORGANON SWO crown parameters (groups 1..18) ----------------------------
## (identical to cch_organon.py)
.MCWPAR <- list(
 `1`=c(4.6366,1.6078,-0.009625,88.52), `2`=c(6.1880,1.0069,0.0,999.99),
 `3`=c(3.4835,1.343,-0.0082544,81.35), `4`=c(4.6600546,1.0701859,0.0,999.99),
 `5`=c(3.2837,1.2031,-0.0071858,83.71), `6`=c(4.5652,1.4147,0.0,999.99),
 `7`=c(4.0,1.65,0.0,999.99), `8`=c(4.5652,1.4147,0.0,999.99),
 `9`=c(3.4298629,1.3532302,0.0,999.99), `10`=c(2.9793895,1.5512443,-0.01416129,54.77),
 `11`=c(4.4443,1.7040,0.0,999.99), `12`=c(4.4443,1.7040,0.0,999.99),
 `13`=c(4.0953,2.3849,-0.011630,102.53), `14`=c(3.0785639,1.9242211,0.0,999.99),
 `15`=c(3.3625,2.0303,-0.0073307,138.93), `16`=c(8.0,1.53,0.0,999.99),
 `17`=c(2.9793895,1.5512443,-0.01416129,54.77), `18`=c(2.9793895,1.5512443,-0.01416129,54.77))
.LCWPAR <- list(
 `1`=c(0.0,0.00371834,0.808121), `2`=c(0.0,0.00308402,0.0), `3`=c(0.355532,0.0,0.0),
 `4`=c(0.0,0.00339675,0.532418), `5`=c(-0.251389,0.00692512,0.985922), `6`=c(0.0,0.0,0.0),
 `7`=c(-0.251389,0.00692512,0.985922), `8`=c(0.0,0.0,0.0), `9`=c(0.118621,0.00384872,0.0),
 `10`=c(0.0,0.0,1.161440), `11`=c(0.0,0.0111972,0.0), `12`=c(0.0,0.0207676,0.0),
 `13`=c(0.0,0.0,1.47018), `14`=c(0.364811,0.0,0.0), `15`=c(0.0,0.0,1.27196),
 `16`=c(0.3227140,0.0,0.0), `17`=c(0.0,0.0,1.161440), `18`=c(0.0,0.0,1.161440))
.CWAPAR <- list(
 `1`=c(0.929973,-0.135212,-0.0157579), `2`=c(0.999291,0.0,-0.0314603),
 `3`=c(0.755583,0.0,0.0), `4`=c(0.755583,0.0,0.0), `5`=c(0.629785,0.0,0.0),
 `6`=c(0.629785,0.0,0.0), `7`=c(0.629785,0.0,0.0), `8`=c(0.629785,0.0,0.0),
 `9`=c(0.5,0.0,0.0), `10`=c(0.5,0.0,0.0), `11`=c(0.5,0.0,0.0), `12`=c(0.5,0.0,0.0),
 `13`=c(0.5,0.0,0.0), `14`=c(0.5,0.0,0.0), `15`=c(0.5,0.0,0.0), `16`=c(0.5,0.0,0.0),
 `17`=c(0.5,0.0,0.0), `18`=c(0.5,0.0,0.0))
.DACBPAR <- list(`1`=0.062,`2`=0.028454,`3`=0.05,`4`=0.05,`5`=0.20,`6`=0.209806,
 `7`=0.20,`8`=0.209806,`9`=0.0,`10`=0.0,`11`=0.0,`12`=0.0,`13`=0.0,`14`=0.0,
 `15`=0.0,`16`=0.0,`17`=0.0,`18`=0.0)

## ---- affine map onto the gompit CCH (CCH1) scale ----------------------------
## from validate_cch.R on 117k trees / 4000 plots: CCH1 ~ 0.062 + 0.0036*cch_hat
CCH_A <- 0.062; CCH_B <- 0.0036

## ---- FIA SPCD -> ORGANON SWO crown group (refined genus crosswalk) ----------
## ported verbatim from refine_cch.py grp() (MODE=refined)
.WHITE_PINES <- c(101,113,117,119,129)
.LIVE_OAKS   <- c(801,805,807,808,810,843,846,838)
grp_organon <- function(spcd){
  s <- as.integer(spcd)
  g <- integer(length(s))
  for(i in seq_along(s)){
    x <- s[i]
    if(is.na(x)){ g[i]<-14L; next }
    if(x < 300L){                                  # conifers
      if(x==202L || x==201L)            g[i]<-1L
      else if(x>=10L && x<=19L)         g[i]<-2L
      else if(x>=90L && x<=99L)         g[i]<-2L
      else if(x>=260L && x<=269L)       g[i]<-6L
      else if(x %in% c(211L,212L,221L,222L)) g[i]<-6L
      else if(x %in% c(241L,242L))      g[i]<-7L
      else if(x %in% c(231L,251L))      g[i]<-8L
      else if((x>=40L && x<=69L) || x==81L) g[i]<-5L
      else if(x>=70L && x<=79L)         g[i]<-3L
      else if(x>=100L && x<=140L)       g[i]<- if(x %in% .WHITE_PINES) 4L else 3L
      else                              g[i]<-1L
    } else {                                        # hardwoods
      if(x>=310L && x<=329L)            g[i]<-13L
      else if(x>=350L && x<=360L)       g[i]<-16L
      else if(x==361L)                  g[i]<-9L
      else if(x>=370L && x<=379L)       g[i]<-16L
      else if(x==431L)                  g[i]<-10L
      else if(x %in% c(491L,492L))      g[i]<-17L
      else if(x==631L)                  g[i]<-11L
      else if(x>=740L && x<=749L)       g[i]<-16L
      else if(x>=920L && x<=928L)       g[i]<-18L
      else if(x>=800L && x<=849L)       g[i]<- if(x %in% .LIVE_OAKS) 12L else 14L
      else                              g[i]<-14L
    }
  }
  g
}

## ---- ORGANON crown geometry (vectorised-free, scalar like the Fortran) ------
.mcw_f  <- function(g,D,H){p<-.MCWPAR[[as.character(g)]];d<-min(D,p[4]);if(H<4.501) H/4.5*p[1] else p[1]+p[2]*d+p[3]*d*d}
.hlcw_f <- function(g,H,CR){H-(1-.DACBPAR[[as.character(g)]])*CR*H}
.lcw_f  <- function(g,M,CR,D,H){p<-.LCWPAR[[as.character(g)]];CL<-CR*H;M*CR^(p[1]+p[2]*CL+p[3]*(D/H))}
.cw_f   <- function(g,HL,LC,H,D,XL){p<-.CWAPAR[[as.character(g)]];rp<-(H-XL)/(H-HL);if(rp<=0) return(NA_real_);LC*rp^(p[1]+p[2]*sqrt(rp)+p[3]*(H/D))}

## stand_cch: one stand. trees given imperial (DBH in, HT ft, CR 0-1, EXPAN TPA),
## with crown group g. Returns per-tree cch on the gompit (CCH1) scale.
stand_cch <- function(DBH,HT,CR,EXPAN,g){
  n<-length(DBH); if(n<1) return(numeric(0))
  ok0 <- is.finite(DBH)&is.finite(HT)&is.finite(CR)&DBH>0&HT>0&CR>0&CR<=1
  if(!any(ok0)) return(rep(0,n))
  top<-max(HT[ok0]); cch<-numeric(41); cch[41]<-top
  for(i in which(ok0)){
    gi<-g[i];D<-DBH[i];H<-HT[i];cr<-CR[i];E<-EXPAN[i]
    CL<-cr*H;HCB<-H-CL;M<-.mcw_f(gi,D,H);LC<-.lcw_f(gi,M,cr,D,H);HL<-.hlcw_f(gi,H,cr)
    if(!is.finite(LC)||!is.finite(HL)||!is.finite(HCB)) next
    thr<-max(HCB,HL)
    for(ii in 40:1){xl<-(ii-1)*(top/40); cw<-0
      if(xl<=thr){ cw<- if(HCB<=HL) LC else { v<-.cw_f(gi,HL,LC,H,D,max(xl,HCB)); if(is.finite(v)) v else LC } }
      else if(xl<H){ v<-.cw_f(gi,HL,LC,H,D,xl); if(is.finite(v)) cw<-v }
      if(!is.finite(cw)) cw<-0
      cch[ii]<-cch[ii]+(cw^2)*(0.001803*E)}
  }
  ## interpolate raw cch_hat to each tree's tip
  raw<-sapply(seq_len(n),function(k){
    h<-HT[k]
    if(!is.finite(h)||!is.finite(top)||top<=0||h>=top) return(0)
    xi<-40*(h/top); idx<-as.integer(xi)+1
    if(idx>=40) return(cch[40]*(40-xi))
    if(idx<2) idx<-2
    xxi<-(idx+1)-1; v<-cch[idx+1]+(cch[idx]-cch[idx+1])*(xxi-xi)
    if(is.finite(v)) v else 0})
  ## affine map onto the gompit CCH1 scale; clamp >=0
  pmax(CCH_A + CCH_B*raw, 0)
}
