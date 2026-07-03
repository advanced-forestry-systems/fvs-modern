prep_dg <- function(dat){
  # kuehne_v8 diameter-growth predictor (log-normal, cm/yr). Covariates other
  # than the species term are shared across all arms, so any minor knot /
  # softwood-centering approximation cancels in the pure_sf vs hybrid vs legA
  # comparison. Matches v8_holdout_predict.R eta construction (main-effect part;
  # the z_L1_bgi and gamma_site BGI-slope refinements are RE/interaction terms
  # not carried in this shared-base benchmark, consistent with the other five
  # components which also use point-estimate fixed effects + standard REs).
  dat[, resp:=(DBH2-DBH1)/YEARS]; dat[, sqrt_years:=sqrt(YEARS)]
  dat[, ln_dbh:=log(DBH1)]; dat[, ln_cr_adj:=log((CR1+0.2)/1.2)]
  dat[, ln_bal_sw_adj:=log(BAL_SW1+0.01)]
  dat[, rd_additive:=sdi_additive1/SDImax_brms]
  dat[, sdi_complexity:=sdi_additive1/pmax(SDI1,1.0)]
  dat[, ba_metric:=BA1*0.2296]
  dat[, ba_x_rd:=ba_metric*rd_additive]; dat[, balsw_x_rd:=BAL_SW1*rd_additive]
  dat[, bgi_rd:=bgi*rd_additive]; dat[, bgi_lndbh:=bgi*ln_dbh]; dat[, bgi_lncradj:=bgi*ln_cr_adj]
  # softwood centering from the (globally loaded) traits table; species-level
  sw <- setNames(traits$softwood, as.integer(traits$SPCD))
  swv <- ifelse(as.integer(dat$SPCD) %in% names(sw), sw[as.character(dat$SPCD)], 0)
  swv[is.na(swv)] <- 0
  sw_mean <- mean(unique(data.table(SPCD=as.integer(dat$SPCD), sw=as.numeric(swv)))$sw, na.rm=TRUE)
  dat[, softwood_c:=as.numeric(swv)-sw_mean]; dat[, bgi_softwood_c:=bgi*softwood_c]
  # 3-piece BGI basis: knots at 25th/75th percentile of bgi (per stan comment)
  k1 <- as.numeric(quantile(dat$bgi,0.25,na.rm=TRUE)); k2 <- as.numeric(quantile(dat$bgi,0.75,na.rm=TRUE))
  dat[, bgi_b2:=pmax(bgi-k1,0)]; dat[, bgi_b3:=pmax(bgi-k2,0)]
  dat <- dat[is.finite(DBH1)&DBH1>=2.54 & is.finite(DBH2) & is.finite(CR1)&CR1>0&CR1<=1 &
    is.finite(YEARS)&YEARS>=1&YEARS<=20 & TREESTATUS1==1 & TREESTATUS2==1 &
    is.finite(BAL_SW1)&BAL_SW1>=0 & is.finite(BAL_HW1)&BAL_HW1>=0 &
    is.finite(rd_additive)&rd_additive>0&rd_additive<3.0 &
    is.finite(sdi_complexity)&sdi_complexity>0&sdi_complexity<10 &
    is.finite(BA1)&BA1>=0 & is.finite(bgi) &
    !is.na(EPA_L1_CODE)&EPA_L1_CODE!="" & !is.na(EPA_L2_CODE)&EPA_L2_CODE!="" &
    !is.na(EPA_L3_CODE)&EPA_L3_CODE!="" & !is.na(FORTYPCD_cond1)&FORTYPCD_cond1>0 &
    resp>0.01 & resp<5.0]
  gk <- function(f,k) if(k %in% names(f)) f[[k]] else 0
  list(dat=dat, gfun=function(y)log(y), hfun=function(e)exp(pmin(e,20)),
       valid=function(y)y>0.01,
       cov=function(d,f)
         gk(f,"b1")*d$ln_dbh + gk(f,"b2")*d$DBH1 + gk(f,"b3")*d$ln_cr_adj +
         gk(f,"b4")*d$ln_bal_sw_adj + gk(f,"b5")*d$BAL_HW1 +
         gk(f,"b6")*d$bgi + gk(f,"b9a")*d$bgi_b2 + gk(f,"b9b")*d$bgi_b3 +
         gk(f,"b7")*d$ba_x_rd + gk(f,"b8")*d$balsw_x_rd + gk(f,"b11")*d$sdi_complexity +
         gk(f,"b12")*d$bgi_rd + gk(f,"b13")*d$bgi_lndbh + gk(f,"b14")*d$bgi_softwood_c +
         gk(f,"b15")*d$bgi_lncradj,
       interval=function(e,d,f){ s<-f[["sigma"]]/d$sqrt_years
         list(l=exp(pmin(e-1.96*s,20)), u=exp(pmin(e+1.96*s,20))) })
}
