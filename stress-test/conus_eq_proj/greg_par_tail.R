## ========================== PARALLEL PROJECTION (mclapply) ==========================
## Identical equations/forms/params as the sequential driver. ONLY the stand loop is
## parallelized: each stand is fully independent (NO RNG inside the per-stand loop; the
## only set.seed is the stand-subsample above), so per-stand outputs are deterministic
## regardless of execution order. Tree cache for a state is built ONCE in the parent
## before forking, so workers share it copy-on-write (no re-reads). Results are
## collected and rbind'd in the SAME stand order as the sequential driver.
suppressPackageStartupMessages({ library(parallel) })
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset=NA))
if(is.na(NCORES) || NCORES<1){ ntsk<-suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS",unset=NA))); NCORES<-if(!is.na(ntsk)&&ntsk>1) ntsk else max(1L, detectCores()-1L) }
cat(sprintf("STEP 2: projecting %d stands x %d cycles using mclapply mc.cores=%d ...\n", nrow(si), NCYC, NCORES))

project_one <- function(ix, spl){
  cn<-si$STAND_CN[ix]; trows<-spl[[cn]]; if(is.null(trows)||!nrow(trows)) return(NULL)
  cov<-si[ix]; st<-FIPS[as.character(cov$STATE)]; inv_year<-cov$INV_YEAR
  tl<-mk_tl(trows,cov); if(is.null(tl)||length(tl$dbh_in)<1) return(NULL)
  tl<-recompute_comp(tl)
  elev<-tl$ELEV; emt<-tl$EMT; td<-tl$TD
  tcw<-tl$TPA; tct<-sum(tcw)
  sfb_dg<-sum(tcw[!tl$has_dg])/tct; sfb_hg<-sum(tcw[!tl$has_hg])/tct; sfb_mo<-sum(tcw[!tl$has_mo])/tct
  fb_dg_n<-sum(tcw[!tl$has_dg]); fb_hg_n<-sum(tcw[!tl$has_hg]); fb_mo_n<-sum(tcw[!tl$has_mo])
  fbrow<-data.table(STAND_CN=cn,VARIANT=VARIANT,fb_dg=sfb_dg,fb_hg=sfb_hg,fb_mo=sfb_mo)
  miss<-!is.finite(tl$HT); if(any(miss)) tl$HT[miss]<-4.5*tl$dbh_in[miss]^0.5+4.5
  loc_rows<-vector("list",NCYC+1L); loc_tre<-vector("list",NCYC+1L); k<-0L
  emit<-function(cy){sm<-stand_metrics(tl); py<-cy*CYCLEN; k<<-k+1L
    cchm<-if(!is.null(tl$cch_last)&&length(tl$cch_last)==length(tl$TPA)&&sum(tl$TPA)>0) sum(tl$cch_last*tl$TPA)/sum(tl$TPA) else NA_real_
    loc_rows[[k]]<<-data.table(STAND_CN=cn,STATE=st,YEAR=inv_year+py,PROJ_YEAR=py,VARIANT=VARIANT,CONFIG=CONFIG,AGB_TONS_AC=NA_real_,BA_FT2AC=sm$BA,QMD_IN=sm$QMD,TPH=sm$TPH,CCH_MEAN=cchm)
    loc_tre[[k]]<<-data.table(STAND_CN=cn,CONFIG=CONFIG,PROJ_YEAR=py,SPCD=tl$SPCD,DBH_IN=tl$dbh_in,HT_M=tl$HT*0.3048,TPA=tl$TPA)}
  emit(0)
  for(cy in 1:NCYC){
    gg<-grp_organon(tl$SPCD)
    cch<-stand_cch(tl$dbh_in, tl$HT, tl$CR, tl$TPA, gg); tl$cch_last<-cch
    surv_cyc<-rep(1.0,length(tl$dbh_in))
    for(yr in 1:CYCLEN){
      dgan<-dg_annual(tl$GSPCD, tl$dbh_in, tl$CR, tl$HT, tl$CCFL, elev, emt)
      hgan<-hg_annual(tl$GSPCD_H, tl$HT, tl$CR, tl$CCFL, cch, elev, td, emt)
      psa <-gomp_surv_annual(tl$SPCD, tl$CR, cch)
      tl$dbh_in<-pmin(pmax(tl$dbh_in+dgan,0.1),200)
      tl$HT<-pmin(tl$HT+hgan,400)
      surv_cyc<-surv_cyc*psa
    }
    rd<-tl$rd_add; st_mult<-pmin(1+pmax(rd-1.0,0)*2.0, 2.0)
    H_cyc<- -log(pmin(pmax(surv_cyc,1e-9),1)); surv<-exp(-H_cyc*st_mult); surv[!is.finite(surv)]<-1
    tl$TPA<-tl$TPA*surv; tl$TPA[!is.finite(tl$TPA)]<-0
    BA_m<-tl$BA*0.2296; HT2m<-tl$HT*0.3048
    tl$CR<-cr_update(tl$SPCD,tl$CR,tl$HT*0.3048,HT2m,BA_m,tl$BAL,tl$BA,if(is.na(tl$cspi))1.0 else tl$cspi,tl$cr_zsp)
    keep<-tl$TPA>1e-4
    if(any(!keep)){ vn<-length(tl$dbh_in); for(nm in names(tl)) if(length(tl[[nm]])==vn) tl[[nm]]<-tl[[nm]][keep] }
    tl<-add_recruits(tl, IG_TPA*CYCLEN)
    tl$HT[!is.finite(tl$HT)]<-4.5*tl$dbh_in[!is.finite(tl$HT)]^0.5+4.5
    tl<-recompute_comp(tl); emit(cy)
  }
  list(rows=rbindlist(loc_rows[seq_len(k)]),
       tre =rbindlist(loc_tre[seq_len(k)]),
       fb  =fbrow,
       fb_dg_n=fb_dg_n, fb_hg_n=fb_hg_n, fb_mo_n=fb_mo_n, tc=tct)
}

all_res<-vector("list", length(state_groups)); gi<-0L
for(stcode in names(state_groups)){
  spl<-get_state_trees(as.integer(stcode))
  gi<-gi+1L
  if(is.null(spl)) next
  if(is.list(spl)&&length(spl)==1&&is.null(spl[[1]])) next
  idxs<-state_groups[[stcode]]
  res<-mclapply(idxs, project_one, spl=spl, mc.cores=NCORES, mc.preschedule=TRUE)
  errs<-vapply(res, function(x) inherits(x,"try-error"), logical(1))
  if(any(errs)){ cat("  WORKER ERRORS in state",stcode,":",sum(errs),"-- first:\n"); print(res[[which(errs)[1]]]); stop("mclapply worker error") }
  all_res[[gi]]<-res
}
flat<-unlist(all_res, recursive=FALSE)
flat<-flat[!vapply(flat, is.null, logical(1))]
nproj<-length(flat)
out<-rbindlist(lapply(flat, `[[`, "rows"))
tre<-rbindlist(lapply(flat, `[[`, "tre"))
psf<-rbindlist(lapply(flat, `[[`, "fb"))
fb_dg_acc<-sum(vapply(flat, `[[`, numeric(1), "fb_dg_n"))
fb_hg_acc<-sum(vapply(flat, `[[`, numeric(1), "fb_hg_n"))
fb_mo_acc<-sum(vapply(flat, `[[`, numeric(1), "fb_mo_n"))
tc_acc<-sum(vapply(flat, `[[`, numeric(1), "tc"))

otag<-sprintf("conus_eq_%s_%s",tolower(VARIANT),CONFIG)
fwrite(out,file.path(OUTD,paste0(otag,"_metrics.csv"))); fwrite(tre,file.path(OUTD,paste0(otag,"_treelists.csv")))
fwrite(psf,file.path(OUTD,paste0(otag,"_fallback.csv")))
cat("\n  stands projected:",nproj,"\n")
cat(sprintf("  Greg coverage (tree-count weighted): DG=%.1f%% (fallback %.1f%%)  HG=%.1f%% (fallback %.1f%%)  MORT=%.1f%% (fallback %.1f%%)\n",
  100*(1-fb_dg_acc/tc_acc),100*fb_dg_acc/tc_acc,100*(1-fb_hg_acc/tc_acc),100*fb_hg_acc/tc_acc,100*(1-fb_mo_acc/tc_acc),100*fb_mo_acc/tc_acc))
cat("Wrote:",file.path(OUTD,paste0(otag,"_metrics.csv")),"(",nrow(out),"rows )\n")
cat("Wrote:",file.path(OUTD,paste0(otag,"_treelists.csv")),"(",nrow(tre),"rows )\n")
cat("Wrote:",file.path(OUTD,paste0(otag,"_fallback.csv")),"(",nrow(psf),"stands )\n")
cat("\n=== year 0/50/100 stand-mean metrics (",CONFIG,") ===\n")
print(out[PROJ_YEAR %in% c(0,50,100),.(BA=mean(BA_FT2AC,na.rm=TRUE),QMD=mean(QMD_IN,na.rm=TRUE),TPH=mean(TPH,na.rm=TRUE),n=.N),by=PROJ_YEAR])
cat("DONE.\n")
