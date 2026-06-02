#!/usr/bin/env Rscript
# Re-fit Greg Johnson's CONUS per-species mortality (Gompit on crown ratio +
# crown closure at tree tip) on the v2 remeasurement panel, which already
# carries CCH1/CCH2 (crown closure at tree tip). Per-species NLL fit with an
# exposure offset for variable interval length. See the model note in
# calibration/R/35_fit_greg_mortality_cch.R.
#
# Usage:
#   Rscript fit_greg_mort_conus.R --data data/conus_remeasurement_pairs_metric_cond_v2.rds \
#       --out /fs/scratch/PUOM0008/crsfaaron/conus_mort/out --min-obs 5000 \
#       --cch-mode start [--subsample 0]

suppressWarnings(suppressMessages(library(data.table)))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d=NULL){ i<-grep(paste0("^--",f,"="),args,value=TRUE); if(!length(i)) return(d); sub(paste0("^--",f,"="),"",i[1]) }
DATA   <- getarg("data", "data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT    <- getarg("out",  "out")
MINOBS <- as.integer(getarg("min-obs","5000"))
CCHMODE<- getarg("cch-mode","start")      # start = CCH1 ; mean = (CCH1+CCH2)/2
SUB    <- as.integer(getarg("subsample","0"))
SLIMOUT<- getarg("slim-out","")     # if set: write slim rds and exit (prep step)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("loading", DATA, "...\n")
d <- as.data.table(readRDS(DATA))
cat("nrow:", nrow(d), " ncol:", ncol(d), "\n")

SLIMCOLS <- c("SPCD","cr","cch","alive","years")
if (all(SLIMCOLS %in% names(d))) {
  # already a slim panel (from the prep step) -> skip mapping/filter
  cat("input is already slim (", paste(SLIMCOLS,collapse=","), ")\n")
  d <- d[, ..SLIMCOLS]
} else {

pick <- function(cands){ h<-intersect(cands, names(d)); if(length(h)) h[1] else NA_character_ }
c_spcd  <- pick(c("SPCD","spcd","SPECIES"))
c_cr    <- pick(c("CR1","cr1","CR","CRATIO1","CRATIO","cr"))
c_cch1  <- pick(c("CCH1","cch1","CCH","cch"))
c_cch2  <- pick(c("CCH2","cch2"))
c_stat  <- pick(c("TREESTATUS2","STATUS2","STATUSCD2","alive","died"))
c_yrs   <- pick(c("YEARS","years","interval","MEAS_INTERVAL","REMPER","remper","INTERVAL"))
cat("columns -> SPCD:",c_spcd," CR:",c_cr," CCH1:",c_cch1," CCH2:",c_cch2,
    " STATUS:",c_stat," YEARS:",c_yrs,"\n")
stopifnot(!is.na(c_spcd), !is.na(c_cr), !is.na(c_cch1), !is.na(c_stat), !is.na(c_yrs))

d[, SPCD := as.integer(get(c_spcd))]
d[, cr   := as.numeric(get(c_cr))]
if (CCHMODE=="mean" && !is.na(c_cch2)) d[, cch := (as.numeric(get(c_cch1))+as.numeric(get(c_cch2)))/2] else d[, cch := as.numeric(get(c_cch1))]
# alive = survived the interval. Map common encodings.
sv <- d[[c_stat]]
if (c_stat %in% c("died")) d[, alive := as.integer(1L - as.integer(sv))]      else
if (c_stat %in% c("alive")) d[, alive := as.integer(sv)]                        else
d[, alive := as.integer(as.integer(sv)==1L)]   # FIA STATUSCD: 1=live
d[, years := as.numeric(get(c_yrs))]

d <- d[is.finite(cr)&is.finite(cch)&is.finite(years)&cr>0&cr<=1&cch>=0&years>=1&years<=20&alive %in% c(0,1)]
cat("after filter:", nrow(d), " surv rate:", round(mean(d$alive),4),
    " cch range [",round(min(d$cch),3),",",round(max(d$cch),3),"]\n")
if (SUB>0 && nrow(d)>SUB) d <- d[sample(.N, SUB)]
# slim to only the columns the fit needs (the source table is 172 cols x 8.2M
# rows; keeping it whole OOMs the mclapply forks). Keep 5 cols, drop the rest.
d <- d[, .(SPCD, cr, cch, alive, years)]; gc()
}  # end mapping/filter branch
cat("slim data:", nrow(d), "rows x", ncol(d), "cols\n")

if (nzchar(SLIMOUT)) {            # prep mode: persist slim panel and stop
  saveRDS(d, SLIMOUT); cat("wrote slim panel ->", SLIMOUT, "\n"); quit(save="no")
}

# Profiled fit: the linear params (b0,b1,b3) are solved exactly each step by a
# cloglog GLM with offset log(years) [cloglog(P_die_T)=b0+b1*X1+b3*X2+log(T)],
# while the two exponents (b2,b4) are searched in 2D. This breaks the b0/b1<->b2
# identifiability ridge a raw 5-param Nelder-Mead falls into, giving stable,
# Greg-comparable coefficients.
# Exponent bounds = Greg's published support. b2 is held away from 0: as b2->0
# the cr term (cr+0.01)^b2 -> 1 (constant), making b0,b1 collinear and blowing
# them up. Restricting b2 to the interior keeps the linear solve identified.
B2_LO <- -1.7; B2_HI <- -0.2
B4_LO <-  0.01; B4_HI <- 3.0
COEF_CAP <- 30          # reject separated/divergent GLM solutions

inner_glm <- function(b2,b4,died,cr,cch,years){
  X1<-(cr+0.01)^b2; X2<-ifelse(cch>0,cch^b4,0)
  fit<-tryCatch(suppressWarnings(glm(died~X1+X2, family=binomial(link="cloglog"),
        offset=log(years), control=list(maxit=50))), error=function(e) NULL)
  if(is.null(fit)||!fit$converged) return(NULL)
  cf<-coef(fit)
  if(any(!is.finite(cf)) || max(abs(cf))>COEF_CAP) return(NULL)   # separation guard
  list(nll=as.numeric(-logLik(fit)), coef=cf)
}
prof_nll <- function(p,died,cr,cch,years){
  r<-inner_glm(p[1],p[2],died,cr,cch,years); if(is.null(r)||!is.finite(r$nll)) 1e12 else r$nll
}
nll_base <- function(cr,cch,alive,years){ optimize(function(b0){H<-exp(min(max(b0,-30),30));HT<-H*years;-sum(ifelse(alive==1,-HT,log1p(-exp(-HT))))},c(-15,5))$objective }
base_intercept <- function(died,years){
  fit<-tryCatch(suppressWarnings(glm(died~1,family=binomial("cloglog"),offset=log(years))),error=function(e)NULL)
  if(is.null(fit)) return(NA_real_); as.numeric(coef(fit)[1])
}
fit_one <- function(dd){
  died<-1L-dd$alive; cr<-dd$cr; cch<-dd$cch; years<-dd$years
  starts<-list(c(-0.6,0.5),c(-1.0,0.3),c(-0.4,1.0)); best<-NULL
  for(s in starts){
    f<-tryCatch(optim(s,prof_nll,died=died,cr=cr,cch=cch,years=years,
          method="L-BFGS-B",lower=c(B2_LO,B4_LO),upper=c(B2_HI,B4_HI),
          control=list(maxit=200)),error=function(e)NULL)
    if(!is.null(f)&&is.finite(f$value)&&f$value<1e11&&(is.null(best)||f$value<best$value)) best<-f
  }
  if(!is.null(best)){
    r<-inner_glm(best$par[1],best$par[2],died,cr,cch,years)
    if(!is.null(r)) return(list(
      par=c(r$coef[["(Intercept)"]], r$coef[["X1"]], best$par[1], r$coef[["X2"]], best$par[2]),
      value=best$value, convergence=best$convergence, fallback=FALSE))
  }
  # fallback: base-rate-only (no cr/cch effect) so the species still has a
  # usable, sane parameter set; flagged so we know it did not take the full form.
  b0<-base_intercept(died,years); if(!is.finite(b0)) return(NULL)
  list(par=c(b0, 0, -0.5, 0, 0.5),
       value=nll_base(cr,cch,1L-died,years), convergence=0L, fallback=TRUE)
}

sp_tab <- sort(table(d$SPCD),decreasing=TRUE); keep<-as.integer(names(sp_tab[sp_tab>=MINOBS]))
cat("species >=",MINOBS,"obs:",length(keep),"\n")
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1"))
setkey(d, SPCD)   # keyed lookup so d[.(sp)] is fast and low-memory per worker
fit_sp <- function(sp){ dd<-d[.(sp)]; if(is.null(dd)||nrow(dd)==0) return(NULL); f<-fit_one(dd); if(is.null(f)) return(NULL)
  bn<-nll_base(dd$cr,dd$cch,dd$alive,dd$years); th<-f$par
  data.table(SPCD=sp,n=nrow(dd),b0=th[1],b1=th[2],b2=th[3],b3=th[4],b4=th[5],
             nll=f$value,nll_baserate=bn,improved=f$value<bn,
             convergence=f$convergence,fallback=f$fallback) }
res <- parallel::mclapply(keep, fit_sp, mc.cores=ncores)
co<-rbindlist(Filter(Negate(is.null), res))
cat("fit", nrow(co), "species;", sum(co$improved), "improved over base rate\n"); fwrite(co,file.path(OUT,"greg_mortality_coefficients.csv"))
fwrite(data.table(n_species=nrow(co),n_improved=sum(co$improved),median_b2=median(co$b2),median_b4=median(co$b4),total_nll=sum(co$nll),total_nll_base=sum(co$nll_baserate)),file.path(OUT,"greg_mortality_fit_summary.csv"))
cat("DONE -> ",file.path(OUT,"greg_mortality_coefficients.csv"),"\n")
