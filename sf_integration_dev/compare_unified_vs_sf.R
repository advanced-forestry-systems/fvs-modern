##=============================================================================
## compare_unified_vs_sf.R
## Unified (species-free hierarchy + species RE z_sp) vs species-free (Leg B).
##   (1) LOO compare (does z_sp improve held-in predictive density?)
##   (2) RE-scale decomposition: does adding sigma_sp shrink the ecoregion /
##       forest-type variances? (the site-productivity confounding question)
## Matched rows: both fits use same data, seed 42, subsample, min_sp -> log_lik
## columns align 1:1 and loo_compare is valid.
##=============================================================================
suppressMessages({library(data.table); library(loo); library(cmdstanr); library(posterior)})
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d=NULL){m<-grep(paste0("^--",n,"="),args,value=TRUE); if(!length(m)) return(d); sub(paste0("^--",n,"="),"",m[1])}
COMP<-ga("comp","hcb"); UNI_LOO<-ga("unified_loo"); UNI_SUMM<-ga("unified_summary")
SF_FIT<-ga("sf_fit"); SF_LOO<-ga("sf_loo"); SF_SUMM<-ga("sf_summary"); OUT<-ga("out",paste0("compare_",COMP))
cat("== compare_unified_vs_sf.R ==  comp:", COMP, "\n")
uni_loo <- readRDS(UNI_LOO)
if (!is.null(SF_LOO) && file.exists(SF_LOO)) { sf_loo <- readRDS(SF_LOO) } else {
  cat("Computing species-free LOO from fit (memory heavy) ...\n"); flush.console()
  sf_fit <- readRDS(SF_FIT); ll <- sf_fit$draws("log_lik", format="draws_matrix")
  sf_loo <- loo::loo(ll); rm(ll, sf_fit); gc(); saveRDS(sf_loo, paste0(OUT,"_sf_loo.rds"))
}
n_uni<-dim(uni_loo$pointwise)[1]; n_sf<-dim(sf_loo$pointwise)[1]
cat(sprintf("N obs: unified=%d  species_free=%d  %s\n", n_uni, n_sf, ifelse(n_uni==n_sf,"(matched)","(MISMATCH - loo_compare invalid!)")))
cmp <- loo::loo_compare(list(unified=uni_loo, species_free=sf_loo))
cat("\n=== LOO compare (positive elpd_diff favors top row) ===\n"); print(cmp)
fwrite(as.data.table(cmp, keep.rownames="model"), paste0(OUT,"_loo_compare.csv"))
get_sig <- function(f){s<-fread(f); s[grepl("^sigma|^phi",variable), .(variable,mean,q5,q95,rhat)]}
su<-get_sig(UNI_SUMM); setnames(su,c("mean","q5","q95","rhat"),c("uni_mean","uni_q5","uni_q95","uni_rhat"))
ss<-get_sig(SF_SUMM); setnames(ss,c("mean","q5","q95","rhat"),c("sf_mean","sf_q5","sf_q95","sf_rhat"))
m<-merge(ss,su,by="variable",all=TRUE); m[,pct_change:=round(100*(uni_mean-sf_mean)/sf_mean,1)]
setcolorder(m,c("variable","sf_mean","uni_mean","pct_change","sf_rhat","uni_rhat"))
cat("\n=== RE-scale decomposition (species-free -> unified; neg pct = shrinkage) ===\n"); print(m)
fwrite(m, paste0(OUT,"_sigma_decomp.csv"))
cat("\nWrote:", paste0(OUT, c("_loo_compare.csv","_sigma_decomp.csv")), "\n")
