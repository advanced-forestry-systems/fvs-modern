##=============================================================================
## 31_fit_dg_organon_v2.R  -- TRACTABLE threaded fit (DO NOT replace v1)
## Same ORGANON-form (Hann SWO) DG model, species x ecodivision hierarchy,
## K1/K2 estimated, DIRECT PAI likelihood -- threaded via reduce_sum +
## stratified subsample + fewer iters so the CONUS fit finishes in hours and
## yields per-species/ecodivision POINT ESTIMATES (posterior means) for arm-6.
## Saves dg_organon_cspi_traits1_fit.rds (extractor-compatible) plus compact
## fixed_summary / species_intercepts / ecodiv_intercepts CSVs directly.
##=============================================================================
library(tidyverse); library(cmdstanr); library(posterior)
K4 <- 2.7; MIN_OBS_SPECIES <- 5000
STAN_THREADED <- "calibration/stan/organon_dg_conus_threaded.stan"
OUT_DIR <- "calibration/output/conus"; dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

prepare_dg_data <- function(dat, site_var="climate_si") {
  message("Preparing DG data... site=", site_var)
  dat <- dat %>% mutate(DBH1_cm=DBH1*2.54, DBH2_cm=DBH2*2.54, HT1_m=HT1*0.3048,
    BA1_metric=BA1*0.2296, BAL1_metric=BAL1*0.2296, QMD1_metric=QMD1*2.54, BA2_metric=BA2*0.2296, BAL2_metric=BAL2*0.2296)
  dat <- dat %>% mutate(dg_obs=(DBH2_cm-DBH1_cm), ln_cr_adj=log((CR1+0.2)/1.2),
    ln_site_prod=log(pmax(.data[[site_var]],1.0)), bal_comp=BAL1_metric/log(DBH1_cm+K4),
    sqrt_ba=sqrt(BA1_metric), rd=BA1_metric/sqrt(QMD1_metric),
    ln_bal=log(BAL1_metric+5.0), ln_bal_5=log(BAL1_metric+5.0))
  dat <- dat %>% filter(dg_obs > -1.27, dg_obs < 2.54*5.0*YEARS, DBH1_cm>=2.54,
    CR1>0 & CR1<=1.0, .data[[site_var]]>0, BA1_metric>0, QMD1_metric>0,
    YEARS>=1 & YEARS<=20)
  sp_counts <- dat %>% count(SPCD) %>% filter(n>=MIN_OBS_SPECIES)
  dat <- dat %>% filter(SPCD %in% sp_counts$SPCD)
  sp_levels <- sort(unique(dat$SPCD)); eco_levels <- sort(unique(dat$ecodiv_code))
  dat <- dat %>% mutate(species_idx=match(SPCD,sp_levels), ecodiv_idx=match(ecodiv_code,eco_levels))
  message("  Trees: ", nrow(dat), " | Species: ", length(sp_levels),
          " | Ecodivisions: ", length(eco_levels), " | Mean interval: ", round(mean(dat$YEARS),1)," yr")
  list(data=dat, species=sp_levels, ecodiv=eco_levels)
}

stratified_subsample <- function(dat, n_target, seed=42, min_per=30) {
  set.seed(seed); if (nrow(dat) <= n_target) return(dat)
  dat$.strat <- paste(dat$species_idx, dat$ecodiv_idx, sep="_")
  strata <- split(seq_len(nrow(dat)), dat$.strat)
  floor_idx <- unlist(lapply(strata, function(ix) if (length(ix)<=min_per) ix else sample(ix,min_per)), use.names=FALSE)
  floor_idx <- unique(floor_idx)
  remaining <- setdiff(seq_len(nrow(dat)), floor_idx)
  n_extra <- max(0, n_target-length(floor_idx))
  extra_idx <- if (n_extra>0 && length(remaining)>0) sample(remaining, min(n_extra,length(remaining))) else integer(0)
  keep <- sort(c(floor_idx, extra_idx)); out <- dat[keep,,drop=FALSE]; out$.strat <- NULL
  message("  Stratified subsample: ", nrow(out), " rows across ", length(strata),
          " species x ecodiv strata (floor ", min_per, "/stratum)"); out
}

make_stan_data_threaded <- function(prep, grainsize) {
  dat <- prep$data
  list(N=nrow(dat), N_species=length(prep$species), N_ecodiv=length(prep$ecodiv),
    grainsize=as.integer(grainsize), dg_obs=dat$dg_obs, dbh=dat$DBH1_cm,
    ln_cr_adj=dat$ln_cr_adj, ln_site_prod=dat$ln_site_prod,
    bal_ratio=dat$BAL1_metric/log(dat$DBH1_cm+2.7), ln_bal=dat$ln_bal, sqrt_ba=dat$sqrt_ba,
    clim1=if("clim_pca1"%in%names(dat))dat$clim_pca1 else rep(0,nrow(dat)),
    clim2=if("clim_pca2"%in%names(dat))dat$clim_pca2 else rep(0,nrow(dat)),
    rd=dat$rd, years=dat$YEARS, species_id=dat$species_idx, ecodiv_id=dat$ecodiv_idx)
}

init_fn_factory <- function(N_species, N_ecodiv, seed) function(chain_id) {
  set.seed(seed+chain_id)
  list(mu_b0=rnorm(1,-2.0,0.3), sigma_sp=0.30, sigma_eco=0.30,
    z_sp=rnorm(N_species,0,0.10), z_eco=rnorm(N_ecodiv,0,0.10), K1=1.0, K2=0.8,
    b1=rnorm(1,0.40,0.10), b2=rnorm(1,-0.02,0.01), b3=rnorm(1,0.80,0.10),
    b4=rnorm(1,0.30,0.10), b5=rnorm(1,-0.005,0.002), b6=rnorm(1,-0.03,0.01),
    b7=rnorm(1,0,0.005), b8=rnorm(1,0,0.005), b9=rnorm(1,0,0.005),
    b10=rnorm(1,0,0.005), sigma=0.5)
}

predict_dg_annualized <- function(params, dat) {
  b0_total<-params$b0_total; b1<-params$b1; b2<-params$b2; b3<-params$b3
  b4<-params$b4; b5<-params$b5; b6<-params$b6; b7<-params$b7; b8<-params$b8
  b9<-params$b9; b10<-params$b10; K1<-params$K1; K2<-params$K2
  N<-nrow(dat); max_years<-max(dat$YEARS)
  d_curr<-dat$DBH1_cm; bal_curr<-dat$BAL1_metric; cr_curr<-dat$CR1
  ba_curr<-dat$BA1_metric; qmd_curr<-dat$QMD1_metric
  bal_rate<-(dat$BAL2_metric-dat$BAL1_metric)/dat$YEARS
  cr_rate<-(dat$CR2-dat$CR1)/dat$YEARS; ba_rate<-(dat$BA2_metric-dat$BA1_metric)/dat$YEARS
  clim1<-if("clim_pca1"%in%names(dat))dat$clim_pca1 else rep(0,N)
  clim2<-if("clim_pca2"%in%names(dat))dat$clim_pca2 else rep(0,N)
  for (t in seq_len(max_years)) {
    rd_curr<-ba_curr/sqrt(qmd_curr); ln_bal_curr<-log(bal_curr+5.0)
    ln_dg<-b0_total + b1*log(d_curr+K1) + b2*d_curr^K2 + b3*log((cr_curr+0.2)/1.2) +
      b4*dat$ln_site_prod + b5*bal_curr/log(d_curr+K4) + b6*sqrt(ba_curr) +
      b7*clim1 + b8*clim2 + b9*rd_curr + b10*rd_curr*ln_bal_curr
    dg_annual<-exp(pmin(pmax(ln_dg,-30),20))
    cradj<-ifelse(cr_curr<=0.17, 1.0-exp(-(25.0*cr_curr)^2), 1.0); dg_annual<-dg_annual*cradj
    d_curr<-d_curr+ifelse(t<=dat$YEARS, dg_annual, 0.0)
    bal_curr<-bal_curr+bal_rate; cr_curr<-cr_curr+cr_rate; ba_curr<-ba_curr+ba_rate
  }
  d_curr - dat$DBH1_cm
}

point_means <- function(fit, vars) {
  if (inherits(fit,"CmdStanMLE")) { s<-fit$summary(variables=vars); setNames(s$estimate,s$variable) }
  else { s<-fit$summary(variables=vars); setNames(s$mean,s$variable) } }
point_means_indexed <- function(fit, base) { d<-fit$draws(variables=base, format="draws_matrix"); apply(d,2,mean) }

if (sys.nframe()==0) {
  args <- commandArgs(trailingOnly=TRUE)
  av <- function(flag, default=NULL, coerce=identity) { i<-which(args==flag)
    if (length(i)>0 && i+1<=length(args)) coerce(args[i+1]) else default }
  site_var <- if ("--bgi"%in%args) "bgi" else if ("--site"%in%args) args[which(args=="--site")+1] else "climate_si"
  engine   <- av("--engine","sample")
  n_sub    <- av("--n",180000L,as.integer); min_per <- av("--min_per",30L,as.integer)
  threads  <- av("--threads_per_chain",12L,as.integer); grainsz <- av("--grainsize",0L,as.integer)
  max_td   <- av("--max_treedepth",10L,as.integer); adapt_d <- av("--adapt_delta",0.85,as.numeric)
  n_chain  <- av("--chains",4L,as.integer); n_warm <- av("--warmup",500L,as.integer)
  n_samp   <- av("--sampling",500L,as.integer); rseed <- av("--seed",42L,as.integer)
  data_file<- av("--data","calibration/data/conus_remeasurement_pairs.rds")
  out_dir  <- av("--out","calibration/output/conus/dg"); stan_file<- av("--stan_file",STAN_THREADED)
  holdout  <- av("--holdout",50000L,as.integer); model_id <- av("--model_id","dg_organon_cspi_traits1")
  message("=== ORGANON DG CONUS v2 (tractable/threaded) ===")
  message("engine=",engine," site=",site_var," n=",n_sub," chains=",n_chain," threads/chain=",threads,
          " warmup/samp=",n_warm,"/",n_samp," adapt_delta=",adapt_d," max_td=",max_td)
  message("stan=",stan_file); message("data=",data_file); message("out=",out_dir)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  if (!file.exists(data_file)) stop("Data file not found: ", data_file)
  raw <- readRDS(data_file); prep_full <- prepare_dg_data(raw, site_var=site_var)
  set.seed(rseed+999)
  ho_idx <- if (holdout>0 && nrow(prep_full$data)>holdout) sample(nrow(prep_full$data),holdout) else integer(0)
  holdout_dat <- if (length(ho_idx)>0) prep_full$data[ho_idx,] else prep_full$data
  train_dat   <- if (length(ho_idx)>0) prep_full$data[-ho_idx,] else prep_full$data
  prep <- prep_full
  prep$data <- if (n_sub>0L) stratified_subsample(train_dat,n_sub,rseed,min_per) else train_dat
  message("Fitting on N=", nrow(prep$data))
  if (grainsz<=0L) grainsz <- max(50L, as.integer(ceiling(nrow(prep$data)/(n_chain*threads*5))))
  message("grainsize=", grainsz)
  mod <- cmdstan_model(stan_file, cpp_options=list(stan_threads=TRUE))
  stan_data <- make_stan_data_threaded(prep, grainsz)
  initf <- init_fn_factory(stan_data$N_species, stan_data$N_ecodiv, rseed)
  t0 <- Sys.time()
  if (engine=="sample") {
    fit <- mod$sample(data=stan_data, chains=n_chain, parallel_chains=n_chain,
      threads_per_chain=threads, iter_warmup=n_warm, iter_sampling=n_samp, seed=rseed,
      max_treedepth=max_td, adapt_delta=adapt_d, init=lapply(seq_len(n_chain),initf), refresh=50)
    method_used <- sprintf("MCMC reduce_sum (%d chains x %d threads, %d/%d warmup/samp, N=%d stratified)",
                           n_chain, threads, n_warm, n_samp, stan_data$N)
  } else if (engine=="optimize") {
    fit <- mod$optimize(data=stan_data, threads=threads, seed=rseed, init=list(initf(1)), jacobian=TRUE)
    method_used <- sprintf("Penalized MAP optimize (threads=%d, N=%d stratified)", threads, stan_data$N)
  } else if (engine=="variational") {
    fit <- mod$variational(data=stan_data, threads=threads, seed=rseed, init=list(initf(1)))
    method_used <- sprintf("ADVI variational (threads=%d, N=%d stratified)", threads, stan_data$N)
  } else stop("Unknown --engine: ", engine)
  elapsed <- as.numeric(difftime(Sys.time(),t0,units="mins"))
  message(sprintf("Fit wall time: %.1f min", elapsed))
  scalar_vars <- c("mu_b0", paste0("b",1:10), "K1","K2","sigma_sp","sigma_eco","sigma")
  if (engine=="sample") {
    summ <- fit$summary(variables=scalar_vars); print(summ)
    rb <- summ %>% filter(rhat>1.05); message("  scalars Rhat>1.05: ", nrow(rb))
    diag <- fit$diagnostic_summary()
    message("  Divergent transitions: ", sum(diag$num_divergent),
            " | max_treedepth hits: ", sum(diag$num_max_treedepth))
  } else if (engine=="optimize") { message("  MAP lp__: ", round(fit$lp(),1))
  } else { message("  ADVI ELBO in fit$output()") }
  fit$save_object(file.path(out_dir, sprintf("%s_fit.rds", model_id)))
  saveRDS(list(site_var=site_var, method=method_used,
    prep_meta=list(species=prep$species, ecodiv=prep$ecodiv), n_obs=stan_data$N,
    engine=engine, threads=threads, chains=n_chain, warmup=n_warm, sampling=n_samp,
    grainsize=grainsz, elapsed_min=elapsed), file.path(out_dir, sprintf("%s_meta.rds", model_id)))
  fe_mean <- point_means(fit, scalar_vars)
  fe_tbl <- tibble(variable=names(fe_mean), mean=as.numeric(fe_mean))
  if (engine=="sample") {
    summ2 <- fit$summary(variables=scalar_vars)
    qd <- fit$summary(variables=scalar_vars, ~quantile(.x, c(0.025,0.5,0.975)))
    names(qd) <- c("variable","q025","q500","q975")
    fe_tbl <- summ2 %>% select(variable,mean,sd,rhat,ess_bulk,ess_tail) %>%
      left_join(qd, by="variable") %>% select(variable,mean,sd,q025,q500,q975,rhat,ess_bulk,ess_tail)
  }
  write_csv(fe_tbl, file.path(out_dir, sprintf("%s_fixed_summary.csv", model_id)))
  sp_mean <- point_means_indexed(fit,"b0_sp")
  write_csv(tibble(idx=seq_along(prep$species), SPCD=prep$species, b0_sp=sp_mean[seq_along(prep$species)]),
            file.path(out_dir, sprintf("%s_species_intercepts.csv", model_id)))
  eco_mean <- point_means_indexed(fit,"b0_eco")
  write_csv(tibble(idx=seq_along(prep$ecodiv), ecodiv=prep$ecodiv, b0_eco=eco_mean[seq_along(prep$ecodiv)]),
            file.path(out_dir, sprintf("%s_ecodiv_intercepts.csv", model_id)))
  pm<-point_means(fit,scalar_vars); spm<-point_means_indexed(fit,"b0_sp"); ecm<-point_means_indexed(fit,"b0_eco")
  needcols <- c("DBH1_cm","BAL1_metric","BAL2_metric","BA1_metric","BA2_metric",
               "QMD1_metric","CR1","CR2","YEARS","ln_site_prod","dg_obs",
               "species_idx","ecodiv_idx")
  holdout_df <- as.data.frame(holdout_dat)
  hd <- holdout_df[stats::complete.cases(holdout_df[, needcols]), ]
  message("Held-out rows after complete-case filter: ", nrow(hd), " of ", nrow(holdout_dat))
  b0_total <- pm["mu_b0"] + spm[hd$species_idx] + ecm[hd$ecodiv_idx]
  params <- list(b0_total=b0_total, b1=pm["b1"],b2=pm["b2"],b3=pm["b3"],b4=pm["b4"],b5=pm["b5"],
                 b6=pm["b6"],b7=pm["b7"],b8=pm["b8"],b9=pm["b9"],b10=pm["b10"],K1=pm["K1"],K2=pm["K2"])
  dg_pred <- predict_dg_annualized(params, hd)
  ann_inc_in <- (dg_pred/hd$YEARS)/2.54; resid <- hd$dg_obs - dg_pred
  ok <- is.finite(ann_inc_in) & is.finite(resid)
  ai <- ann_inc_in[ok]; rr <- resid[ok]; oo <- hd$dg_obs[ok]
  stats <- tibble(n=length(rr), n_dropped=sum(!ok),
    bias_cm=mean(rr), rmse_cm=sqrt(mean(rr^2)),
    pseudo_r2=1-sum(rr^2)/sum((oo-mean(oo))^2),
    ann_inc_in_min=min(ai), ann_inc_in_p50=median(ai), ann_inc_in_mean=mean(ai),
    ann_inc_in_p95=as.numeric(quantile(ai,0.95)), ann_inc_in_max=max(ai),
    pct_negative=100*mean(ai<0), pct_over_0p5=100*mean(ai>0.5),
    pct_explosive=100*mean(ai>2.0))
  message("=== HELD-OUT SANITY (annual DBH increment, in/yr) ==="); print(stats)
  write_csv(stats, file.path(out_dir, sprintf("%s_heldout_sanity.csv", model_id)))
  message("METHOD: ", method_used); message("=== Done -> ", out_dir, " ===")
}
