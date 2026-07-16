##=============================================================================
## fit_dg_sf_v8_forest_eco.R
## Fit the species-free v8 forest_eco diameter-growth model
## (dg_kuehne2022_speciesfree_v8_forest_eco.stan) for arm 4.
##
## Authored 2026-06-20 by combining the species-free prep of
## 32_fit_dg_kuehne_speciesfree.R with the v8 covariate construction of
## 32_fit_dg_kuehne_v8.R. Covariate transforms verified RAW-scale against the
## Stan priors (b7 ~ N(-0.001,0.005) on raw CCFL1; b6 ~ N(0.4,0.3) on ln(SICOND)).
##
## Data block (26 fields) supplied; ln_sicond_sq and sicond_x_rdadd are computed
## in the Stan transformed-data block, species_site_slope and z_* are internal.
##
## Usage on Cardinal (cwd = ~/fvs-conus):
##   module reset; module load gcc/12.3.0 R/4.4.0
##   Rscript fit_dg_sf_v8_forest_eco.R [--smoke] [--subsample=N]
##=============================================================================
suppressPackageStartupMessages({library(data.table); library(cmdstanr)})
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){m=grep(paste0("^--",n,"="),args,value=TRUE);if(!length(m))return(d);sub(paste0("^--",n,"="),"",m[1])}
hf <- function(n) any(grepl(paste0("^--",n,"$"),args))
STAN_FILE <- ga("stan_file","stan/dg_kuehne2022_speciesfree_v8_forest_eco.stan")
OUT_DIR   <- ga("outdir","/fs/scratch/PUOM0008/crsfaaron/sf_arm4/dg")  # scratch (home over quota)
OUT_NAME  <- ga("outname","dg_sf_v8_forest_eco")
DATA_FILE <- ga("data","data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds")
TRAITS_FILE <- ga("traits","traits/species_traits.rds")
SUBSAMPLE <- as.integer(ga("subsample", NA_character_))
SMOKE <- hf("smoke")
MIN_OBS_SPECIES <- 5000
CSPI_SHIFT <- 1.0
trait_cols <- c("wood_specific_gravity","shade_tolerance_num","softwood",
                "leaf_longevity_months","max_ht_m","max_dbh_cm",
                "vulnerability_score","sensitivity")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
cat("== fit_dg_sf_v8_forest_eco.R ==\nStan:",STAN_FILE,"\nData:",DATA_FILE,"\nSmoke:",SMOKE,"\n\n")
stopifnot(file.exists(STAN_FILE), file.exists(DATA_FILE), file.exists(TRAITS_FILE))

dat <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat("loaded", nrow(dat), "rows\n")

## derived response + predictors (RAW scale, verified vs Stan priors)
dat[, dg_obs_a := (DBH2 - DBH1)/YEARS]
dat[, sqrt_years := sqrt(YEARS)]
dat[, ln_dbh := log(DBH1)]
dat[, ln_cr_adj := log((CR1 + 0.2)/1.2)]
dat[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
dat[, ln_sicond := log(pmax(SICOND, 1))]
dat[, ln_elev := log(pmax(ELEV, 1))]
dat[, rd_additive := sdi_additive1 / SDImax_brms]
dat[, sdi_complexity := sdi_additive1 / pmax(SDI1, 1.0)]
dat[, is_plant := as.numeric(is_plantation)]

dat <- dat[
  is.finite(dg_obs_a) & dg_obs_a > -0.5 & dg_obs_a < 5.0 &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(SICOND) & SICOND > 0 & is.finite(ELEV) & ELEV > 0 &
  is.finite(CCFL1) & CCFL1 >= 0 & is.finite(is_plant) &
  is.finite(rd_additive) & rd_additive > 0 & rd_additive < 3.0 &
  is.finite(sdi_complexity) & sdi_complexity > 0 & sdi_complexity < 10 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" & !is.na(EPA_L2_CODE) & EPA_L2_CODE != "" &
  !is.na(EPA_L3_CODE) & EPA_L3_CODE != "" & !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
]
cat("after filters:", nrow(dat), "rows\n")

sp_counts <- dat[, .N, by=SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("after species filter (n>=",MIN_OBS_SPECIES,"):", nrow(dat), "rows;", nrow(sp_counts), "species\n")

sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))
FT_levels <- sort(unique(as.integer(dat$FORTYPCD_cond1)))
dat[, sp_idx := match(SPCD, sp_levels)]
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
dat[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]
cat("N_sp/L1/L2/L3/FT =", length(sp_levels), length(L1_levels), length(L2_levels),
    length(L3_levels), length(FT_levels), "\n")

## trait matrix W (standardized), aligned to sp_levels
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with=FALSE]
W <- as.matrix(traits_sub[, trait_cols, with=FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[,j]); if (any(na)) W[na,j] <- median(W[!na,j], na.rm=TRUE)
  W[,j] <- (W[,j] - mean(W[,j]))/sd(W[,j])
}

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) { set.seed(42); dat <- dat[sort(sample.int(nrow(dat), SUBSAMPLE))]; cat("subsampled to", nrow(dat), "\n") }

stan_data <- list(
  N_obs=nrow(dat), N_sp=length(sp_levels), N_L1=length(L1_levels),
  N_L2=length(L2_levels), N_L3=length(L3_levels), N_FT=length(FT_levels),
  P_trait=ncol(W),
  dg_obs_a=dat$dg_obs_a, sqrt_years=dat$sqrt_years,
  ln_dbh=dat$ln_dbh, dbh=dat$DBH1, ln_cr_adj=dat$ln_cr_adj,
  ln_bal_sw_adj=dat$ln_bal_sw_adj, bal_hw=dat$BAL_HW1,
  ln_sicond=dat$ln_sicond, ccfl1=dat$CCFL1, is_plantation=dat$is_plant,
  ln_elev=dat$ln_elev, sdi_complexity=dat$sdi_complexity, rd_additive=dat$rd_additive,
  sp_idx=dat$sp_idx, L1_idx=dat$L1_idx, L2_idx=dat$L2_idx, L3_idx=dat$L3_idx,
  FT_idx=dat$FT_idx, W=W)

mod <- cmdstan_model(STAN_FILE)
WARM <- as.integer(ga("warmup", NA_character_)); SAMP <- as.integer(ga("sampling", NA_character_))
if (SMOKE) { iw<-50; is<-50; ch<-2; cat("SMOKE 50+50 x2\n")
} else { iw <- if(!is.na(WARM)) WARM else 1000; is <- if(!is.na(SAMP)) SAMP else 1000; ch<-4
         cat(sprintf("PRODUCTION %d+%d x4\n", iw, is)) }
t0 <- Sys.time()
fit <- mod$sample(data=stan_data, chains=ch, parallel_chains=ch,
                  iter_warmup=iw, iter_sampling=is, seed=42,
                  adapt_delta=0.9, max_treedepth=10, refresh=100)
cat("wall min:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "\n")

fit$save_object(file.path(OUT_DIR, paste0(OUT_NAME,"_fit.rds")))
vars <- c(paste0("b",0:12), paste0("gamma[",seq_len(ncol(W)),"]"),
          paste0("gamma_site[",seq_len(ncol(W)),"]"),
          "sigma_L1","sigma_L2","sigma_L3","sigma_FT","sigma_L1_csi","sigma",
          paste0("trait_effect[",seq_along(sp_levels),"]"),
          paste0("species_site_slope[",seq_along(sp_levels),"]"))
summ <- fit$summary(variables=vars, "mean","median","sd",~quantile(.x,c(.05,.95)),"rhat","ess_bulk")
fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME,"_summary.csv")))
saveRDS(list(form="dg_sf_v8_forest_eco", sp=sp_levels, L1=L1_levels, L2=L2_levels,
             L3=L3_levels, FT=FT_levels, trait_cols=trait_cols, W_mean_sd="standardized",
             cspi_shift=CSPI_SHIFT, n_obs=nrow(dat)),
        file.path(OUT_DIR, paste0(OUT_NAME,"_meta.rds")))
cat("sigma:\n"); print(fit$summary("sigma","mean","sd","rhat"))
cat("max rhat (b/sigma):", round(max(summ[grepl("^b|sigma",variable), rhat], na.rm=TRUE),3), "\n")
cat("DONE\n")
