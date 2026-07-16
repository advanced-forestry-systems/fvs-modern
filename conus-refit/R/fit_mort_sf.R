##=============================================================================
## fit_mort_sf.R
## Arm 4 species-free MORTALITY fit for the v8_forest_eco campaign.
## Stan: stan/gompit_mortality_speciesfree.stan (gompit/cloglog survival).
##
## Authored 2026-06-21 by mirroring R/fit_dg_sf_v8_forest_eco.R (arm 4 prep,
## cspiv6 data, output/sf_arm4) and the species-free mortality prep of
## R/34_fit_mortality_speciesfree.R, whose data block already matches this Stan
## file exactly: alive, T_years, dbh, dbh_sq, cr_init, ln_csi, bal_metric,
## sqrt_ba_rd, indices, W. Response is per-period survival (annualized in the
## projector via the T_years exposure).
##
## Usage (cwd = ~/fvs-conus):
##   module reset; module load gcc/12.3.0 R/4.4.0; export PYTHONNOUSERSITE=1
##   Rscript R/fit_mort_sf.R [--smoke] [--subsample=N]
##=============================================================================
suppressPackageStartupMessages({library(data.table); library(cmdstanr)})
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){m=grep(paste0("^--",n,"="),args,value=TRUE);if(!length(m))return(d);sub(paste0("^--",n,"="),"",m[1])}
hf <- function(n) any(grepl(paste0("^--",n,"$"),args))
STAN_FILE   <- ga("stan_file","stan/gompit_mortality_speciesfree.stan")
OUT_DIR     <- ga("outdir","output/sf_arm4/mort")
OUT_NAME    <- ga("outname","mort_sf")
DATA_FILE   <- ga("data","data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds")
TRAITS_FILE <- ga("traits","traits/species_traits.rds")
SUBSAMPLE   <- as.integer(ga("subsample", NA_character_))
SMOKE       <- hf("smoke")
MIN_OBS_SPECIES <- 5000
trait_cols <- c("wood_specific_gravity","shade_tolerance_num","softwood",
                "leaf_longevity_months","max_ht_m","max_dbh_cm",
                "vulnerability_score","sensitivity")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
cat("== fit_mort_sf.R ==\nStan:",STAN_FILE,"\nData:",DATA_FILE,"\nSmoke:",SMOKE,"\n\n")
stopifnot(file.exists(STAN_FILE), file.exists(DATA_FILE), file.exists(TRAITS_FILE))

dat <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat("loaded", nrow(dat), "rows\n")

## survival flag + predictors (mirrors 34_fit_mortality_speciesfree.R)
dat[, alive := as.integer(TREESTATUS2 == 1)]
if ("climate_si" %in% names(dat)) {
  med <- median(dat$climate_si, na.rm=TRUE)
  dat[!is.finite(climate_si), climate_si := med]
  dat[, ln_csi := log(pmax(climate_si, 0.1))]
} else dat[, ln_csi := 0]
dat[!is.finite(ln_csi), ln_csi := 0]
dat[, rd_ratio := sdi_additive1 / SDImax_brms]
dat[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]

dat <- dat[
  TREESTATUS1 == 1 &
  !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1,2) &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  is.finite(BA1) & BA1 >= 0 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(rd_ratio) & rd_ratio >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
]
cat("after filters:", nrow(dat), "rows\n")
cat("mortality rate:", round(1 - mean(dat$alive), 4), "\n")

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
  alive=dat$alive, T_years=dat$YEARS,
  dbh=dat$DBH1, dbh_sq=dat$DBH1^2, cr_init=dat$CR1, ln_csi=dat$ln_csi,
  bal_metric=(dat$BAL_SW1 + dat$BAL_HW1), sqrt_ba_rd=dat$sqrt_ba_rd,
  sp_idx=dat$sp_idx, L1_idx=dat$L1_idx, L2_idx=dat$L2_idx, L3_idx=dat$L3_idx,
  FT_idx=dat$FT_idx, W=W)

mod <- cmdstan_model(STAN_FILE)
if (SMOKE) { iw<-50; is<-50; ch<-2; cat("SMOKE 50+50 x2\n") } else { iw<-1000; is<-1000; ch<-4; cat("PRODUCTION 1000+1000 x4\n") }
t0 <- Sys.time()
fit <- mod$sample(data=stan_data, chains=ch, parallel_chains=ch,
                  iter_warmup=iw, iter_sampling=is, seed=42,
                  adapt_delta=0.9, max_treedepth=10, refresh=100)
cat("wall min:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "\n")

fit$save_object(file.path(OUT_DIR, paste0(OUT_NAME,"_fit.rds")))
core <- c("b0",paste0("b",1:6),
          paste0("gamma[",seq_len(ncol(W)),"]"),
          "sigma_L1","sigma_L2","sigma_L3","sigma_FT")
percol <- tryCatch(paste0("trait_effect[",seq_along(sp_levels),"]"), error=function(e) character(0))
vars <- c(core, percol)
summ <- tryCatch(fit$summary(variables=vars,"mean","median","sd",~quantile(.x,c(.05,.95)),"rhat","ess_bulk"),
                 error=function(e) fit$summary(variables=core,"mean","median","sd",~quantile(.x,c(.05,.95)),"rhat","ess_bulk"))
fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME,"_summary.csv")))
saveRDS(list(form="gompit_mortality_speciesfree", sp=sp_levels, L1=L1_levels,
             L2=L2_levels, L3=L3_levels, FT=FT_levels, trait_cols=trait_cols,
             n_obs=nrow(dat)),
        file.path(OUT_DIR, paste0(OUT_NAME,"_meta.rds")))
cat("=== b coefs ===\n"); print(summ[grepl("^b[0-9]",variable)])
cat("max rhat (b/sigma):", round(max(summ[grepl("^b|sigma",variable), rhat], na.rm=TRUE),3), "\n")
cat("DONE\n")
