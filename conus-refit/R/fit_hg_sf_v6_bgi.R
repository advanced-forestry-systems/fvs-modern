##=============================================================================
## fit_hg_sf_v6_bgi.R
## Arm 4 species-free HEIGHT GROWTH fit for the v8_forest_eco campaign.
## Stan: stan/hg_organon_speciesfree_v6_bgi.stan (Organon HG, species-free, BGI).
##
## Authored 2026-06-21 by mirroring R/fit_dg_sf_v8_forest_eco.R (arm 4 prep,
## cspiv6 data, output/sf_arm4) and the species-free HG prep of
## R/32_fit_hg_speciesfree_v8.R MINUS the v8 CCH term (the v6_bgi data block
## carries no cch field). Covariates verified against the read v6_bgi data
## block: hg_obs_a, sqrt_years, ln_dbh, ln_ht, ln_cr_adj, bal_log, bgi,
## ba_metric, slope_pct, cos_aspect, indices, W, bgi_knot1, bgi_knot2.
## Response is annual height increment (m/yr); likelihood normal(mu_a, sigma/sqrt_years).
##
## Usage (cwd = ~/fvs-conus):
##   module reset; module load gcc/12.3.0 R/4.4.0; export PYTHONNOUSERSITE=1
##   Rscript R/fit_hg_sf_v6_bgi.R [--smoke] [--subsample=N]
##=============================================================================
suppressPackageStartupMessages({library(data.table); library(cmdstanr)})
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n,d=NULL){m=grep(paste0("^--",n,"="),args,value=TRUE);if(!length(m))return(d);sub(paste0("^--",n,"="),"",m[1])}
hf <- function(n) any(grepl(paste0("^--",n,"$"),args))
STAN_FILE   <- ga("stan_file","stan/hg_organon_speciesfree_v6_bgi.stan")
OUT_DIR     <- ga("outdir","output/sf_arm4/hg")
OUT_NAME    <- ga("outname","hg_sf_v6_bgi")
DATA_FILE   <- ga("data","data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds")
TRAITS_FILE <- ga("traits","traits/species_traits.rds")
SUBSAMPLE   <- as.integer(ga("subsample", NA_character_))
SMOKE       <- hf("smoke")
MIN_OBS_SPECIES <- 5000
trait_cols <- c("wood_specific_gravity","shade_tolerance_num","softwood",
                "leaf_longevity_months","max_ht_m","max_dbh_cm",
                "vulnerability_score","sensitivity")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
cat("== fit_hg_sf_v6_bgi.R ==\nStan:",STAN_FILE,"\nData:",DATA_FILE,"\nSmoke:",SMOKE,"\n\n")
stopifnot(file.exists(STAN_FILE), file.exists(DATA_FILE), file.exists(TRAITS_FILE))

dat <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat("loaded", nrow(dat), "rows\n")
if (!"bgi" %in% names(dat)) stop("DATA lacks 'bgi' column; cspiv6 does not carry BGI. Add BGI to the arm-4 pairs or point --data at the v2 file that has it.")
if (!"HT2" %in% names(dat)) stop("DATA lacks 'HT2'; cannot form height increment.")

## response + predictors (mirrors 32_fit_hg_speciesfree_v8.R, minus cch)
dat[, hg_obs_a := (HT2 - HT1)/YEARS]
dat[, sqrt_years := sqrt(YEARS)]
dat[, ln_dbh := log(DBH1)]
dat[, ln_ht := log(pmax(HT1, 1.5))]
dat[, ln_cr_adj := log((CR1 + 0.2)/1.2)]
dat[, bal_log := log((BAL_SW1 + BAL_HW1) + 5)]
if (!"SLOPE" %in% names(dat)) dat[, SLOPE := 0]
if (!"ASPECT" %in% names(dat)) dat[, ASPECT := 0]
dat[!is.finite(SLOPE), SLOPE := 0]; dat[!is.finite(ASPECT), ASPECT := 0]
dat[, slope_pct := as.numeric(SLOPE)]
dat[, cos_aspect := cos(as.numeric(ASPECT) * pi/180)]
dat[!is.finite(slope_pct), slope_pct := 0]; dat[!is.finite(cos_aspect), cos_aspect := 0]

dat <- dat[
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(HT1) & HT1 > 1.5 & is.finite(HT2) & HT2 > 1.5 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(bgi) & is.finite(BA1) & BA1 >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0 &
  hg_obs_a > -0.5 & hg_obs_a < 5.0
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

traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with=FALSE]
W <- as.matrix(traits_sub[, trait_cols, with=FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[,j]); if (any(na)) W[na,j] <- median(W[!na,j], na.rm=TRUE)
  W[,j] <- (W[,j] - mean(W[,j]))/sd(W[,j])
}

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) { set.seed(42); dat <- dat[sort(sample.int(nrow(dat), SUBSAMPLE))]; cat("subsampled to", nrow(dat), "\n") }

bgi_knots <- quantile(dat$bgi, c(0.33, 0.67), na.rm=TRUE)
cat("BGI knots:", round(bgi_knots,3), "\n")

stan_data <- list(
  N_obs=nrow(dat), N_sp=length(sp_levels), N_L1=length(L1_levels),
  N_L2=length(L2_levels), N_L3=length(L3_levels), N_FT=length(FT_levels),
  P_trait=ncol(W),
  hg_obs_a=dat$hg_obs_a, sqrt_years=dat$sqrt_years,
  ln_dbh=dat$ln_dbh, ln_ht=dat$ln_ht, ln_cr_adj=dat$ln_cr_adj,
  bal_log=dat$bal_log, bgi=dat$bgi, ba_metric=dat$BA1*0.2296,
  slope_pct=dat$slope_pct, cos_aspect=dat$cos_aspect,
  sp_idx=dat$sp_idx, L1_idx=dat$L1_idx, L2_idx=dat$L2_idx, L3_idx=dat$L3_idx,
  FT_idx=dat$FT_idx, W=W,
  bgi_knot1=unname(bgi_knots[1]), bgi_knot2=unname(bgi_knots[2]))

mod <- cmdstan_model(STAN_FILE)
if (SMOKE) { iw<-50; is<-50; ch<-2; cat("SMOKE 50+50 x2\n") } else { iw<-1000; is<-1000; ch<-4; cat("PRODUCTION 1000+1000 x4\n") }
t0 <- Sys.time()
fit <- mod$sample(data=stan_data, chains=ch, parallel_chains=ch,
                  iter_warmup=iw, iter_sampling=is, seed=42,
                  adapt_delta=0.9, max_treedepth=10, refresh=100)
cat("wall min:", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "\n")

fit$save_object(file.path(OUT_DIR, paste0(OUT_NAME,"_fit.rds")))
core <- c("a0","a1","a2","a2_quad","a3","a4","a5","a6","a7","a8","a9a","a9b","a10",
          paste0("gamma[",seq_len(ncol(W)),"]"),
          paste0("gamma_site[",seq_len(ncol(W)),"]"),
          "sigma_L1","sigma_L2","sigma_L3","sigma_FT","sigma_L1_bgi","sigma")
percol <- tryCatch(c(paste0("trait_effect[",seq_along(sp_levels),"]"),
                     paste0("species_site_slope[",seq_along(sp_levels),"]")),
                   error=function(e) character(0))
vars <- c(core, percol)
summ <- tryCatch(fit$summary(variables=vars,"mean","median","sd",~quantile(.x,c(.05,.95)),"rhat","ess_bulk"),
                 error=function(e) fit$summary(variables=core,"mean","median","sd",~quantile(.x,c(.05,.95)),"rhat","ess_bulk"))
fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME,"_summary.csv")))
saveRDS(list(form="hg_organon_speciesfree_v6_bgi", sp=sp_levels, L1=L1_levels,
             L2=L2_levels, L3=L3_levels, FT=FT_levels, trait_cols=trait_cols,
             bgi_knots=unname(bgi_knots), n_obs=nrow(dat)),
        file.path(OUT_DIR, paste0(OUT_NAME,"_meta.rds")))
cat("sigma:\n"); print(fit$summary("sigma","mean","sd","rhat"))
cat("max rhat (a/sigma):", round(max(summ[grepl("^a|sigma",variable), rhat], na.rm=TRUE),3), "\n")
cat("DONE\n")
