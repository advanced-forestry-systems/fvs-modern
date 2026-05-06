## =============================================================================
## 32d_preflight_hg_plantcspi.R
##
## HG preflight that adds plantation main effect AND plantation-by-CSPI
## interaction on top of the centered quadratic CSPI form.
##
## Differs from 32_fit_hg_organon.R in two ways:
##   1. Reads conus_remeasurement_pairs_metric_cond.rds (COND-enriched)
##      instead of conus_remeasurement_pairs_metric.rds, so is_plantation
##      is already attached.
##   2. Packs is_plantation into stan_data and points at the
##      plantcspi.stan file.
##
## Usage:
##   Rscript calibration/R/32d_preflight_hg_plantcspi.R \
##       --data       calibration/data/conus_remeasurement_pairs_metric_cond.rds \
##       --stan_file  calibration/stan/hg_organon_fixedK_plantcspi.stan \
##       --traits_rds calibration/traits/species_traits.rds \
##       --site       cspi --traits 1 \
##       --out        calibration/output/conus/hg/preflight_plantcspi \
##       --n          200000 --chains 2 --warmup 500 --sampling 500 \
##       --seed       42 --min_sp_n 5000 --save_draws \
##       --K1h 1.0 --K2h 0.5 --K4h 2.7
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
  library(dplyr)
  library(readr)
})

DEFAULT_TRAIT_COLS <- c("wood_specific_gravity", "shade_tolerance_num",
                        "softwood", "leaf_longevity_months",
                        "max_ht_m", "max_dbh_cm",
                        "vulnerability_score", "sensitivity",
                        "sla_mm2_mg", "seed_mass_mg")

opt_list <- list(
  make_option("--data",       type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond.rds"),
  make_option("--stan_file",  type = "character",
              default = "calibration/stan/hg_organon_fixedK_plantcspi.stan"),
  make_option("--site",       type = "character", default = "cspi"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--traits",     type = "integer",   default = 1L),
  make_option("--out",        type = "character"),
  make_option("--n",          type = "integer",   default = 0L),
  make_option("--chains",     type = "integer",   default = 2L),
  make_option("--warmup",     type = "integer",   default = 500L),
  make_option("--sampling",   type = "integer",   default = 500L),
  make_option("--seed",       type = "integer",   default = 42L),
  make_option("--min_sp_n",   type = "integer",   default = 5000L),
  make_option("--K1h",        type = "double",    default = 1.0),
  make_option("--K2h",        type = "double",    default = 0.5),
  make_option("--K4h",        type = "double",    default = 2.7),
  make_option("--max_treedepth", type = "integer", default = 10L),
  make_option("--adapt_delta",   type = "double",  default = 0.9),
  make_option("--init_scale",    type = "double",  default = 0.1),
  make_option("--save_draws", action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opts$out))

cat("==============================================================\n")
cat("FVS-CONUS HG plantation x CSPI interaction preflight\n")
cat("Stan file:", opts$stan_file, "\n")
cat("Data    :", opts$data, "\n")
cat("Output  :", opts$out, "\n")
cat("Subsample:", opts$n, " chains:", opts$chains, "\n")
cat("==============================================================\n")

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

dat_in <- readRDS(opts$data)
setDT(dat_in)

if (opts$n > 0L && opts$n < nrow(dat_in)) {
  set.seed(opts$seed)
  dat_in <- dat_in[sample(.N, opts$n)]
  message("Sub-sampled to ", format(nrow(dat_in), big.mark = ","), " rows")
}

site_var <- opts$site
dat_in[, hg_obs_a := (HT2 - HT1) / YEARS]

dat <- dat_in[
  !is.na(HT1) & HT1 >= 1.37 &
  !is.na(HT2) & HT2 > HT1 &
  hg_obs_a > 0.001 & hg_obs_a < 2.0 &
  !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
  !is.na(BA1) & BA1 > 0 &
  !is.na(CCFL1) & CCFL1 >= 0 &
  !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  YEARS >= 1 & YEARS <= 20 &
  !is.na(get(site_var)) & get(site_var) > 0
]

# Coerce missing is_plantation to 0 (treat as natural by default; the
# COND join was 100% match so this should affect very few rows)
dat[, is_plantation := fifelse(is.na(is_plantation), 0L,
                                as.integer(is_plantation))]

sp_n <- dat[, .N, by = SPCD]
rare <- sp_n[N < opts$min_sp_n, SPCD]
dat[SPCD %in% rare, SPCD := -1L]

cspi_shift <- 1.0 - min(dat[[site_var]], na.rm = TRUE)
cspi_shift <- max(cspi_shift, 0.01)
dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]

dat[, log_hg_obs_a := log(hg_obs_a)]
dat[, ln_cr_adj    := log((CR1 + 0.2) / 1.2)]
dat[, sqrt_ba      := sqrt(BA1)]
dat[, ba_x_rd      := BA1   * rd_add]
dat[, ccfl_x_rd    := CCFL1 * rd_add]
dat[, sqrt_years   := sqrt(YEARS)]

sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(dat$EPA_L1_CODE))
L2_levels <- sort(unique(dat$EPA_L2_CODE))
L3_levels <- sort(unique(dat$EPA_L3_CODE))
dat[, sp_idx := match(SPCD,        sp_levels)]
dat[, L1_idx := match(EPA_L1_CODE, L1_levels)]
dat[, L2_idx := match(EPA_L2_CODE, L2_levels)]
dat[, L3_idx := match(EPA_L3_CODE, L3_levels)]

cat(sprintf("Obs after filter: %s\n", format(nrow(dat), big.mark = ",")))
cat(sprintf("Plantation share: %.3f\n", mean(dat$is_plantation)))
cat(sprintf("Species: %d  L1: %d  L2: %d  L3: %d\n",
            length(sp_levels), length(L1_levels),
            length(L2_levels), length(L3_levels)))

## ---- Trait matrix ----
build_W <- function(sp_levels, traits_path) {
  if (!file.exists(traits_path)) return(list(W = matrix(0, length(sp_levels), 0), P_trait = 0L))
  tr <- readRDS(traits_path); setDT(tr)
  keep <- intersect(DEFAULT_TRAIT_COLS, names(tr))
  tr_sp <- tr[match(sp_levels, SPCD), c("SPCD", keep), with = FALSE]
  tr_sp[, SPCD := sp_levels]
  informative <- character(0)
  for (cc in keep) {
    v <- tr_sp[[cc]]
    med <- median(v, na.rm = TRUE); if (is.na(med)) med <- 0
    v[is.na(v)] <- med
    s <- sd(v)
    if (!is.finite(s) || s < 1e-12) {
      cat("  dropping zero-variance trait:", cc, "\n")
      tr_sp[, (cc) := NULL]; next
    }
    tr_sp[, (cc) := (v - mean(v)) / s]
    informative <- c(informative, cc)
  }
  W <- as.matrix(tr_sp[, informative, with = FALSE])
  cat("Traits kept:", paste(informative, collapse = ", "), "\n")
  list(W = W, P_trait = ncol(W))
}
tmat <- if (opts$traits == 1L) build_W(sp_levels, opts$traits_rds) else
        list(W = matrix(0, length(sp_levels), 0), P_trait = 0L)
W <- tmat$W
if (ncol(W) == 0L) W <- matrix(0, length(sp_levels), 1)

stan_data <- list(
  N_obs         = nrow(dat),
  N_sp          = length(sp_levels),
  N_L1          = length(L1_levels),
  N_L2          = length(L2_levels),
  N_L3          = length(L3_levels),
  P_trait       = tmat$P_trait,
  log_hg_obs_a  = dat$log_hg_obs_a,
  sqrt_years    = dat$sqrt_years,
  ht            = dat$HT1,
  ln_cr_adj     = dat$ln_cr_adj,
  ln_cspi_shift = dat$ln_cspi_shift,
  ccfl          = dat$CCFL1,
  sqrt_ba       = dat$sqrt_ba,
  ba_x_rd       = dat$ba_x_rd,
  ccfl_x_rd     = dat$ccfl_x_rd,
  is_plantation = as.numeric(dat$is_plantation),
  K1h           = opts$K1h,
  K2h           = opts$K2h,
  K4h           = opts$K4h,
  sp_idx        = dat$sp_idx,
  L1_idx        = dat$L1_idx,
  L2_idx        = dat$L2_idx,
  L3_idx        = dat$L3_idx,
  W             = W
)

mod <- cmdstan_model(opts$stan_file)

# Prior-centered init list per chain
init_fn <- function(chain_id) {
  set.seed(opts$seed + chain_id)
  list(
    a0           = rnorm(1, -2.5, 0.3),
    a1           = rnorm(1,  0.5, 0.1),
    a2           = rnorm(1, -0.03, 0.01),
    a3           = rnorm(1,  0.8, 0.1),
    a4           = rnorm(1,  0.4, 0.1),
    a4b          = rnorm(1,  0.0, 0.05),
    a_plant      = rnorm(1,  0.0, 0.05),
    a_plant_cspi = rnorm(1,  0.0, 0.05),
    a5           = rnorm(1, -0.005, 0.002),
    a6           = rnorm(1, -0.03, 0.01),
    a7           = rnorm(1,  0.0, 0.005),
    a8           = rnorm(1,  0.0, 0.005),
    gamma        = rnorm(stan_data$P_trait, 0, 0.05),
    z_sp_raw     = rnorm(stan_data$N_sp, 0, 0.05),
    z_L1_raw     = rnorm(stan_data$N_L1, 0, 0.05),
    z_L2_raw     = rnorm(stan_data$N_L2, 0, 0.05),
    z_L3_raw     = rnorm(stan_data$N_L3, 0, 0.05),
    sigma_sp = 0.1, sigma_L1 = 0.15,
    sigma_L2 = 0.1, sigma_L3 = 0.1,
    sigma    = 0.3
  )
}
inits <- lapply(seq_len(opts$chains), init_fn)

t0 <- Sys.time()
fit <- mod$sample(
  data            = stan_data,
  chains          = opts$chains,
  parallel_chains = opts$chains,
  iter_warmup     = opts$warmup,
  iter_sampling   = opts$sampling,
  seed            = opts$seed,
  adapt_delta     = opts$adapt_delta,
  max_treedepth   = opts$max_treedepth,
  init            = inits,
  refresh         = 100
)
t1 <- Sys.time()
cat(sprintf("\nWall time: %.1f min\n", as.numeric(t1 - t0, units = "mins")))

key_vars <- c("a0","a1","a2","a3","a4","a4b","a_plant",
              "a5","a6","a7","a8",
              "sigma_sp","sigma_L1","sigma_L2","sigma_L3","sigma")
summ <- fit$summary(variables = key_vars)
print(summ)
write_csv(summ, file.path(opts$out, "hg_plantcspi_summary.csv"))

if (opts$save_draws) {
  fit$save_object(file.path(opts$out, "hg_plantcspi_fit.rds"))
}

saveRDS(list(
  stan_file     = opts$stan_file,
  K1h = opts$K1h, K2h = opts$K2h, K4h = opts$K4h,
  prep_meta     = list(sp = sp_levels, L1 = L1_levels,
                       L2 = L2_levels, L3 = L3_levels,
                       cspi_shift = cspi_shift),
  summary       = summ,
  n_obs         = nrow(dat),
  plantation_share = mean(dat$is_plantation),
  wall_min      = as.numeric(t1 - t0, units = "mins")
), file.path(opts$out, "hg_plantcspi_meta.rds"), compress = "xz")

cat("\nDone. Artifacts in", opts$out, "\n")
