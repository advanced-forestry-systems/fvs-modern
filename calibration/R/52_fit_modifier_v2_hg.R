## =============================================================================
## 52_fit_modifier_v2_hg.R
##
## End-to-end Phase 3 v2 fit for HG: re-extracts residuals using the v2
## (continuous-time) enriched RDS so the new kernel columns are
## available, then fits modifier_v2.stan. Self-contained; does not
## depend on a pre-saved residual bundle.
##
## v2 brings:
##   - acute_kernel_5yr / 10yr (fire, wind, harvest)
##   - gamma_kernel_3yr / 5yr (insect, disease, cutting)
##   - plant_age_taper_30 / 50 (plantation effect attenuates with STDAGE)
##   - mid-interval year imputation when DSTRBYR / TRTYR is missing
##   - stacked disturbance support (event 2 and 3)
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
  library(readr)
})

opts <- OptionParser(option_list = list(
  make_option("--fit",       type = "character",
              default = "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds"),
  make_option("--meta",      type = "character",
              default = "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_meta.rds"),
  make_option("--pairs",     type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"),
  make_option("--stan_file", type = "character",
              default = "calibration/stan/modifier_v2.stan"),
  make_option("--acute_lambda", type = "integer", default = 5L,
              help = "5 or 10: which acute kernel column to use"),
  make_option("--gamma_lambda", type = "integer", default = 3L,
              help = "3 or 5: which gamma kernel column to use"),
  make_option("--plant_taper",  type = "integer", default = 30L,
              help = "30 or 50: which plantation age taper to use"),
  make_option("--chains",    type = "integer", default = 2L),
  make_option("--warmup",    type = "integer", default = 250L),
  make_option("--sampling",  type = "integer", default = 250L),
  make_option("--seed",      type = "integer", default = 42L),
  make_option("--max_treedepth", type = "integer", default = 8L),
  make_option("--subsample", type = "integer", default = 30000L),
  make_option("--out",       type = "character",
              default = "calibration/output/conus/hg/modifier_v2"),
  make_option("--save_draws", action = "store_true", default = FALSE)
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS Phase 3 v2 modifier fit (HG, continuous-time kernels)\n")
cat("  fit  :", opts$fit, "\n")
cat("  pairs:", opts$pairs, "\n")
cat("  acute_lambda:", opts$acute_lambda, "yr\n")
cat("  gamma_lambda:", opts$gamma_lambda, "yr\n")
cat("  plant_taper :", opts$plant_taper, "\n")
cat("==============================================================\n")

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

fit  <- readRDS(opts$fit)
meta <- readRDS(opts$meta)
sp_levels <- meta$prep_meta$sp
L1_levels <- meta$prep_meta$L1
L2_levels <- meta$prep_meta$L2
L3_levels <- meta$prep_meta$L3
cspi_shift <- meta$prep_meta$cspi_shift

pairs <- readRDS(opts$pairs)
setDT(pairs)
cat(sprintf("v2 pairs: %s rows, %d cols\n",
            format(nrow(pairs), big.mark = ","), ncol(pairs)))

## ---- HG plausibility filter (mirror 32_fit_hg_organon.R) ----
pairs[, hg_obs_a := (HT2 - HT1) / YEARS]
d <- pairs[
  !is.na(HT1) & HT1 >= 1.37 &
  !is.na(HT2) & HT2 > HT1 &
  hg_obs_a > 0.001 & hg_obs_a < 2.0 &
  !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
  !is.na(BA1) & BA1 > 0 &
  !is.na(CCFL1) & CCFL1 >= 0 &
  !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  YEARS >= 1 & YEARS <= 20 &
  !is.na(cspi) & cspi > 0
]

## Map to base-fit levels
rare <- setdiff(unique(d$SPCD), sp_levels)
if (length(rare) > 0) d[SPCD %in% rare, SPCD := -1L]
d <- d[SPCD %in% sp_levels &
       EPA_L1_CODE %in% L1_levels &
       EPA_L2_CODE %in% L2_levels &
       EPA_L3_CODE %in% L3_levels]

d[, sp_idx := match(SPCD,        sp_levels)]
d[, L1_idx := match(EPA_L1_CODE, L1_levels)]
d[, L2_idx := match(EPA_L2_CODE, L2_levels)]
d[, L3_idx := match(EPA_L3_CODE, L3_levels)]

## Derived covariates for eta computation
d[, log_hg_obs_a := log(hg_obs_a)]
d[, ln_cr_adj    := log((CR1 + 0.2) / 1.2)]
d[, sqrt_ba      := sqrt(BA1)]
d[, ba_x_rd      := BA1 * rd_add]
d[, ccfl_x_rd    := CCFL1 * rd_add]
d[, sqrt_years   := sqrt(YEARS)]
d[, ln_cspi_shift := log(cspi + cspi_shift)]
K1h <- meta$K1h; K2h <- meta$K2h; K4h <- meta$K4h
d[, ln_ht_k1 := log(HT1 + K1h)]
d[, ht_k2    := HT1 ^ K2h]
d[, comp     := CCFL1 / log(HT1 + K4h)]

## Pull posterior-mean coefficients
coefs <- fit$summary(variables = NULL, mean = mean)
get_mean <- function(v) {
  hit <- coefs[coefs$variable == v, ]; if (nrow(hit) == 0) NA_real_ else hit$mean
}
vec_mean <- function(base) {
  rx <- paste0("^", base, "\\[[0-9]+\\]$")
  hit <- coefs[grepl(rx, coefs$variable), ]
  if (nrow(hit) == 0) return(numeric(0))
  idx <- as.integer(regmatches(hit$variable,
                               regexpr("(?<=\\[)[0-9]+", hit$variable, perl = TRUE)))
  out <- numeric(max(idx)); out[idx] <- hit$mean; out
}
a <- sapply(0:8, function(i) get_mean(paste0("a", i)))
trait_effect <- vec_mean("trait_effect")
z_sp <- vec_mean("z_sp")
z_L1 <- vec_mean("z_L1"); z_L2 <- vec_mean("z_L2"); z_L3 <- vec_mean("z_L3")
if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))

eta <- a[1] +
  trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
  z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
  a[2] * d$ln_ht_k1 + a[3] * d$ht_k2 + a[4] * d$ln_cr_adj +
  a[5] * d$ln_cspi_shift + a[6] * d$comp + a[7] * d$sqrt_ba +
  a[8] * d$ba_x_rd + a[9] * d$ccfl_x_rd

residual <- d$log_hg_obs_a - eta
cat(sprintf("Residuals: n=%s mean=%.4f sd=%.4f\n",
            format(sum(is.finite(residual)), big.mark = ","),
            mean(residual, na.rm = TRUE),
            sd(residual, na.rm = TRUE)))

d[, residual := residual]
d <- d[is.finite(residual)]

if (opts$subsample > 0L && opts$subsample < nrow(d)) {
  set.seed(opts$seed)
  d <- d[sample(.N, opts$subsample)]
  cat(sprintf("Subsampled to %s rows\n",
              format(nrow(d), big.mark = ",")))
}

## ---- Build stan_data with v2 kernel columns ----
acute_col <- sprintf("acute_kernel_%dyr", opts$acute_lambda)
gamma_col <- sprintf("gamma_kernel_%dyr", opts$gamma_lambda)
trt_acute <- sprintf("trt_acute_kernel_%dyr", opts$acute_lambda)
trt_gamma <- sprintf("trt_gamma_kernel_%dyr", opts$gamma_lambda)
taper_col <- sprintf("plant_age_taper_%d", opts$plant_taper)

# Coerce missing kernel values to 0 (no event = no contribution)
for (c in c(acute_col, gamma_col, trt_acute, trt_gamma)) {
  if (!c %in% names(d)) stop("Missing column ", c, " in v2 pairs")
  v <- d[[c]]; v[is.na(v)] <- 0; d[[c]] <- v
}
# plant_age_taper is 1 (full effect) when STDAGE missing
v <- d[[taper_col]]; v[is.na(v)] <- 1; d[[taper_col]] <- v

# Indicator columns
for (c in c("is_plantation","had_fire_t1","had_wind_t1","had_harvest_t1",
            "had_insect_t1","had_disease_t1","had_cutting_t1","had_site_prep_t1")) {
  v <- d[[c]]; v[is.na(v)] <- 0L; d[[c]] <- as.integer(v)
}

stan_data <- list(
  N_obs              = nrow(d),
  N_L1               = length(L1_levels),
  residual           = d$residual,
  weight             = d$sqrt_years,
  is_plantation      = as.numeric(d$is_plantation),
  plant_age_taper    = d[[taper_col]],
  d_fire             = as.numeric(d$had_fire_t1),
  d_wind             = as.numeric(d$had_wind_t1),
  d_harvest          = as.numeric(d$had_harvest_t1),
  d_insect           = as.numeric(d$had_insect_t1),
  d_disease          = as.numeric(d$had_disease_t1),
  acute_kernel_dstrb = d[[acute_col]],
  gamma_kernel_dstrb = d[[gamma_col]],
  t_cutting          = as.numeric(d$had_cutting_t1),
  t_site_prep        = as.numeric(d$had_site_prep_t1),
  acute_kernel_trt   = d[[trt_acute]],
  gamma_kernel_trt   = d[[trt_gamma]],
  L1_idx             = as.integer(d$L1_idx)
)

cat(sprintf("\nStan data: N_obs=%s  N_L1=%d\n",
            format(stan_data$N_obs, big.mark = ","), stan_data$N_L1))
cat(sprintf("plantation prev=%.3f  fire/insect/disease/wind/harvest=%.4f %.4f %.4f %.4f %.4f\n",
            mean(stan_data$is_plantation),
            mean(stan_data$d_fire), mean(stan_data$d_insect),
            mean(stan_data$d_disease), mean(stan_data$d_wind),
            mean(stan_data$d_harvest)))
cat(sprintf("cutting prev=%.4f  site_prep prev=%.4f\n",
            mean(stan_data$t_cutting), mean(stan_data$t_site_prep)))
cat(sprintf("plant_age_taper mean=%.3f (1.0 = no attenuation)\n",
            mean(stan_data$plant_age_taper)))

mod <- cmdstan_model(opts$stan_file)

init_fn <- function(chain_id) {
  set.seed(opts$seed + chain_id)
  list(
    alpha_0 = 0,
    alpha_plant = 0, alpha_fire = 0, alpha_wind = 0, alpha_harvest = 0,
    alpha_insect = 0, alpha_disease = 0,
    alpha_cutting = 0, alpha_siteprep = 0,
    z_L1_raw = rnorm(stan_data$N_L1, 0, 0.05),
    sigma_L1 = 0.05,
    sigma_resid = max(0.05, sd(stan_data$residual, na.rm = TRUE) * 0.9)
  )
}
inits <- lapply(seq_len(opts$chains), init_fn)

t0 <- Sys.time()
fit_v2 <- mod$sample(
  data            = stan_data,
  chains          = opts$chains,
  parallel_chains = opts$chains,
  iter_warmup     = opts$warmup,
  iter_sampling   = opts$sampling,
  seed            = opts$seed,
  max_treedepth   = opts$max_treedepth,
  init            = inits,
  refresh         = 100
)
t1 <- Sys.time()
cat(sprintf("\nWall time: %.1f min\n",
            as.numeric(t1 - t0, units = "mins")))

key_vars <- c("alpha_0","alpha_plant",
              "alpha_fire","alpha_wind","alpha_harvest",
              "alpha_insect","alpha_disease",
              "alpha_cutting","alpha_siteprep",
              "sigma_L1","sigma_resid")
summ <- fit_v2$summary(variables = key_vars)
print(summ)

tag <- sprintf("hg_modifier_v2_a%d_g%d_t%d",
               opts$acute_lambda, opts$gamma_lambda, opts$plant_taper)
write_csv(summ, file.path(opts$out, sprintf("%s_summary.csv", tag)))
saveRDS(list(
  acute_lambda = opts$acute_lambda,
  gamma_lambda = opts$gamma_lambda,
  plant_taper  = opts$plant_taper,
  N_obs        = stan_data$N_obs,
  summary      = summ,
  wall_min     = as.numeric(t1 - t0, units = "mins")
), file.path(opts$out, sprintf("%s_meta.rds", tag)), compress = "xz")

if (opts$save_draws)
  fit_v2$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))

cat("\nDone. Artifacts in", opts$out, "\n")
