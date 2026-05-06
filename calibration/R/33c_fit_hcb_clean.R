## =============================================================================
## 33c_fit_hcb_clean.R
##
## Clean HCB driver written from scratch in the modern pattern (matches
## 32_fit_hg_organon.R style). Replaces the wrapper approach in 33b
## which depended on functions from the legacy 33_fit_hcb_organon.R
## that don't return standard stan_data fields.
##
## Reads pairs metric RDS, builds stan_data with N_sp/N_L1/N_L2/N_L3/W
## the same way HG does, and fits hcb_organon.stan.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
  library(readr)
  library(dplyr)
})

DEFAULT_TRAIT_COLS <- c("wood_specific_gravity", "shade_tolerance_num",
                        "softwood", "leaf_longevity_months",
                        "max_ht_m", "max_dbh_cm",
                        "vulnerability_score", "sensitivity",
                        "sla_mm2_mg", "seed_mass_mg")

opt_list <- list(
  make_option("--data",       type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds"),
  make_option("--stan_file",  type = "character",
              default = "calibration/stan/hcb_organon.stan"),
  make_option("--site",       type = "character", default = "cspi"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--traits",     type = "integer", default = 1L),
  make_option("--out",        type = "character",
              default = "calibration/output/conus/hcb"),
  make_option("--n",          type = "integer", default = 0L),
  make_option("--chains",     type = "integer", default = 2L),
  make_option("--warmup",     type = "integer", default = 500L),
  make_option("--sampling",   type = "integer", default = 1000L),
  make_option("--seed",       type = "integer", default = 42L),
  make_option("--min_sp_n",   type = "integer", default = 5000L),
  make_option("--max_treedepth", type = "integer", default = 10L),
  make_option("--adapt_delta",   type = "double",  default = 0.9),
  make_option("--save_draws", action = "store_true", default = TRUE)
)
opts <- parse_args(OptionParser(option_list = opt_list))

cat("==============================================================\n")
cat("FVS-CONUS HCB clean driver\n")
cat("Data:", opts$data, "\n")
cat("Out :", opts$out, "\n")
cat("==============================================================\n")

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

dat_in <- readRDS(opts$data)
setDT(dat_in)
if (opts$n > 0L && opts$n < nrow(dat_in)) {
  set.seed(opts$seed)
  dat_in <- dat_in[sample(.N, opts$n)]
  cat("Subsampled to", format(nrow(dat_in), big.mark = ","), "rows\n")
}

site_var <- opts$site

## HCB plausibility filter (cross-sectional, derived from CR1)
dat <- dat_in[
  !is.na(HT1) & HT1 >= 1.37 &
  !is.na(CR1) & CR1 > 0 & CR1 < 1.0 &
  !is.na(DBH1) & DBH1 > 0 &
  !is.na(BA1) & BA1 > 0 &
  !is.na(BAL1) & BAL1 >= 0 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  !is.na(get(site_var)) & get(site_var) > 0
]
dat[, hcb1  := (1 - CR1) * HT1]
dat[, ratio := pmin(pmax(hcb1 / HT1, 1e-4), 1 - 1e-4)]

## Rare species pooling
sp_n <- dat[, .N, by = SPCD]
rare <- sp_n[N < opts$min_sp_n, SPCD]
dat[SPCD %in% rare, SPCD := -1L]

cspi_shift <- 1.0 - min(dat[[site_var]], na.rm = TRUE)
cspi_shift <- max(cspi_shift, 0.01)
dat[, ln_cspi_shift := log(get(site_var) + cspi_shift)]

dat[, ln_ht       := log(HT1)]
dat[, ln_dbh      := log(DBH1)]
dat[, bal_over_ht := BAL1 / (HT1 + 1)]
dat[, sqrt_ba     := sqrt(BA1)]

sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(dat$EPA_L1_CODE))
L2_levels <- sort(unique(dat$EPA_L2_CODE))
L3_levels <- sort(unique(dat$EPA_L3_CODE))
dat[, sp_idx := match(SPCD,        sp_levels)]
dat[, L1_idx := match(EPA_L1_CODE, L1_levels)]
dat[, L2_idx := match(EPA_L2_CODE, L2_levels)]
dat[, L3_idx := match(EPA_L3_CODE, L3_levels)]

cat(sprintf("Obs: %s  N_sp: %d  N_L1: %d  N_L2: %d  N_L3: %d\n",
            format(nrow(dat), big.mark = ","),
            length(sp_levels), length(L1_levels),
            length(L2_levels), length(L3_levels)))

## Trait matrix
build_W <- function(sp_levels, traits_path) {
  if (!file.exists(traits_path))
    return(list(W = matrix(0, length(sp_levels), 0), P_trait = 0L))
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
  ratio         = dat$ratio,
  ln_ht         = dat$ln_ht,
  ln_dbh        = dat$ln_dbh,
  bal_over_ht   = dat$bal_over_ht,
  sqrt_ba       = dat$sqrt_ba,
  ln_cspi_shift = dat$ln_cspi_shift,
  sp_idx        = dat$sp_idx,
  L1_idx        = dat$L1_idx,
  L2_idx        = dat$L2_idx,
  L3_idx        = dat$L3_idx,
  W             = W
)

mod <- cmdstan_model(opts$stan_file)

init_fn <- function(chain_id) {
  set.seed(opts$seed + chain_id)
  list(
    h0 = rnorm(1, -1.0, 0.3), h1 = rnorm(1, 0.5, 0.1),
    h2 = rnorm(1, -0.3, 0.1), h3 = rnorm(1, 0.5, 0.1),
    h4 = rnorm(1, -0.05, 0.05), h5 = rnorm(1, -0.2, 0.1),
    gamma    = if (tmat$P_trait > 0) rnorm(tmat$P_trait, 0, 0.05) else array(0, dim = c(0)),
    z_sp_raw = rnorm(stan_data$N_sp, 0, 0.05),
    z_L1_raw = rnorm(stan_data$N_L1, 0, 0.05),
    z_L2_raw = rnorm(stan_data$N_L2, 0, 0.05),
    z_L3_raw = rnorm(stan_data$N_L3, 0, 0.05),
    sigma_sp = 0.1, sigma_L1 = 0.15,
    sigma_L2 = 0.1, sigma_L3 = 0.1,
    phi      = 50
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

if (opts$save_draws) {
  sd_path <- file.path(opts$out, "hcb_organon_cspi_traits1_fit.rds")
  fit$save_object(sd_path)
  cat("Saved fit:", sd_path, "\n")
}

## Targeted summary: HCB Stan uses h0..h5 fixed effects + traits + variance
## components. Skip the full fit$summary() (5+ GB scan, 30+ min) by asking
## cmdstan to summarize only the variables we actually report.
fixef_pat <- "^(h[0-9]+|gamma\\[|trait_effect\\[)"
all_vars  <- fit$metadata()$model_params
keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                      all_vars %in% c("sigma_sp", "sigma_L1",
                                      "sigma_L2", "sigma_L3", "phi")]
summ <- fit$summary(variables = keep_vars)
write_csv(summ, file.path(opts$out, "hcb_organon_cspi_traits1_summary.csv"))
saveRDS(list(
  stan_file = opts$stan_file,
  prep_meta = list(sp = sp_levels, L1 = L1_levels,
                   L2 = L2_levels, L3 = L3_levels,
                   cspi_shift = cspi_shift),
  summary   = summ,
  n_obs     = nrow(dat),
  wall_min  = as.numeric(t1 - t0, units = "mins")
), file.path(opts$out, "hcb_organon_cspi_traits1_meta.rds"), compress = "xz")

cat("Done. HCB clean driver finished.\n")
