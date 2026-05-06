## =============================================================================
## 34b_fit_mortality_logit.R
##
## FVS-CONUS: Individual tree mortality driver for mort_logit_simple.stan
## (logit link with simplified L1 + species random intercept hierarchy).
##
## This is a lightweight adaptation of 34_fit_mortality.R that drops L2/L3
## random effects from the Stan data list to match the simplified model.
## All data prep (filtering, covariate construction, trait matrix) is
## shared via sourcing 34_fit_mortality.R.
##
## Usage:
##   Rscript scripts/34b_fit_mortality_logit.R --site cspi
##   Rscript scripts/34b_fit_mortality_logit.R --site cspi --n 50000
##
## Author: A. Weiskittel, 2026-04-18
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(optparse)
})

## Source the shared functions from 34 (data prep, trait builder)
## Suppress auto-execution of the driver block by setting an env var
Sys.setenv(FVS_CONUS_SKIP_DRIVER = "1")
## Resolve path relative to the project root (working directory on Cardinal
## is ~/fvs-modern, so the file lives at calibration/R/34_fit_mortality.R).
## Also check the local scripts/ path as a fallback for development.
mort_source <- if (file.exists("calibration/R/34_fit_mortality.R")) {
  "calibration/R/34_fit_mortality.R"
} else if (file.exists("scripts/34_fit_mortality.R")) {
  "scripts/34_fit_mortality.R"
} else {
  stop("Cannot find 34_fit_mortality.R in calibration/R/ or scripts/")
}
source(mort_source, local = TRUE)
Sys.setenv(FVS_CONUS_SKIP_DRIVER = "")

## ---- Stan data packer (simplified: L1 + species only) --------------------
make_stan_data_mort_simple <- function(prep, tmat) {
  d <- prep$data
  W <- tmat$W
  if (ncol(W) == 0L) W <- matrix(0, nrow = length(prep$sp), ncol = 1)
  list(
    N_obs         = nrow(d),
    N_sp          = length(prep$sp),
    N_L1          = length(prep$L1),
    P_trait       = tmat$P_trait,
    alive         = d$alive,
    log_years     = d$log_years,
    dbh           = d$DBH1,
    dbh2          = d$dbh2,
    bal_over_ba   = d$bal_over_ba,
    cr            = d$CR1,
    sqrt_ba_rd    = d$sqrt_ba_rd,
    ln_cspi_shift = d$ln_cspi_shift,
    sp_idx        = d$sp_idx,
    L1_idx        = d$L1_idx,
    W             = W
  )
}

## ---- Driver ---------------------------------------------------------------
opts <- OptionParser(option_list = list(
  make_option("--data",       type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds"),
  make_option("--stan_file",  type = "character",
              default = "calibration/stan/mort_logit_simple.stan"),
  make_option("--site",       type = "character", default = "cspi"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--traits",     type = "integer",   default = 1L),
  make_option("--out",        type = "character",
              default = "calibration/output/conus/mortality"),
  make_option("--n",          type = "integer",   default = 0L),
  make_option("--chains",     type = "integer",   default = 4L),
  make_option("--warmup",     type = "integer",   default = 1000L),
  make_option("--sampling",   type = "integer",   default = 1000L),
  make_option("--seed",       type = "integer",   default = 42L),
  make_option("--min_sp_n",   type = "integer",   default = 5000L),
  make_option("--save_draws", action = "store_true", default = FALSE)
)) |> parse_args()

stopifnot("Stan file not found" = file.exists(opts$stan_file))
dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("FVS-CONUS Mortality (logit, simple L1+sp) fit\n")
cat("Stan file:    ", opts$stan_file, "\n")
cat("Site variable:", opts$site, "\n")
cat("Trait prior:  ", if (opts$traits == 1) "enabled" else "disabled", "\n")
cat("Data:         ", opts$data, "\n")
cat("Chains:       ", opts$chains, "  warmup:", opts$warmup,
    "  sampling:", opts$sampling, "\n")
cat("==============================================================\n\n")

dat_in <- readRDS(opts$data)
setDT(dat_in)
if (opts$n > 0 && opts$n < nrow(dat_in)) {
  set.seed(opts$seed)
  dat_in <- dat_in[sample(.N, opts$n)]
  message("Sub-sampled to ", format(nrow(dat_in), big.mark = ","), " rows")
}

prep  <- prepare_mort_data(dat_in, site_var = opts$site,
                           min_sp_n = opts$min_sp_n)
tmat  <- build_trait_matrix(prep, opts$traits_rds,
                            enable = (opts$traits == 1L))
sdata <- make_stan_data_mort_simple(prep, tmat)

message("Stan data: N_obs=", sdata$N_obs, "  N_sp=", sdata$N_sp,
        "  N_L1=", sdata$N_L1, "  P_trait=", sdata$P_trait)

mod <- cmdstan_model(opts$stan_file)

t0  <- Sys.time()
fit <- mod$sample(
  data            = sdata,
  chains          = opts$chains,
  parallel_chains = opts$chains,
  iter_warmup     = opts$warmup,
  iter_sampling   = opts$sampling,
  seed            = opts$seed,
  adapt_delta     = 0.90,         # logit is better behaved than cloglog
  max_treedepth   = 10,
  refresh         = 200
)
t1 <- Sys.time()
cat("\nWall time:", round(as.numeric(t1 - t0, units = "mins"), 1), "min\n")

fixef_pat <- "^(m[0-9]+|gamma\\[|trait_effect\\[)"
## Targeted summary: avoid the 30+ min stall on full fit$summary().
all_vars  <- fit$metadata()$model_params
keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                      all_vars %in% c("sigma_sp", "sigma_L1")]
summ_fx   <- fit$summary(variables = keep_vars)

tag <- sprintf("mort_logit_simple_%s_traits%d", opts$site, opts$traits)
write_csv(summ_fx, file.path(opts$out, sprintf("%s_summary.csv", tag)))
saveRDS(list(
  site_var   = opts$site,
  traits     = opts$traits == 1L,
  trait_cols = tmat$trait_cols,
  stan_file  = opts$stan_file,
  prep_meta  = list(sp = prep$sp, L1 = prep$L1,
                    cspi_shift = prep$cspi_shift),
  summary    = summ_fx,
  n_obs      = nrow(prep$data),
  wall_min   = as.numeric(t1 - t0, units = "mins")
), file.path(opts$out, sprintf("%s_meta.rds", tag)))

if (opts$save_draws) {
  fit$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))
}

cat("\nDone. Outputs under:", opts$out, "\n")
