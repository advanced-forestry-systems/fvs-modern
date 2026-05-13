## =============================================================================
## 33b_fit_hcb_wrapper.R
##
## Wrapper around the existing HCB script that ensures the production
## fit is saved as an RDS. The original 33_fit_hcb_organon.R reads
## hardcoded data paths and does not expose --save_draws, which is why
## three days of HCB compute produced no usable posterior. This wrapper
## sources the function definitions, runs the fit, and saves the result.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
})

opts <- OptionParser(option_list = list(
  make_option("--data",       type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds"),
  make_option("--site",       type = "character", default = "cspi"),
  make_option("--out",        type = "character",
              default = "calibration/output/conus/hcb"),
  make_option("--n",          type = "integer",   default = 0L),
  make_option("--chains",     type = "integer",   default = 2L),
  make_option("--warmup",     type = "integer",   default = 500L),
  make_option("--sampling",   type = "integer",   default = 1000L),
  make_option("--seed",       type = "integer",   default = 42L),
  make_option("--min_sp_n",   type = "integer",   default = 5000L),
  make_option("--max_treedepth", type = "integer", default = 10L),
  make_option("--adapt_delta",   type = "double",  default = 0.9)
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS HCB resubmit (wrapper, saves draws)\n")
cat("Data:", opts$data, "\n")
cat("Out :", opts$out, "\n")
cat("==============================================================\n")

## Source the existing function definitions (sys.nframe() != 0 prevents
## the auto-execute block from running)
source("calibration/R/33_fit_hcb_organon.R", local = FALSE)

## Override hardcoded constants
OUT_DIR    <- opts$out
STAN_MODEL <- "calibration/stan/hcb_organon.stan"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Load and prepare data
raw_data <- readRDS(opts$data)
setDT(raw_data)
if (opts$n > 0L && opts$n < nrow(raw_data)) {
  set.seed(opts$seed)
  raw_data <- raw_data[sample(.N, opts$n)]
  cat("Subsampled to", format(nrow(raw_data), big.mark = ","), "rows\n")
}

prep <- prepare_hcb_data(raw_data, site_var = opts$site)

## Fit with our parameters
mod <- cmdstan_model(STAN_MODEL)
stan_data <- make_stan_data_hcb(prep)

## Prior-centered inits (HG playbook)
N_sp <- stan_data$N_sp; N_L1 <- stan_data$N_L1
N_L2 <- stan_data$N_L2; N_L3 <- stan_data$N_L3
P_trait <- if (is.null(stan_data$P_trait)) 0L else stan_data$P_trait
if (P_trait == 0L && is.null(stan_data$W)) {
  # Stan model still expects a W matrix even when traits are disabled
  stan_data$W <- matrix(0, N_sp, 1)
  stan_data$P_trait <- 0L
}
init_fn <- function(chain_id) {
  set.seed(opts$seed + chain_id)
  list(
    h0 = rnorm(1, -1.0, 0.3), h1 = rnorm(1, 0.5, 0.1),
    h2 = rnorm(1, -0.3, 0.1), h3 = rnorm(1, 0.5, 0.1),
    h4 = rnorm(1, -0.05, 0.05), h5 = rnorm(1, -0.2, 0.1),
    gamma    = if (P_trait > 0) rnorm(P_trait, 0, 0.05)
             else array(0, dim = c(0)),
    z_sp_raw = rnorm(N_sp, 0, 0.05),
    z_L1_raw = rnorm(N_L1, 0, 0.05),
    z_L2_raw = rnorm(N_L2, 0, 0.05),
    z_L3_raw = rnorm(N_L3, 0, 0.05),
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
cat(sprintf("\nWall time: %.1f min\n",
            as.numeric(t1 - t0, units = "mins")))

## Save the fit + summary
sd_path <- file.path(OUT_DIR, "hcb_organon_cspi_traits1_fit.rds")
fit$save_object(sd_path)
cat("Saved fit:", sd_path, "\n")

## Targeted summary: HCB Stan uses h0..h5 fixed effects + traits + variance
## components. Skip the slow full fit$summary().
fixef_pat <- "^(h[0-9]+|gamma\\[|trait_effect\\[)"
all_vars  <- fit$metadata()$model_params
keep_vars <- all_vars[grepl(fixef_pat, all_vars) |
                      all_vars %in% c("sigma_sp","sigma_L1",
                                      "sigma_L2","sigma_L3","phi")]
summ <- fit$summary(variables = keep_vars)
write_csv(summ, file.path(OUT_DIR, "hcb_organon_cspi_traits1_summary.csv"))

saveRDS(list(
  stan_file = STAN_MODEL,
  prep_meta = list(sp = prep$species, L1 = prep$L1,
                   L2 = prep$L2, L3 = prep$L3,
                   cspi_shift = prep$cspi_shift),
  summary   = summ,
  n_obs     = nrow(prep$data),
  wall_min  = as.numeric(t1 - t0, units = "mins")
), file.path(OUT_DIR, "hcb_organon_cspi_traits1_meta.rds"), compress = "xz")

cat("Done. HCB posterior banked.\n")
