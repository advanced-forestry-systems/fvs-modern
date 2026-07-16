## =============================================================================
## 36c_preflight_htdbh_wykoff_simple.R
##
## FVS-CONUS: HT-DBH Wykoff preflight with simplified L1 + species hierarchy.
## Sources 36_fit_ht_dbh.R for data prep and trait matrix, but overrides the
## Stan data packer to drop L2/L3 indices.
##
## Usage:
##   Rscript calibration/R/36c_preflight_htdbh_wykoff_simple.R --n 60000
##
## Author: A. Weiskittel, 2026-04-19
## =============================================================================

Sys.setenv(FVS_CONUS_SKIP_DRIVER = "1")

## Locate helpers using same pattern as 36b
find_htdbh_helpers <- function() {
  here <- normalizePath(
    sub("--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
    mustWork = FALSE)
  same_dir <- if (!is.na(here)) file.path(dirname(here),
                                           "36_fit_ht_dbh.R") else NA
  candidates <- c(
    same_dir,
    "scripts/36_fit_ht_dbh.R",
    "calibration/R/36_fit_ht_dbh.R",
    "./36_fit_ht_dbh.R"
  )
  hit <- candidates[!is.na(candidates) & file.exists(candidates)]
  if (!length(hit))
    stop("Cannot locate 36_fit_ht_dbh.R (tried: ",
         paste(candidates, collapse = ", "), ")")
  hit[1]
}
source(find_htdbh_helpers(), chdir = FALSE, echo = FALSE, local = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(optparse)
})

## Simplified data packer: L1 + species only
pack_wykoff_simple <- function(prep, tmat) {
  d <- prep$data
  W <- tmat$W
  if (ncol(W) == 0L) W <- matrix(0, nrow = length(prep$sp), ncol = 1)
  list(
    N_obs         = nrow(d),
    N_sp          = length(prep$sp),
    N_L1          = length(prep$L1),
    P_trait       = tmat$P_trait,
    ht_obs        = d$HT1,
    dbh           = d$DBH1,
    bal           = d$BAL1,
    sqrt_ba       = d$sqrt_ba,
    ln_cspi_shift = d$ln_cspi_shift,
    ba_x_rd       = d$ba_x_rd,
    bal_x_rd      = d$bal_x_rd,
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
              default = "calibration/stan/ht_dbh_wykoff_simple.stan"),
  make_option("--site",       type = "character", default = "cspi"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--traits",     type = "integer",   default = 1L),
  make_option("--out",        type = "character",
              default = "calibration/output/conus/ht_dbh/preflight"),
  make_option("--n",          type = "integer",   default = 60000L),
  make_option("--chains",     type = "integer",   default = 4L),
  make_option("--warmup",     type = "integer",   default = 1000L),
  make_option("--sampling",   type = "integer",   default = 1000L),
  make_option("--seed",       type = "integer",   default = 42L),
  make_option("--min_sp_n",   type = "integer",   default = 2000L)
)) |> parse_args()

stopifnot("Stan file not found" = file.exists(opts$stan_file))
dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

cat("==============================================================\n")
cat("FVS-CONUS HT-DBH Wykoff (simplified L1+sp) preflight\n")
cat("Stan file:    ", opts$stan_file, "\n")
cat("Site variable:", opts$site, "\n")
cat("Chains:       ", opts$chains, "  warmup:", opts$warmup,
    "  sampling:", opts$sampling, "\n")
cat("==============================================================\n\n")

dat_in <- readRDS(opts$data); setDT(dat_in)
set.seed(opts$seed)
if (opts$n > 0 && opts$n < nrow(dat_in)) {
  dat_in <- dat_in[sample(.N, opts$n)]
  message("Preflight subsample: ", format(nrow(dat_in), big.mark = ","), " rows")
}

prep <- prepare_htdbh_data(dat_in, site_var = opts$site,
                           min_sp_n = opts$min_sp_n)
tmat <- build_trait_matrix(prep, opts$traits_rds, enable = (opts$traits == 1L))
sdata <- pack_wykoff_simple(prep, tmat)

message("Stan data: N_obs=", sdata$N_obs, "  N_sp=", sdata$N_sp,
        "  N_L1=", sdata$N_L1, "  P_trait=", sdata$P_trait)

mod <- cmdstan_model(opts$stan_file)

t0 <- Sys.time()
fit <- mod$sample(
  data            = sdata,
  chains          = opts$chains,
  parallel_chains = opts$chains,
  iter_warmup     = opts$warmup,
  iter_sampling   = opts$sampling,
  seed            = opts$seed,
  adapt_delta     = 0.9,
  max_treedepth   = 12,
  refresh         = 100
)
t1 <- Sys.time()
cat("\nWall time:", round(as.numeric(t1 - t0, units = "mins"), 1), "min\n")

fixef_pat <- "^(a_|b0|b1|s0|s1|gamma\\[)"
summ_fx <- fit$summary() |>
  dplyr::filter(grepl(fixef_pat, variable) |
                variable %in% c("sigma_sp", "sigma_L1"))

tag <- sprintf("htdbh_wykoff_simple_%s_traits%d", opts$site, opts$traits)
write_csv(summ_fx, file.path(opts$out, sprintf("%s_summary.csv", tag)))
try(saveRDS(loo::loo(fit$draws("log_lik", format="draws_matrix")), file.path(opts$out, sprintf("%s_loo.rds", tag))))
saveRDS(list(
  form       = "wykoff_simple",
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

cat("\nDone. Outputs under:", opts$out, "\n")
