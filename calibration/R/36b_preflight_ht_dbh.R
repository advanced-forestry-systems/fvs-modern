## =============================================================================
## 36b_preflight_ht_dbh.R
##
## Small-sample preflight for the three H-D candidate forms (chapman,
## wykoff, schnute). Matches the style of 31b_preflight_dg.R: sources
## 36_fit_ht_dbh.R for its helpers (prepare_htdbh_data, build_trait_matrix,
## stan_data_packer_ht) and runs a short sample with 500 warmup / 500
## sampling on a stratified subsample before committing a full run.
##
## Usage:
##   Rscript scripts/36b_preflight_ht_dbh.R
##   Rscript scripts/36b_preflight_ht_dbh.R --n 100000 --traits 0
##
## Author: A. Weiskittel, 2026-04-15
## =============================================================================

Sys.setenv(FVS_CONUS_SKIP_DRIVER = "1")

## Locate 36_fit_ht_dbh.R using the same pattern as 31b: try the directory
## this script lives in first, then fall through to common layout candidates.
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

opts <- OptionParser(option_list = list(
  make_option("--data",       type = "character",
              default = "calibration/data/conus_remeasurement_pairs.rds"),
  make_option("--stan_dir",   type = "character", default = "scripts/stan"),
  make_option("--site",       type = "character", default = "cspi"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--traits",     type = "integer",   default = 1L),
  make_option("--out",        type = "character",
              default = "calibration/output/conus/ht_dbh/preflight"),
  make_option("--n",          type = "integer",   default = 60000L),
  make_option("--chains",     type = "integer",   default = 2L),
  make_option("--warmup",     type = "integer",   default = 500L),
  make_option("--sampling",   type = "integer",   default = 500L),
  make_option("--seed",       type = "integer",   default = 42L),
  make_option("--min_sp_n",   type = "integer",   default = 2000L),
  make_option("--dbh_ref",    type = "double",    default = 25.0),
  make_option("--form",       type = "character", default = "all",
              help = paste("H-D form to run (or 'all'):",
                           paste(KNOWN_FORMS, collapse = " | ")))
)) |> parse_args()

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

dat_in <- readRDS(opts$data); setDT(dat_in)
set.seed(opts$seed)
if (opts$n > 0 && opts$n < nrow(dat_in)) {
  dat_in <- dat_in[sample(.N, opts$n)]
  message("Preflight subsample: ", format(nrow(dat_in), big.mark=","), " rows")
}

prep <- prepare_htdbh_data(dat_in, site_var = opts$site,
                           min_sp_n = opts$min_sp_n)
tmat <- build_trait_matrix(prep, opts$traits_rds, enable = (opts$traits == 1L))

forms_to_run <- if (opts$form == "all") KNOWN_FORMS else opts$form
stopifnot(all(forms_to_run %in% KNOWN_FORMS))

results <- list()
for (form in forms_to_run) {
  cat("\n==============================================================\n")
  cat("Preflight:", form, "\n")
  cat("==============================================================\n")

  stan_file <- file.path(opts$stan_dir,
                         switch(form,
                                chapman = "ht_dbh_chapman.stan",
                                wykoff  = "ht_dbh_wykoff.stan",
                                schnute = "ht_dbh_schnute.stan"))
  packer <- stan_data_packer_ht(form)
  sdata  <- if (form == "schnute") packer(prep, tmat, dbh_ref = opts$dbh_ref)
            else                    packer(prep, tmat)

  ad <- if (form == "schnute") 0.95 else 0.9
  td <- if (form == "schnute") 13 else 12

  mod <- cmdstan_model(stan_file)
  t0  <- Sys.time()
  fit <- mod$sample(
    data            = sdata,
    chains          = opts$chains,
    parallel_chains = opts$chains,
    iter_warmup     = opts$warmup,
    iter_sampling   = opts$sampling,
    seed            = opts$seed,
    adapt_delta     = ad,
    max_treedepth   = td,
    refresh         = 100
  )
  t1  <- Sys.time()

  diag <- tryCatch(fit$diagnostic_summary(), error = function(e) list())
  summ <- fit$summary() |>
    dplyr::filter(grepl("^(a0|a_|b0|b1|b_rate|c_shape|a_rate|b_shape|s0|s1|sigma_)",
                        variable)) |>
    dplyr::select(variable, mean, median, sd, q5, q95, rhat, ess_bulk)

  write_csv(summ, file.path(opts$out,
                             sprintf("preflight_htdbh_%s_traits%d.csv",
                                     form, opts$traits)))
  results[[form]] <- list(
    form      = form,
    wall_min  = as.numeric(t1 - t0, units = "mins"),
    divergent = sum(diag$num_divergent %||% 0),
    max_td    = max(diag$num_max_treedepth %||% 0),
    ebfmi_min = min(diag$ebfmi %||% NA_real_),
    rhat_max  = max(summ$rhat, na.rm = TRUE),
    ess_min   = min(summ$ess_bulk, na.rm = TRUE)
  )
  print(results[[form]])
}

saveRDS(results, file.path(opts$out,
                           sprintf("preflight_htdbh_traits%d_summary.rds",
                                   opts$traits)))

cat("\n==============================================================\n")
cat("HT-DBH preflight done.  Results saved to ", opts$out, "\n")
cat("==============================================================\n")
