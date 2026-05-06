## =============================================================================
## 31b_preflight_dg.R
##
## FVS-CONUS: Short pre-flight fit for DG model candidates. Validates that
## each Stan file compiles, samples cleanly, and produces sane posterior
## summaries plus a held-out plot MAB on a small subsample, before
## committing to the multi day production run of 31_fit_dg_multimodel.R.
##
## Capabilities:
##   * Run a single form:          --form weiskittel
##   * Run all four in sequence:   --form all
##   * Plot-level held-out split for out-of-sample MAB / RMSE
##   * Convergence diagnostics (Rhat, ESS bulk / tail)
##   * Per-form timing and summary table
##
## Typical wall time on Cardinal cpu partition (n=25,000, 200+200, 1 chain):
##   organon     ~20 min
##   weiskittel  ~15 min
##   kuehne      ~15 min
##   compound    ~25 min
##
## Usage:
##   Rscript scripts/31b_preflight_dg.R                         # all 4 forms
##   Rscript scripts/31b_preflight_dg.R --form kuehne --n 25000
##   Rscript scripts/31b_preflight_dg.R --form all --n 50000 --chains 2
##
## Author: A. Weiskittel, 2026-04-15
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(optparse)
})

KNOWN_FORMS <- c("organon", "weiskittel", "kuehne", "compound", "organon_gj")

opts <- OptionParser(option_list = list(
  make_option("--form",     type = "character", default = "all",
              help = paste("Model form (or 'all'):",
                           paste(KNOWN_FORMS, collapse = " | "))),
  make_option("--data",     type = "character",
              default = "calibration/data/conus_remeasurement_pairs.rds"),
  make_option("--stan_dir", type = "character",
              default = "scripts/stan"),
  make_option("--site",     type = "character", default = "cspi"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--traits",   type = "integer",  default = 1L,
              help = "1 = enable trait-informed species prior, 0 = disable"),
  make_option("--out",      type = "character",
              default = "calibration/output/conus/preflight"),
  make_option("--n",        type = "integer",  default = 25000L),
  make_option("--warmup",   type = "integer",  default = 200L),
  make_option("--sampling", type = "integer",  default = 200L),
  make_option("--chains",   type = "integer",  default = 1L),
  make_option("--seed",     type = "integer",  default = 42L),
  make_option("--holdout_plots", type = "double", default = 0.15,
              help = "Fraction of plots to hold out for validation MAB"),
  make_option("--min_sp_n", type = "integer",  default = 2000L)
)) |> parse_args()

forms_to_run <- if (opts$form == "all") KNOWN_FORMS else opts$form
stopifnot(all(forms_to_run %in% KNOWN_FORMS))

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

## Source the multi-model driver for prepare_dg_data + stan_data_packer.
## Set the sentinel env var so 31 stops before its own driver block runs.
Sys.setenv(FVS_CONUS_SKIP_DRIVER = "1")

## Locate 31_fit_dg_multimodel.R. Try candidates in order so the same script
## works from both the local workspace (scripts/) and Cardinal (calibration/R/).
find_helpers <- function() {
  here <- normalizePath(
    sub("--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
    mustWork = FALSE)
  same_dir <- if (!is.na(here)) file.path(dirname(here),
                                          "31_fit_dg_multimodel.R") else NA
  candidates <- c(
    same_dir,
    "scripts/31_fit_dg_multimodel.R",
    "calibration/R/31_fit_dg_multimodel.R",
    "./31_fit_dg_multimodel.R"
  )
  hit <- candidates[!is.na(candidates) & file.exists(candidates)]
  if (!length(hit))
    stop("Cannot locate 31_fit_dg_multimodel.R (tried: ",
         paste(candidates, collapse = ", "), ")")
  hit[1]
}

helpers_path <- find_helpers()
message("Sourcing helpers: ", helpers_path)
source(helpers_path, chdir = FALSE, echo = FALSE, local = FALSE)

cat("==============================================================\n")
cat("FVS-CONUS DG Pre-flight (multi-model)\n")
cat("==============================================================\n")
cat("Forms:         ", paste(forms_to_run, collapse = ", "), "\n")
cat("Site variable: ", opts$site, "\n")
cat("Trait prior:   ", if (opts$traits == 1) "enabled" else "disabled", "\n")
cat("Sample N:      ", opts$n, "\n")
cat("Hold-out plots:", opts$holdout_plots, "\n")
cat("Warmup/Samp:   ", opts$warmup, "/", opts$sampling,
    "  Chains:", opts$chains, "\n")
cat("Start:         ", format(Sys.time()), "\n")
cat("==============================================================\n\n")

## ---- Load + stratified subsample ------------------------------------------
cat("Loading:", opts$data, "\n")
raw <- readRDS(opts$data)
setDT(raw)
cat("  Full data:", format(nrow(raw), big.mark = ","), "rows\n")

raw <- raw[!is.na(get(opts$site))]
cat("  With ", opts$site, ":", format(nrow(raw), big.mark = ","), "rows\n")

set.seed(opts$seed)
sub <- raw[sample(.N, min(opts$n, .N))]
cat("  Subsampled:", format(nrow(sub), big.mark = ","), "rows\n\n")

## Plot-level hold-out split (before prepare_dg_data so splits survive
## rare-species pooling inside prepare).
plot_keys <- unique(sub$plot_key)
set.seed(opts$seed + 1)
holdout_keys <- sample(plot_keys,
                       ceiling(length(plot_keys) * opts$holdout_plots))
sub[, in_holdout := plot_key %in% holdout_keys]

cat("Fit plots:     ", length(plot_keys) - length(holdout_keys), "\n")
cat("Hold-out plots:", length(holdout_keys), "\n\n")

## Prepare using the shared function from 31_fit_dg_multimodel.R
prep <- prepare_dg_data(sub[in_holdout == FALSE],
                        site_var = opts$site,
                        min_sp_n = opts$min_sp_n)

tmat <- build_trait_matrix(prep, opts$traits_rds,
                           enable = (opts$traits == 1L))

prep_ho <- prepare_dg_data(sub[in_holdout == TRUE],
                           site_var = opts$site,
                           min_sp_n = opts$min_sp_n)

## Hold-out factor index alignment: species / L1 / L2 / L3 levels in prep_ho
## may not match the levels used to fit prep. Map hold-out codes against
## prep$*. Any new code maps to NA and is dropped from the MAB tally with a
## warning (plot-level split typically keeps these small).
align_idx <- function(codes, levels) match(codes, levels)

prep_ho$data[, sp_idx := align_idx(SPCD,         prep$sp)]
prep_ho$data[, L1_idx := align_idx(EPA_L1_CODE, prep$L1)]
prep_ho$data[, L2_idx := align_idx(EPA_L2_CODE, prep$L2)]
prep_ho$data[, L3_idx := align_idx(EPA_L3_CODE, prep$L3)]

n_ho_drop <- sum(is.na(prep_ho$data$sp_idx) | is.na(prep_ho$data$L1_idx) |
                 is.na(prep_ho$data$L2_idx) | is.na(prep_ho$data$L3_idx))
if (n_ho_drop > 0) {
  cat("  Hold-out rows dropped (unseen factor level):", n_ho_drop, "\n")
}

ho_valid <- prep_ho$data[
  !is.na(sp_idx) & !is.na(L1_idx) & !is.na(L2_idx) & !is.na(L3_idx)
]

## ---- Form-specific R-side eta for hold-out prediction ---------------------
##
## Uses posterior medians of fixed effects and variance components. Random
## effect z_* are pulled by index. This is a quick MAB probe, not a proper
## posterior predictive.
## ---------------------------------------------------------------------------
predict_eta <- function(form, draws, d) {
  m <- function(v) median(draws[[v]])
  zget <- function(prefix, idx) {
    nm <- paste0(prefix, "[", idx, "]")
    vapply(nm, function(x) median(draws[[x]]), numeric(1))
  }
  z_sp <- zget("z_sp", d$sp_idx)
  z_L1 <- zget("z_L1", d$L1_idx)
  z_L2 <- zget("z_L2", d$L2_idx)
  z_L3 <- zget("z_L3", d$L3_idx)
  re   <- z_sp + z_L1 + z_L2 + z_L3

  switch(form,
    organon = {
      m("a0") + re +
        m("a1") * log(d$DBH1 + m("K1")) +
        m("a2") * d$DBH1 ^ m("K2") +
        m("a3") * d$ln_cr_adj +
        m("a4") * d$ln_cspi_shift +
        m("a5") * d$BAL1 / log(d$DBH1 + m("K4")) +
        m("a6") * sqrt(d$BA1) +
        m("a7") * d$ba_x_rd +
        m("a8") * d$bal_x_rd
    },
    weiskittel = {
      m("g0") + re +
        m("g1") * d$ln_dbh + m("g2") * d$DBH1 +
        m("g3") * d$ln_cr_adj +
        m("g4") * d$BAL_SW1 + m("g5") * d$BAL_HW1 +
        m("g6") * d$rd_add +
        m("g7") * d$ln_cspi_shift +
        m("g8") * d$ba_x_rd + m("g9") * d$bal_x_rd
    },
    kuehne = {
      m("b0") + re +
        m("b1") * d$ln_dbh + m("b2") * d$DBH1 +
        m("b3") * d$ln_cr_adj +
        m("b4") * d$ln_bal_sw_adj +
        m("b5") * d$BAL_HW1 +
        m("b6") * d$ln_cspi_shift +
        m("b7") * d$ba_x_rd + m("b8") * d$bal_x_rd
    },
    compound = {
      ## Reparameterized 2026-04-16: separate c_num, c_den on log-transforms
      ln_num <- log(d$DBH1 + 1)
      ln_den <- log(d$CR1 * pmax(d$HT1, 1.37) + 1)
      m("c0") + re +
        m("c_num") * ln_num +
        m("c_den") * ln_den +
        m("c3") * d$BAL_SW1 + m("c4") * d$BAL_HW1 +
        m("c5") * d$rd_add +
        m("c6") * d$ln_cspi_shift +
        m("c7") * d$ba_x_rd + m("c8") * d$bal_x_rd
    },
    organon_gj = {
      m("s0") + re +
        m("s1") * log(d$DBH1 + 1) +
        m("s2") * d$DBH1 +
        m("s3") * d$ln_cr_adj +
        m("s4") * d$BAL1 / log(d$DBH1 + 1) +
        m("s5") * sqrt(d$BA1)
    }
  )
}

## ---- Loop over requested forms --------------------------------------------
results <- list()

for (form in forms_to_run) {
  cat("\n==============================================================\n")
  cat(">>> Pre-flight form:", form, "\n")
  cat("==============================================================\n")

  stan_file <- switch(form,
    organon    = file.path(opts$stan_dir, "dg_organon.stan"),
    weiskittel = file.path(opts$stan_dir, "dg_weiskittel2016.stan"),
    kuehne     = file.path(opts$stan_dir, "dg_kuehne2022.stan"),
    compound   = file.path(opts$stan_dir, "dg_compound_gj.stan"),
    organon_gj = file.path(opts$stan_dir, "dg_organon_gj.stan"))
  stopifnot(file.exists(stan_file))

  ## After 2026-04-16 linearization all forms use the same tuning
  form_adapt_delta   <- 0.9
  form_max_treedepth <- 12

  packer <- stan_data_packer(form)
  sdata  <- packer(prep, tmat)

  mod <- cmdstan_model(stan_file)

  t0  <- Sys.time()
  fit <- mod$sample(
    data            = sdata,
    chains          = opts$chains,
    parallel_chains = opts$chains,
    iter_warmup     = opts$warmup,
    iter_sampling   = opts$sampling,
    seed            = opts$seed,
    adapt_delta     = form_adapt_delta,
    max_treedepth   = form_max_treedepth,
    refresh         = 0,
    show_messages   = FALSE
  )
  t1 <- Sys.time()
  wall_min <- as.numeric(t1 - t0, units = "mins")
  cat(sprintf("  Wall time: %.1f min\n", wall_min))

  ## Convergence summary on structural params
  summ <- fit$summary(variables = c("sigma", "sigma_sp",
                                    "sigma_L1", "sigma_L2", "sigma_L3"))
  summ$rhat_ok <- summ$rhat < 1.05
  summ$ess_ok  <- summ$ess_bulk > 100
  print(summ)

  rhat_max <- max(summ$rhat,     na.rm = TRUE)
  ess_min  <- min(summ$ess_bulk, na.rm = TRUE)
  cat(sprintf("  Rhat max: %.3f   ESS_bulk min: %.0f\n", rhat_max, ess_min))
  if (rhat_max > 1.05 || ess_min < 100) {
    cat("  WARNING: convergence flags raised\n")
  }

  ## Hold-out MAB / RMSE via posterior-median point prediction
  ho_mab <- NA; ho_rmse <- NA; ho_bias <- NA
  if (nrow(ho_valid) > 10) {
    draws_df <- as_draws_df(fit$draws())
    eta_ho   <- predict_eta(form, draws_df, ho_valid)
    pred_a   <- exp(eta_ho)
    resid    <- ho_valid$dg_obs_a - pred_a
    ho_bias  <- mean(resid,       na.rm = TRUE)
    ho_mab   <- mean(abs(resid),  na.rm = TRUE)
    ho_rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
    cat(sprintf("  Hold-out: bias=%.4f  MAB=%.4f  RMSE=%.4f  (n=%d)\n",
                ho_bias, ho_mab, ho_rmse, nrow(ho_valid)))
  }

  results[[form]] <- list(
    form       = form,
    stan_file  = stan_file,
    wall_min   = wall_min,
    rhat_max   = rhat_max,
    ess_min    = ess_min,
    summ       = summ,
    n_obs      = nrow(prep$data),
    n_holdout  = nrow(ho_valid),
    ho_bias    = ho_bias,
    ho_mab     = ho_mab,
    ho_rmse    = ho_rmse
  )
}

## ---- Summary across forms -------------------------------------------------
cat("\n==============================================================\n")
cat("Pre-flight summary\n")
cat("==============================================================\n")
cmp <- rbindlist(lapply(results, function(r) data.table(
  form     = r$form,
  n_obs    = r$n_obs,
  wall_min = round(r$wall_min, 1),
  rhat_max = round(r$rhat_max, 3),
  ess_min  = round(r$ess_min, 0),
  ho_mab   = round(r$ho_mab,  4),
  ho_rmse  = round(r$ho_rmse, 4),
  ho_bias  = round(r$ho_bias, 4)
)))
print(cmp)

stamp   <- format(Sys.time(), "%Y%m%d_%H%M")
out_rds <- file.path(opts$out,
                     sprintf("preflight_dg_%s_n%d_%s.rds",
                             opts$form, opts$n, stamp))
saveRDS(list(options = opts, results = results, compare = cmp), out_rds)
fwrite(cmp, sub("\\.rds$", "_compare.csv", out_rds))
cat("\nSaved:", out_rds, "\n")

cat("==============================================================\n")
cat("Pre-flight complete:", format(Sys.time()), "\n")
cat("==============================================================\n")
