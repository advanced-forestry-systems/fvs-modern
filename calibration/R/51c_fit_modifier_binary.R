## =============================================================================
## 51c_fit_modifier_binary.R
##
## Phase 3 Stage B for binary outcomes (mortality). Reads the bundle
## produced by 50b_extract_mort_modifier_data.R and fits
## modifier_binary.stan to recover plantation, disturbance, and
## treatment modifier coefficients on the logit scale.
##
## Posterior-mean alpha_* coefficients plug into the projection engine
## as additive eta corrections:
##   logit(p_total) = eta_base + log_years + delta(alpha_hat)
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

opts <- OptionParser(option_list = list(
  make_option("--bundle",    type = "character"),
  make_option("--stan_file", type = "character",
              default = "calibration/stan/modifier_binary.stan"),
  make_option("--lambda",    type = "integer", default = 10L,
              help = "5 | 10 | 20 — decay envelope to use"),
  make_option("--chains",    type = "integer", default = 2L),
  make_option("--warmup",    type = "integer", default = 500L),
  make_option("--sampling",  type = "integer", default = 500L),
  make_option("--seed",      type = "integer", default = 42L),
  make_option("--max_treedepth", type = "integer", default = 10L),
  make_option("--adapt_delta",   type = "double",  default = 0.9),
  make_option("--save_draws", action = "store_true", default = FALSE),
  make_option("--out",       type = "character"),
  make_option("--subsample", type = "integer", default = 0L)
)) |> parse_args()

stopifnot(!is.null(opts$bundle), file.exists(opts$bundle))
stopifnot(file.exists(opts$stan_file))
stopifnot(!is.null(opts$out))
stopifnot(opts$lambda %in% c(5L, 10L, 20L))

cat("==============================================================\n")
cat("FVS-CONUS Phase 3 Stage B (binary, modifier_binary.stan)\n")
cat("  bundle    :", opts$bundle, "\n")
cat("  stan_file :", opts$stan_file, "\n")
cat("  lambda    :", opts$lambda, "yr\n")
cat("  chains    :", opts$chains, "  warmup:", opts$warmup,
    "  sampling:", opts$sampling, "\n")
cat("  out       :", opts$out, "\n")
cat("==============================================================\n")

bundle <- readRDS(opts$bundle)
stopifnot(bundle$family == "binary")
cat(sprintf("Bundle: model=%s  family=%s  %s rows\n",
            bundle$model, bundle$family,
            format(nrow(bundle$data), big.mark = ",")))

d <- bundle$data
setDT(d)

# Drop rows with non-finite eta_base (numerical edge cases)
before <- nrow(d)
d <- d[is.finite(eta_base) & is.finite(log_years)]
if (nrow(d) < before)
  cat(sprintf("Dropped %s rows with non-finite eta/log_years\n",
              format(before - nrow(d), big.mark = ",")))

if (opts$subsample > 0L && opts$subsample < nrow(d)) {
  set.seed(opts$seed)
  d <- d[sample(.N, opts$subsample)]
  cat(sprintf("Subsampled to %s rows\n",
              format(nrow(d), big.mark = ",")))
}

dstrb_col <- sprintf("dstrb_decay_%dyr", opts$lambda)
trt_col   <- sprintf("trt_decay_%dyr",   opts$lambda)
for (c in c(dstrb_col, trt_col)) {
  if (!c %in% names(d)) stop("bundle missing column ", c)
}

stan_data <- list(
  N_obs         = nrow(d),
  N_L1          = length(bundle$L1_levels),
  alive         = as.integer(d$alive),
  eta_base      = d$eta_base,
  log_years     = d$log_years,
  is_plantation = as.numeric(d$is_plantation),
  d_fire        = as.numeric(d$had_fire_t1     %||% rep(0, nrow(d))),
  d_insect      = as.numeric(d$had_insect_t1   %||% rep(0, nrow(d))),
  d_disease     = as.numeric(d$had_disease_t1  %||% rep(0, nrow(d))),
  d_wind        = as.numeric(d$had_wind_t1     %||% rep(0, nrow(d))),
  d_harvest     = as.numeric(d$had_harvest_t1  %||% rep(0, nrow(d))),
  dstrb_decay   = d[[dstrb_col]],
  t_cutting     = as.numeric(d$had_cutting_t1   %||% rep(0, nrow(d))),
  t_site_prep   = as.numeric(d$had_site_prep_t1 %||% rep(0, nrow(d))),
  trt_decay     = d[[trt_col]],
  L1_idx        = as.integer(d$L1_idx)
)

cat("\nStan data ready:\n")
cat(sprintf("  N_obs = %s, N_L1 = %d, alive_frac = %.3f\n",
            format(stan_data$N_obs, big.mark = ","), stan_data$N_L1,
            mean(stan_data$alive)))
cat(sprintf("  eta_base mean = %.3f, sd = %.3f\n",
            mean(stan_data$eta_base), sd(stan_data$eta_base)))

mod <- cmdstan_model(opts$stan_file)

init_fn <- function(chain_id) {
  set.seed(opts$seed + chain_id)
  list(
    alpha_0        = rnorm(1, 0, 0.01),
    alpha_plant    = rnorm(1, 0, 0.05),
    alpha_fire     = rnorm(1, 0, 0.05),
    alpha_insect   = rnorm(1, 0, 0.05),
    alpha_disease  = rnorm(1, 0, 0.05),
    alpha_wind     = rnorm(1, 0, 0.05),
    alpha_harvest  = rnorm(1, 0, 0.05),
    alpha_cutting  = rnorm(1, 0, 0.05),
    alpha_siteprep = rnorm(1, 0, 0.05),
    z_L1_raw       = rnorm(stan_data$N_L1, 0, 0.05),
    sigma_L1       = 0.05
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

coef_vars <- c("alpha_0", "alpha_plant",
               "alpha_fire", "alpha_insect", "alpha_disease",
               "alpha_wind", "alpha_harvest",
               "alpha_cutting", "alpha_siteprep",
               "sigma_L1")
summ <- fit$summary(variables = coef_vars)
print(summ)

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
tag <- sprintf("%s_modifier_lambda%d", bundle$model, opts$lambda)
write.csv(summ,
          file.path(opts$out, sprintf("%s_summary.csv", tag)),
          row.names = FALSE)

saveRDS(list(
  bundle_meta = list(model = bundle$model, family = bundle$family,
                     response = bundle$response, fit_path = bundle$fit_path),
  lambda      = opts$lambda,
  N_obs       = stan_data$N_obs,
  N_L1        = stan_data$N_L1,
  L1_levels   = bundle$L1_levels,
  summary     = summ,
  wall_min    = as.numeric(t1 - t0, units = "mins")
), file.path(opts$out, sprintf("%s_meta.rds", tag)), compress = "xz")

if (opts$save_draws) {
  fit$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))
}

cat(sprintf("\nDone. Artifacts in %s\n", opts$out))
