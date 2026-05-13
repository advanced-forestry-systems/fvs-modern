## =============================================================================
## 51b_fit_modifier_traitmed.R
##
## Trait-mediated modifier fit. Reads the residual bundle produced by
## 50_extract_base_residuals.R, loads species traits, builds the W matrix
## aligned to the residual bundle's species set, and fits
## modifier_traitmed.stan.
##
## Output: same shape as 51_fit_modifier.R but with gamma_alpha_*
## coefficients in the summary so we can read whether plantation /
## disturbance / treatment effects vary with traits.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
})

DEFAULT_TRAIT_COLS <- c("wood_specific_gravity", "shade_tolerance_num",
                        "softwood", "leaf_longevity_months",
                        "max_ht_m", "max_dbh_cm",
                        "vulnerability_score", "sensitivity",
                        "sla_mm2_mg", "seed_mass_mg")

opts <- OptionParser(option_list = list(
  make_option("--residuals",  type = "character"),
  make_option("--stan_file",  type = "character",
              default = "calibration/stan/modifier_traitmed.stan"),
  make_option("--traits_rds", type = "character",
              default = "calibration/traits/species_traits.rds"),
  make_option("--lambda",     type = "integer", default = 10L),
  make_option("--chains",     type = "integer", default = 2L),
  make_option("--warmup",     type = "integer", default = 250L),
  make_option("--sampling",   type = "integer", default = 250L),
  make_option("--seed",       type = "integer", default = 42L),
  make_option("--max_treedepth", type = "integer", default = 8L),
  make_option("--adapt_delta",   type = "double",  default = 0.9),
  make_option("--save_draws", action = "store_true", default = FALSE),
  make_option("--out",        type = "character"),
  make_option("--subsample",  type = "integer", default = 30000L)
)) |> parse_args()

stopifnot(!is.null(opts$residuals), file.exists(opts$residuals))
stopifnot(!is.null(opts$out))
stopifnot(opts$lambda %in% c(5L, 10L, 20L))
stopifnot(file.exists(opts$stan_file))
stopifnot(file.exists(opts$traits_rds))

cat("==============================================================\n")
cat("FVS-CONUS Phase 3 trait-mediated modifier fit\n")
cat("  residuals :", opts$residuals, "\n")
cat("  stan_file :", opts$stan_file, "\n")
cat("  traits    :", opts$traits_rds, "\n")
cat("  lambda    :", opts$lambda, "yr\n")
cat("  out       :", opts$out, "\n")
cat("==============================================================\n")

bundle <- readRDS(opts$residuals)
cat(sprintf("Residual bundle: model=%s family=%s rows=%s\n",
            bundle$model, bundle$family,
            format(nrow(bundle$data), big.mark = ",")))

d <- bundle$data
setDT(d)
sp_levels <- bundle$sp_levels

before <- nrow(d)
d <- d[is.finite(residual)]
if (nrow(d) < before)
  cat(sprintf("Dropped %s non-finite rows\n",
              format(before - nrow(d), big.mark = ",")))

if (opts$subsample > 0L && opts$subsample < nrow(d)) {
  set.seed(opts$seed)
  d <- d[sample(.N, opts$subsample)]
  cat(sprintf("Subsampled to %s rows\n",
              format(nrow(d), big.mark = ",")))
}

## ---- Build W aligned to bundle species ---------------------------------

tr <- readRDS(opts$traits_rds)
setDT(tr)
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
    cat("  drop zero-variance trait:", cc, "\n")
    tr_sp[, (cc) := NULL]; next
  }
  tr_sp[, (cc) := (v - mean(v)) / s]
  informative <- c(informative, cc)
}
W <- as.matrix(tr_sp[, informative, with = FALSE])
P_trait <- ncol(W)
if (P_trait == 0) {
  W <- matrix(0, length(sp_levels), 1)
  P_trait <- 0
}
cat(sprintf("Traits kept: %s (P_trait = %d)\n",
            paste(informative, collapse = ", "), P_trait))

## ---- Pick decay envelope ----

dstrb_col <- sprintf("dstrb_decay_%dyr", opts$lambda)
trt_col   <- sprintf("trt_decay_%dyr",   opts$lambda)
for (c in c(dstrb_col, trt_col)) {
  if (!c %in% names(d))
    stop("residual bundle missing column ", c)
}

stan_data <- list(
  N_obs         = nrow(d),
  N_L1          = length(bundle$L1_levels),
  N_sp          = length(sp_levels),
  P_trait       = P_trait,
  residual      = d$residual,
  weight        = d$weight,
  is_plantation = as.numeric(d$is_plantation),
  d_fire        = as.numeric(d$had_fire_t1),
  d_insect      = as.numeric(d$had_insect_t1),
  d_disease     = as.numeric(d$had_disease_t1),
  d_wind        = as.numeric(d$had_wind_t1),
  d_harvest     = as.numeric(d$had_harvest_t1),
  dstrb_decay   = d[[dstrb_col]],
  t_cutting     = as.numeric(d$had_cutting_t1),
  t_site_prep   = as.numeric(d$had_site_prep_t1),
  trt_decay     = d[[trt_col]],
  L1_idx        = as.integer(d$L1_idx),
  sp_idx        = as.integer(d$sp_idx),
  W             = W
)

cat(sprintf("\nStan data: N_obs=%s  N_L1=%d  N_sp=%d  P_trait=%d\n",
            format(stan_data$N_obs, big.mark = ","),
            stan_data$N_L1, stan_data$N_sp, stan_data$P_trait))

mod <- cmdstan_model(opts$stan_file)

init_fn <- function(chain_id) {
  set.seed(opts$seed + chain_id)
  list(
    alpha_0       = 0,
    alpha_plant   = 0, alpha_fire = 0, alpha_insect = 0,
    alpha_disease = 0, alpha_wind = 0, alpha_harvest = 0,
    alpha_cutting = 0, alpha_siteprep = 0,
    gamma_alpha_plant    = rnorm(P_trait, 0, 0.01),
    gamma_alpha_fire     = rnorm(P_trait, 0, 0.01),
    gamma_alpha_insect   = rnorm(P_trait, 0, 0.01),
    gamma_alpha_disease  = rnorm(P_trait, 0, 0.01),
    gamma_alpha_wind     = rnorm(P_trait, 0, 0.01),
    gamma_alpha_harvest  = rnorm(P_trait, 0, 0.01),
    gamma_alpha_cutting  = rnorm(P_trait, 0, 0.01),
    gamma_alpha_siteprep = rnorm(P_trait, 0, 0.01),
    z_L1_raw    = rnorm(stan_data$N_L1, 0, 0.05),
    sigma_L1    = 0.05,
    sigma_resid = max(0.05, sd(stan_data$residual, na.rm = TRUE) * 0.9)
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

key_global <- c("alpha_0", "alpha_plant",
                "alpha_fire","alpha_insect","alpha_disease",
                "alpha_wind","alpha_harvest",
                "alpha_cutting","alpha_siteprep",
                "sigma_L1","sigma_resid")
gamma_vars <- c(paste0("gamma_alpha_plant[",   seq_len(P_trait), "]"),
                paste0("gamma_alpha_insect[",  seq_len(P_trait), "]"),
                paste0("gamma_alpha_cutting[", seq_len(P_trait), "]"),
                paste0("gamma_alpha_harvest[", seq_len(P_trait), "]"))

summ_g <- fit$summary(variables = key_global)
cat("\n=== global modifier alphas ===\n")
print(summ_g)

if (P_trait > 0) {
  summ_t <- fit$summary(variables = gamma_vars)
  cat("\n=== trait-mediated gammas (plant, insect, cutting, harvest) ===\n")
  summ_t$trait <- rep(informative, 4)
  print(summ_t)
}

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
tag <- sprintf("%s_traitmed_lambda%d", bundle$model, opts$lambda)
write.csv(summ_g, file.path(opts$out, sprintf("%s_global_summary.csv", tag)),
          row.names = FALSE)
if (P_trait > 0)
  write.csv(summ_t, file.path(opts$out, sprintf("%s_gamma_summary.csv", tag)),
            row.names = FALSE)
saveRDS(list(
  bundle_meta = list(model = bundle$model, family = bundle$family,
                     response = bundle$response, fit_path = bundle$fit_path),
  lambda      = opts$lambda,
  N_obs       = stan_data$N_obs,
  N_L1        = stan_data$N_L1,
  N_sp        = stan_data$N_sp,
  P_trait     = stan_data$P_trait,
  traits_kept = informative,
  global      = summ_g,
  gamma       = if (P_trait > 0) summ_t else NULL,
  wall_min    = as.numeric(t1 - t0, units = "mins")
), file.path(opts$out, sprintf("%s_meta.rds", tag)), compress = "xz")

if (opts$save_draws) {
  fit$save_object(file.path(opts$out, sprintf("%s_fit.rds", tag)))
}

cat(sprintf("\nDone. Artifacts in %s\n", opts$out))
