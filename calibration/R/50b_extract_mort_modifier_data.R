## =============================================================================
## 50b_extract_mort_modifier_data.R
##
## Phase 3 Stage A for the BINARY mortality base fit. Produces the
## modifier-ready bundle that 51c_fit_modifier_binary.R consumes.
##
## Unlike 50_extract_base_residuals.R (which serializes a continuous
## residual = response - eta_base), this script carries the binary
## alive/dead indicator and the precomputed eta_base separately so the
## modifier model can use a Bernoulli logit likelihood. This preserves
## the correct binomial uncertainty rather than approximating logit(0/1)
## as a continuous residual.
##
## Bundle schema:
##   $model      = "mort_logit_simple"
##   $family     = "binary"
##   $sp_levels, $L1_levels
##   $data: data.table with eta_base, log_years, alive, modifier covariates
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

opts <- OptionParser(option_list = list(
  make_option("--fit",   type = "character"),
  make_option("--meta",  type = "character"),
  make_option("--pairs", type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond.rds"),
  make_option("--out",   type = "character")
)) |> parse_args()

stopifnot(!is.null(opts$fit), file.exists(opts$fit))
stopifnot(!is.null(opts$meta), file.exists(opts$meta))
stopifnot(file.exists(opts$pairs))
stopifnot(!is.null(opts$out))

cat("==============================================================\n")
cat("FVS-CONUS Phase 3 Stage A (binary mortality)\n")
cat("  fit   :", opts$fit, "\n")
cat("  meta  :", opts$meta, "\n")
cat("  pairs :", opts$pairs, "\n")
cat("  out   :", opts$out, "\n")
cat("==============================================================\n")

meta <- readRDS(opts$meta)
fit  <- readRDS(opts$fit)

sp_levels <- meta$prep_meta$sp
L1_levels <- meta$prep_meta$L1
cspi_shift <- meta$prep_meta$cspi_shift
cat(sprintf("Levels from meta: N_sp=%d  N_L1=%d  cspi_shift=%.4f\n",
            length(sp_levels), length(L1_levels), cspi_shift))

## ---- 1. Extract posterior means via fit$draws (avoid slow $summary) ------
get_mean <- function(name) {
  draws <- fit$draws(variables = name, format = "draws_matrix")
  mean(as.numeric(draws))
}
vec_mean <- function(name) {
  draws <- fit$draws(variables = name, format = "draws_matrix")
  colMeans(draws)
}

m0 <- get_mean("m0"); m1 <- get_mean("m1"); m2 <- get_mean("m2")
m3 <- get_mean("m3"); m4 <- get_mean("m4"); m5 <- get_mean("m5"); m6 <- get_mean("m6")
sigma_sp <- get_mean("sigma_sp"); sigma_L1 <- get_mean("sigma_L1")

z_sp_raw <- vec_mean("z_sp_raw")
z_L1_raw <- vec_mean("z_L1_raw")
z_sp <- as.numeric(z_sp_raw) * sigma_sp
z_L1 <- as.numeric(z_L1_raw) * sigma_L1
if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
if (length(z_L1) < length(L1_levels)) z_L1 <- rep_len(z_L1, length(L1_levels))

# Trait loadings — present only if traits=1 in the base fit
P_trait <- tryCatch(length(vec_mean("gamma")), error = function(e) 0L)
if (P_trait > 0L) {
  gamma_vec <- vec_mean("gamma")
  W <- meta$W
  if (is.null(W)) W <- matrix(0, nrow = length(sp_levels), ncol = P_trait)
  trait_effect <- as.numeric(W %*% gamma_vec)
  if (length(trait_effect) < length(sp_levels))
    trait_effect <- rep_len(trait_effect, length(sp_levels))
} else {
  trait_effect <- rep(0, length(sp_levels))
}

## ---- 2. Filter pairs to mortality plausibility envelope -------------------
pairs <- readRDS(opts$pairs)
setDT(pairs)
d <- pairs[
  !is.na(STATUS2) & STATUS2 %in% c(1L, 2L) &
  !is.na(DBH1) & DBH1 > 0 &
  !is.na(CR1)  & CR1 > 0  & CR1 <= 1.0 &
  !is.na(BA1)  & BA1 > 0 &
  !is.na(BAL1) & BAL1 >= 0 &
  !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
  !is.na(EPA_L1_CODE) &
  !is.na(cspi) & cspi > 0 &
  YEARS >= 1 & YEARS <= 20
]
d[, alive       := as.integer(STATUS2 == 1L)]
d[, dbh2        := DBH1^2]
d[, bal_over_ba := BAL1 / BA1]
d[, sqrt_ba_rd  := sqrt(BA1 * rd_add)]
d[, log_years   := log(YEARS)]
d[, ln_cspi_shift := log(cspi + cspi_shift)]

# Map pairs to fit level indices
d[, sp_idx := match(SPCD, sp_levels)]
d[, L1_idx := match(EPA_L1_CODE, L1_levels)]
d <- d[!is.na(sp_idx) & !is.na(L1_idx)]
cat(sprintf("After plausibility + level alignment: %s rows  (alive frac %.3f)\n",
            format(nrow(d), big.mark = ","), mean(d$alive)))

## ---- 3. Compute eta_base per observation ---------------------------------
eta <- m0 +
  trait_effect[d$sp_idx] +
  z_sp[d$sp_idx] +
  z_L1[d$L1_idx] +
  m1 * d$DBH1 +
  m2 * d$dbh2 +
  m3 * d$bal_over_ba +
  m4 * d$CR1 +
  m5 * d$sqrt_ba_rd +
  m6 * d$ln_cspi_shift

cat(sprintf("eta_base summary: mean=%.3f  sd=%.3f  p01=%.2f  p99=%.2f\n",
            mean(eta), sd(eta),
            quantile(eta, 0.01), quantile(eta, 0.99)))

## ---- 4. Slim and write the bundle ---------------------------------------
keep_cols <- c(
  "SPCD", "sp_idx", "EPA_L1_CODE", "L1_idx",
  "YEARS", "log_years", "alive",
  "is_plantation",
  "had_fire_t1", "had_insect_t1", "had_disease_t1",
  "had_wind_t1", "had_harvest_t1",
  "had_cutting_t1", "had_site_prep_t1",
  "years_since_dstrb", "years_since_trt",
  "dstrb_decay_5yr", "dstrb_decay_10yr", "dstrb_decay_20yr",
  "trt_decay_5yr", "trt_decay_10yr", "trt_decay_20yr"
)
keep_cols <- intersect(keep_cols, names(d))
out <- d[, ..keep_cols]
out[, eta_base := eta]

# Coerce indicator NAs to 0 / decay NAs to 0
ind_cols <- grep("^(is_plantation|had_)", names(out), value = TRUE)
for (c in ind_cols) {
  v <- out[[c]]; v[is.na(v)] <- 0L; out[[c]] <- as.integer(v)
}
decay_cols <- grep("_decay_", names(out), value = TRUE)
for (c in decay_cols) {
  v <- out[[c]]; v[is.na(v)] <- 0; out[[c]] <- v
}

cat(sprintf("Writing %s ...\n", opts$out))
dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  model      = "mort_logit_simple",
  family     = "binary",
  response   = "alive",
  fit_path   = opts$fit,
  meta_path  = opts$meta,
  pairs_path = opts$pairs,
  sp_levels  = sp_levels,
  L1_levels  = L1_levels,
  data       = out
), opts$out, compress = "xz")
cat(sprintf("Done. %s rows, %d cols\n",
            format(nrow(out), big.mark = ","), ncol(out)))
