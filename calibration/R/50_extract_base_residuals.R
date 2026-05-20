## =============================================================================
## 50_extract_base_residuals.R
##
## Phase 3 Stage A. Given a base-model Stan fit (from 31_fit_dg_organon.R,
## 33_fit_hcb_organon.R, etc.) and the enriched pairs RDS, compute
## per-observation eta_base posterior mean and the residual on the scale
## the base model uses. Save a compact RDS that the Phase 3 modifier
## fit (51_fit_modifier.R) will consume.
##
## Supports four response families:
##   --family log        residual = log(obs)   - eta_base     (DG, HG)
##   --family logit      residual = logit(obs) - eta_base     (HCB, CR, mort)
##   --family identity   residual = obs        - eta_base
##
## The residual output carries only the columns needed by the modifier
## fit, so downstream scripts do not have to re-read the 466 MB metric
## pairs file.
##
## Usage (DG ORGANON example):
##   Rscript calibration/R/50_extract_base_residuals.R \
##       --fit      calibration/output/conus/dg/dg_organon_cspi_traits1_fit.rds \
##       --pairs    calibration/data/conus_remeasurement_pairs_metric_cond.rds \
##       --family   log \
##       --response hg_obs_a \
##       --model    dg_organon \
##       --out      calibration/output/conus/dg/dg_organon_residuals.rds
##
## Notes:
##   * The Stan models compile eta in the `model` block (not
##     `generated quantities`), so we reconstruct eta on the R side
##     using the posterior-mean point estimates of the fixed
##     coefficients and random effects. That's a good approximation
##     for Stage A; the integrated Phase 4 re-fit will handle full
##     posterior uncertainty.
##   * For DG ORGANON: eta = a0 + trait_effect[sp] + z_sp[sp] + z_L1
##     + z_L2 + z_L3 + a1*ln(DBH+K1) + a2*DBH^K2 + a3*ln_cr_adj +
##     a4*ln_cspi_shift + a5*comp + a6*sqrt_ba + a7*ba_x_rd +
##     a8*ccfl_x_rd
##   * Keeps the raw observed response alongside so downstream QC can
##     re-check the transform.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--fit",      type = "character"),
  make_option("--pairs",    type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond.rds"),
  make_option("--family",   type = "character", default = "log",
              help = "log | logit | identity"),
  make_option("--response", type = "character",
              help = "name of the observed response column in pairs"),
  make_option("--model",    type = "character",
              help = "model id: dg_organon | hcb_organon | hg_organon | ..."),
  make_option("--out",      type = "character"),
  make_option("--subsample", type = "integer", default = 0L,
              help = "Optional row subsample; 0 = keep all that match")
)) |> parse_args()

stopifnot(!is.null(opts$fit),      file.exists(opts$fit))
stopifnot(!is.null(opts$pairs),    file.exists(opts$pairs))
stopifnot(!is.null(opts$response), !is.null(opts$out), !is.null(opts$model))
stopifnot(opts$family %in% c("log", "logit", "identity"))

cat("==============================================================\n")
cat("FVS-CONUS Phase 3 Stage A residual extraction\n")
cat("  fit     :", opts$fit, "\n")
cat("  pairs   :", opts$pairs, "\n")
cat("  family  :", opts$family, "\n")
cat("  response:", opts$response, "\n")
cat("  model   :", opts$model, "\n")
cat("  out     :", opts$out, "\n")
cat("==============================================================\n")

fit <- readRDS(opts$fit)
pairs <- readRDS(opts$pairs)
setDT(pairs)

cat(sprintf("Pairs: %s rows, %d cols\n",
            format(nrow(pairs), big.mark = ","), ncol(pairs)))

## ---- 1. Recover species / L1 / L2 / L3 level vectors ---------------------

# The base-fit RDS carries prep_meta (from *_meta.rds) but fit.rds is just the
# CmdStanMCMC object. Look for the companion _meta.rds next to --fit.
meta_path <- sub("_fit\\.rds$", "_meta.rds", opts$fit)
if (!file.exists(meta_path))
  stop("Could not find companion meta file at ", meta_path,
       "\n  50_extract_base_residuals.R needs it for level vectors and K constants.")
meta <- readRDS(meta_path)

sp_levels <- meta$prep_meta$sp
L1_levels <- meta$prep_meta$L1
L2_levels <- meta$prep_meta$L2
L3_levels <- meta$prep_meta$L3

cat(sprintf("Levels from meta: N_sp=%d  N_L1=%d  N_L2=%d  N_L3=%d\n",
            length(sp_levels), length(L1_levels),
            length(L2_levels), length(L3_levels)))

## ---- Posterior-mean accessors (on-demand) -------------------------------
## Defined here at the top so all family branches (HT-DBH and CR compute eta
## inline; DG/HG/HCB compute eta later) can use them. On-demand fit$draws()
## per variable avoids the multi-minute full fit$summary() pass.
get_mean <- function(v) {
  out <- tryCatch(
    fit$draws(variables = v, format = "draws_matrix"),
    error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  mean(as.numeric(out))
}
vec_mean <- function(base, default_n = 0L) {
  out <- tryCatch(
    fit$draws(variables = base, format = "draws_matrix"),
    error = function(e) NULL)
  if (is.null(out)) return(rep(0, default_n))
  colMeans(out)
}

## ---- 2. Filter pairs to the same plausibility envelope the base fit used -

model_family <- tolower(opts$model)

if (startsWith(model_family, "dg_") &&
    !startsWith(model_family, "dg_kuehne") &&
    !startsWith(model_family, "dg_weiskittel") &&
    !startsWith(model_family, "dg_compound")) {
  d <- pairs[
    !is.na(DBH1) & DBH1 > 0 &
    !is.na(DBH2) & DBH2 > DBH1 &
    (DBH2 - DBH1) / YEARS > 0.001 & (DBH2 - DBH1) / YEARS < 5.0 &
    !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(CCFL1) & CCFL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    YEARS >= 1 & YEARS <= 20 &
    !is.na(cspi) & cspi > 0
  ]
  d[, hg_obs_a   := NULL]
  d[, dg_obs_a   := (DBH2 - DBH1) / YEARS]
  d[, ln_cr_adj  := log((CR1 + 0.2) / 1.2)]
  d[, sqrt_ba    := sqrt(BA1)]
  d[, ba_x_rd    := BA1   * rd_add]
  d[, ccfl_x_rd  := CCFL1 * rd_add]
  d[, sqrt_years := sqrt(YEARS)]
  cspi_shift <- meta$prep_meta$cspi_shift
  d[, ln_cspi_shift := log(cspi + cspi_shift)]
  K1 <- meta$K1h %||% meta$K1; K2 <- meta$K2h %||% meta$K2; K4 <- meta$K4h %||% meta$K4
  if (is.null(K1) || is.null(K2) || is.null(K4))
    stop("Meta file missing K1/K2/K4 constants for DG")
  d[, ln_dbh_k1 := log(DBH1 + K1)]
  d[, dbh_k2    := DBH1^K2]
  d[, comp      := BAL1 / log(DBH1 + K4)]
} else if (startsWith(model_family, "dg_kuehne") || startsWith(model_family, "dgkuehne")) {
  ## Kuehne 2022 form. Uses softwood/hardwood BAL split, log-DBH (no K1
  ## offset), and direct DBH coefficient (no K2 power).
  d <- pairs[
    !is.na(DBH1) & DBH1 > 0 &
    !is.na(DBH2) & DBH2 > DBH1 &
    (DBH2 - DBH1) / YEARS > -0.5 & (DBH2 - DBH1) / YEARS < 5.0 &
    !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(BAL_SW1) & BAL_SW1 >= 0 &
    !is.na(BAL_HW1) & BAL_HW1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    YEARS >= 1 & YEARS <= 20 &
    !is.na(cspi) & cspi > 0
  ]
  d[, dg_obs_a    := (DBH2 - DBH1) / YEARS]
  d[, ln_dbh      := log(DBH1)]
  d[, ln_cr_adj   := log((CR1 + 0.2) / 1.2)]
  d[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
  d[, ba_x_rd     := BA1 * rd_add]
  d[, bal_x_rd    := BAL1 * rd_add]
  d[, sqrt_years  := sqrt(YEARS)]
  cspi_shift <- meta$prep_meta$cspi_shift
  d[, ln_cspi_shift := log(cspi + cspi_shift)]
  ## eta / residual computed in Section 5 after sp_idx etc are populated.
} else if (startsWith(model_family, "hcb_")) {
  d <- pairs[
    !is.na(HT1) & HT1 >= 1.37 &
    !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
    !is.na(DBH1) & DBH1 > 0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    !is.na(cspi) & cspi > 0
  ]
  d[, hcb1  := (1 - CR1) * HT1]
  d[, ratio := pmin(pmax(hcb1 / HT1, 1e-4), 1 - 1e-4)]
  d[, ln_ht       := log(HT1)]
  d[, ln_dbh      := log(DBH1)]
  d[, bal_over_ht := BAL1 / (HT1 + 1)]
  d[, sqrt_ba     := sqrt(BA1)]
  cspi_shift <- meta$prep_meta$cspi_shift
  d[, ln_cspi_shift := log(cspi + cspi_shift)]
} else if (startsWith(model_family, "hg_")) {
  d <- pairs[
    !is.na(HT1) & HT1 >= 1.37 &
    !is.na(HT2) & HT2 > HT1 &
    (HT2 - HT1) / YEARS > 0.001 & (HT2 - HT1) / YEARS < 2.0 &
    !is.na(CR1) & CR1 > 0 & CR1 <= 1.0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(CCFL1) & CCFL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    YEARS >= 1 & YEARS <= 20 &
    !is.na(cspi) & cspi > 0
  ]
  d[, hg_obs_a  := (HT2 - HT1) / YEARS]
  d[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
  d[, sqrt_ba   := sqrt(BA1)]
  d[, ba_x_rd   := BA1   * rd_add]
  d[, ccfl_x_rd := CCFL1 * rd_add]
  d[, sqrt_years := sqrt(YEARS)]
  cspi_shift <- meta$prep_meta$cspi_shift
  d[, ln_cspi_shift := log(cspi + cspi_shift)]
  K1h <- meta$K1h; K2h <- meta$K2h; K4h <- meta$K4h
  if (is.null(K1h) || is.null(K2h) || is.null(K4h))
    stop("Meta file missing K1h/K2h/K4h constants for HG")
  d[, ln_ht_k1 := log(HT1 + K1h)]
  d[, ht_k2    := HT1^K2h]
  d[, comp     := CCFL1 / log(HT1 + K4h)]
} else if (startsWith(model_family, "ht_dbh") || startsWith(model_family, "htdbh")) {
  d <- pairs[
    !is.na(HT1) & HT1 > 1.37 &
    !is.na(DBH1) & DBH1 > 0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) &
    !is.na(cspi) & cspi > 0
  ]
  d[, sqrt_ba       := sqrt(BA1)]
  d[, ba_x_rd       := BA1 * rd_add]
  d[, bal_x_rd      := BAL1 * rd_add]
  cspi_shift <- meta$prep_meta$cspi_shift
  d[, ln_cspi_shift := log(cspi + cspi_shift)]
  ## eta / residual computed in Section 5 after sp_idx / L1_idx are populated.
} else if (startsWith(model_family, "cr_") || startsWith(model_family, "crown")) {
  pairs[, HCB1_tmp := HT1 * (1 - CR1)]
  pairs[, HCB2_tmp := HT2 * (1 - CR2)]
  d <- pairs[
    !is.na(HT1) & HT1 >= 1.37 &
    !is.na(HT2) & HT2 >= 1.37 &
    !is.na(CR1) & CR1 > 0 & CR1 < 1.0 &
    !is.na(CR2) & CR2 > 0 & CR2 < 1.0 &
    HCB1_tmp > 0 & HCB2_tmp > 0 &
    !is.na(BA1) & BA1 > 0 &
    !is.na(BAL1) & BAL1 >= 0 &
    !is.na(rd_add) & rd_add > 0 & rd_add < 2.0 &
    !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
    !is.na(cspi) & cspi > 0
  ]
  d[, delta_hcb := pmax(HCB2_tmp - HCB1_tmp, 0)]
  d[, max_dhcb  := HT2 - HCB1_tmp]
  d <- d[max_dhcb > 0]
  d[, cr_ratio := pmin(pmax(delta_hcb / max_dhcb, 1e-3), 1 - 1e-3)]
  d[, ln_cr      := log(CR1)]
  d[, sqrt_ba    := sqrt(BA1)]
  d[, ln_bal_ba  := log(BAL1 / BA1 + 1)]
  d[, cr_over_rd := CR1 / rd_add]
  cspi_shift <- meta$prep_meta$cspi_shift
  d[, ln_cspi_shift := log(cspi + cspi_shift)]
  ## eta / residual computed in Section 5 after sp_idx / L1_idx are populated.
} else {
  stop("Unsupported model family: ", opts$model)
}

cat(sprintf("After plausibility filter: %s rows\n",
            format(nrow(d), big.mark = ",")))

## ---- 3. Map each observation to its fit's level indices ------------------

# Pool rare species into SPCD = -1L to match the fit
rare_set <- setdiff(unique(d$SPCD), sp_levels)
if (length(rare_set) > 0) {
  d[SPCD %in% rare_set, SPCD := -1L]
}
d <- d[SPCD %in% sp_levels]
d <- d[EPA_L1_CODE %in% L1_levels]
## Only filter / map L2 and L3 when the base fit actually used them.
## HT-DBH (and other 2-level fits) leave L2/L3 level vectors empty, so
## an unconditional filter would drop every row.
if (length(L2_levels) > 0) d <- d[EPA_L2_CODE %in% L2_levels]
if (length(L3_levels) > 0) d <- d[EPA_L3_CODE %in% L3_levels]

d[, sp_idx := match(SPCD,        sp_levels)]
d[, L1_idx := match(EPA_L1_CODE, L1_levels)]
if (length(L2_levels) > 0) {
  d[, L2_idx := match(EPA_L2_CODE, L2_levels)]
} else {
  d[, L2_idx := 1L]
}
if (length(L3_levels) > 0) {
  d[, L3_idx := match(EPA_L3_CODE, L3_levels)]
} else {
  d[, L3_idx := 1L]
}

cat(sprintf("After level alignment: %s rows\n",
            format(nrow(d), big.mark = ",")))

if (opts$subsample > 0L && opts$subsample < nrow(d)) {
  set.seed(42)
  d <- d[sample(.N, opts$subsample)]
  cat(sprintf("Subsampled to %s rows\n",
              format(nrow(d), big.mark = ",")))
}

## ---- 4. Pull posterior-mean point estimates ------------------------------
## Accessors get_mean / vec_mean already defined near top; nothing to do here.

## ---- 5. Compute eta_base per observation, by family ----------------------

if (startsWith(model_family, "dg_") &&
    !startsWith(model_family, "dg_kuehne") &&
    !startsWith(model_family, "dg_weiskittel") &&
    !startsWith(model_family, "dg_compound")) {
  a0 <- get_mean("a0"); a1 <- get_mean("a1"); a2 <- get_mean("a2")
  a3 <- get_mean("a3"); a4 <- get_mean("a4"); a5 <- get_mean("a5")
  a6 <- get_mean("a6"); a7 <- get_mean("a7"); a8 <- get_mean("a8")
  trait_effect <- vec_mean("trait_effect", default_n = length(sp_levels))
  z_sp <- vec_mean("z_sp")
  z_L1 <- vec_mean("z_L1"); z_L2 <- vec_mean("z_L2"); z_L3 <- vec_mean("z_L3")
  if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
  eta <- a0 +
    trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
    a1 * d$ln_dbh_k1 + a2 * d$dbh_k2 + a3 * d$ln_cr_adj +
    a4 * d$ln_cspi_shift + a5 * d$comp + a6 * d$sqrt_ba +
    a7 * d$ba_x_rd + a8 * d$ccfl_x_rd
  obs_raw <- d[[opts$response]]
  residual <- log(obs_raw) - eta
  weight <- d$sqrt_years          # normal(eta, sigma / sqrt_years) => weight = sqrt_years
} else if (startsWith(model_family, "hg_")) {
  a0 <- get_mean("a0"); a1 <- get_mean("a1"); a2 <- get_mean("a2")
  a3 <- get_mean("a3"); a4 <- get_mean("a4"); a5 <- get_mean("a5")
  a6 <- get_mean("a6"); a7 <- get_mean("a7"); a8 <- get_mean("a8")
  trait_effect <- vec_mean("trait_effect", default_n = length(sp_levels))
  z_sp <- vec_mean("z_sp")
  z_L1 <- vec_mean("z_L1"); z_L2 <- vec_mean("z_L2"); z_L3 <- vec_mean("z_L3")
  if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
  eta <- a0 +
    trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
    a1 * d$ln_ht_k1 + a2 * d$ht_k2 + a3 * d$ln_cr_adj +
    a4 * d$ln_cspi_shift + a5 * d$comp + a6 * d$sqrt_ba +
    a7 * d$ba_x_rd + a8 * d$ccfl_x_rd
  obs_raw <- d[[opts$response]]
  residual <- log(obs_raw) - eta
  weight <- d$sqrt_years
} else if (startsWith(model_family, "hcb_")) {
  h0 <- get_mean("h0"); h1 <- get_mean("h1"); h2 <- get_mean("h2")
  h3 <- get_mean("h3"); h4 <- get_mean("h4"); h5 <- get_mean("h5")
  z_sp <- vec_mean("z_sp")
  z_L1 <- vec_mean("z_L1"); z_L2 <- vec_mean("z_L2"); z_L3 <- vec_mean("z_L3")
  if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
  eta <- h0 + z_sp[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
    h1 * d$ln_ht + h2 * d$ln_dbh + h3 * d$bal_over_ht +
    h4 * d$sqrt_ba + h5 * d$ln_cspi_shift
  obs_raw <- d[[opts$response]]                  # ratio in (0, 1)
  residual <- qlogis(pmin(pmax(obs_raw, 1e-4), 1 - 1e-4)) - eta
  weight <- rep(1, nrow(d))
} else if (startsWith(model_family, "dg_kuehne") || startsWith(model_family, "dgkuehne")) {
  b0 <- get_mean("b0"); b1 <- get_mean("b1"); b2 <- get_mean("b2")
  b3 <- get_mean("b3"); b4 <- get_mean("b4"); b5 <- get_mean("b5")
  b6 <- get_mean("b6"); b7 <- get_mean("b7"); b8 <- get_mean("b8")
  trait_effect <- vec_mean("trait_effect", default_n = length(sp_levels))
  z_sp <- vec_mean("z_sp")
  z_L1 <- vec_mean("z_L1"); z_L2 <- vec_mean("z_L2"); z_L3 <- vec_mean("z_L3")
  if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
  eta <- b0 +
    trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
    z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
    b1 * d$ln_dbh + b2 * d$DBH1 + b3 * d$ln_cr_adj +
    b4 * d$ln_bal_sw_adj + b5 * d$BAL_HW1 + b6 * d$ln_cspi_shift +
    b7 * d$ba_x_rd + b8 * d$bal_x_rd
  obs_raw <- d$dg_obs_a
  ## lognormal variant: residual on log scale matches HG / DG ORG ln pattern.
  residual <- log(pmax(obs_raw, 0.001)) - eta
  weight <- d$sqrt_years
} else if (startsWith(model_family, "ht_dbh") || startsWith(model_family, "htdbh")) {
  b0 <- get_mean("b0"); b1 <- get_mean("b1")
  a_bal <- get_mean("a_bal"); a_ba <- get_mean("a_ba")
  a_cspi <- get_mean("a_cspi")
  a_bard <- get_mean("a_bard"); a_blrd <- get_mean("a_blrd")
  trait_effect <- vec_mean("trait_effect", default_n = length(sp_levels))
  z_sp <- vec_mean("z_sp")
  z_L1 <- vec_mean("z_L1")
  if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
  eta <- b0 +
    trait_effect[d$sp_idx] + z_sp[d$sp_idx] + z_L1[d$L1_idx] +
    a_bal  * d$BAL1   + a_ba   * d$sqrt_ba +
    a_cspi * d$ln_cspi_shift +
    a_bard * d$ba_x_rd + a_blrd * d$bal_x_rd +
    b1 / (d$DBH1 + 1.0)
  obs_raw  <- d$HT1
  residual <- log(pmax(obs_raw - 1.37, 0.01)) - eta
  weight   <- rep(1, nrow(d))
} else if (startsWith(model_family, "cr_") || startsWith(model_family, "crown")) {
  r0 <- get_mean("r0"); r1 <- get_mean("r1"); r2 <- get_mean("r2")
  r3 <- get_mean("r3"); r4 <- get_mean("r4"); r5 <- get_mean("r5")
  r6 <- get_mean("r6")
  trait_effect <- vec_mean("trait_effect", default_n = length(sp_levels))
  z_sp <- vec_mean("z_sp")
  z_L1 <- vec_mean("z_L1"); z_L2 <- vec_mean("z_L2"); z_L3 <- vec_mean("z_L3")
  if (length(z_sp) < length(sp_levels)) z_sp <- rep_len(z_sp, length(sp_levels))
  eta <- r0 + trait_effect[d$sp_idx] + z_sp[d$sp_idx] +
         z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] +
         r1 * d$ln_cr + r2 * d$CR1 + r3 * d$sqrt_ba +
         r4 * d$ln_bal_ba + r5 * d$cr_over_rd + r6 * d$ln_cspi_shift
  obs_raw <- d$cr_ratio
  residual <- -qlogis(pmin(pmax(obs_raw, 1e-4), 1 - 1e-4)) - eta
  weight   <- rep(1, nrow(d))
} else {
  stop("Unsupported model family: ", opts$model)
}

if (opts$family == "identity") {
  residual <- obs_raw - eta
}

cat(sprintf("Residual summary: n=%s  mean=%.4f  sd=%.4f  p01=%.3f  p99=%.3f\n",
            format(sum(is.finite(residual)), big.mark = ","),
            mean(residual, na.rm = TRUE),
            sd(residual, na.rm = TRUE),
            quantile(residual, 0.01, na.rm = TRUE),
            quantile(residual, 0.99, na.rm = TRUE)))

## ---- 6. Keep only what the modifier fit needs ---------------------------

keep_cols <- c(
  "SPCD", "sp_idx", "EPA_L1_CODE", "L1_idx",
  "YEARS",
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
out[, obs_raw  := obs_raw]
out[, residual := residual]
out[, weight   := weight]

# Coerce NAs on modifier indicators to 0 so they can be safely multiplied
ind_cols <- grep("^(is_plantation|had_)", names(out), value = TRUE)
for (c in ind_cols) {
  v <- out[[c]]
  v[is.na(v)] <- 0L
  out[[c]] <- as.integer(v)
}

decay_cols <- grep("_decay_", names(out), value = TRUE)
for (c in decay_cols) {
  v <- out[[c]]
  v[is.na(v)] <- 0
  out[[c]] <- v
}

cat(sprintf("\nWriting %s ...\n", opts$out))
dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  model      = opts$model,
  family     = opts$family,
  response   = opts$response,
  fit_path   = opts$fit,
  pairs_path = opts$pairs,
  sp_levels  = sp_levels,
  L1_levels  = L1_levels,
  data       = out
), opts$out, compress = "xz")
cat(sprintf("Done. %s rows, %d cols\n",
            format(nrow(out), big.mark = ","), ncol(out)))

`%||%` <- function(a, b) if (!is.null(a)) a else b
