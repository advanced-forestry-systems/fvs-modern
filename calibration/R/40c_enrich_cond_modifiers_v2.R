## =============================================================================
## 40c_enrich_cond_modifiers_v2.R
##
## Phase 2 v2: continuous-time treatment of disturbance and management.
## Reads the v1 COND-enriched pairs and adds:
##
##   1. Mid-interval imputation of DSTRBYR / TRTYR when the code is
##      recorded but the year is missing (t := YEARS / 2)
##   2. Acute decay kernel (current) and gamma kernel (peak-and-decay)
##      at multiple time constants
##   3. Plantation age taper using STDAGE so plantation effect attenuates
##      as the stand matures
##   4. Stacked disturbance support: DSTRBCD2 / DSTRBCD3 indicators with
##      their own time-since columns and kernels
##
## Output columns appended (kept side-by-side with v1 columns):
##
##   imputed_dstrb1_flag, imputed_trt1_flag         indicator imputation
##   years_since_dstrb1_imp, years_since_trt1_imp   imputed time
##   acute_kernel_5yr, acute_kernel_10yr             exp(-t/tau)
##   gamma_kernel_3yr, gamma_kernel_5yr              (t/tau)exp(1-t/tau)
##   plant_age_taper_30, plant_age_taper_50          sigmoid(STDAGE)
##   d_event_2_*, d_event_3_*                        DSTRBCD2 / 3 indicators
##   years_since_dstrb_2, years_since_dstrb_3        time since older events
##   acute_kernel_5yr_2, acute_kernel_5yr_3          decay envelopes for stacked
##
## Output file: conus_remeasurement_pairs_metric_cond_v2.rds
##
## Design note: kernels are applied as data so that the modifier Stan
## model only needs to multiply pre-computed columns. Pivoting to
## different kernel shapes is a one-line change here, not a model edit.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--in_rds",  type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond.rds"),
  make_option("--out_rds", type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"),
  make_option("--tau_acute_short", type = "double", default = 5),
  make_option("--tau_acute_long",  type = "double", default = 10),
  make_option("--tau_gamma_short", type = "double", default = 3),
  make_option("--tau_gamma_long",  type = "double", default = 5),
  make_option("--tau_plant_30",    type = "double", default = 30),
  make_option("--tau_plant_50",    type = "double", default = 50),
  make_option("--scale_plant",     type = "double", default = 10)
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS Phase 2 v2 enrichment (continuous-time kernels)\n")
cat("  in :", opts$in_rds, "\n")
cat("  out:", opts$out_rds, "\n")
cat("==============================================================\n")

stopifnot(file.exists(opts$in_rds))
d <- readRDS(opts$in_rds)
setDT(d)
cat(sprintf("Input: %s rows, %d cols\n",
            format(nrow(d), big.mark = ","), ncol(d)))

## ---- 1. Mid-interval year imputation -----------------------------------

# When code is recorded but year is missing, set t := YEARS / 2 (the
# midpoint of the measurement interval). Marker indicator preserved.

clip_yst <- function(x) pmin(pmax(x, 0), 50)

if (all(c("DSTRBCD1_cond1", "DSTRBYR1_cond1", "INVYR1", "YEARS") %in% names(d))) {
  d[, imputed_dstrb1_flag := as.integer(
    !is.na(DSTRBCD1_cond1) & DSTRBCD1_cond1 > 0 &
    (is.na(DSTRBYR1_cond1) | DSTRBYR1_cond1 <= 0)
  )]
  d[, years_since_dstrb1_imp := fifelse(
    imputed_dstrb1_flag == 1L,
    clip_yst(YEARS / 2),
    clip_yst(INVYR1 - DSTRBYR1_cond1)
  )]
  # NA when no event at all
  d[is.na(DSTRBCD1_cond1) | DSTRBCD1_cond1 <= 0,
    years_since_dstrb1_imp := NA_real_]
  cat(sprintf("DSTRBCD1 events: %s, of which year-imputed: %s\n",
              format(sum(!is.na(d$DSTRBCD1_cond1) & d$DSTRBCD1_cond1 > 0),
                     big.mark = ","),
              format(sum(d$imputed_dstrb1_flag, na.rm = TRUE),
                     big.mark = ",")))
}

if (all(c("TRTCD1_cond1", "TRTYR1_cond1", "INVYR1", "YEARS") %in% names(d))) {
  d[, imputed_trt1_flag := as.integer(
    !is.na(TRTCD1_cond1) & TRTCD1_cond1 > 0 &
    (is.na(TRTYR1_cond1) | TRTYR1_cond1 <= 0)
  )]
  d[, years_since_trt1_imp := fifelse(
    imputed_trt1_flag == 1L,
    clip_yst(YEARS / 2),
    clip_yst(INVYR1 - TRTYR1_cond1)
  )]
  d[is.na(TRTCD1_cond1) | TRTCD1_cond1 <= 0,
    years_since_trt1_imp := NA_real_]
  cat(sprintf("TRTCD1 events: %s, of which year-imputed: %s\n",
              format(sum(!is.na(d$TRTCD1_cond1) & d$TRTCD1_cond1 > 0),
                     big.mark = ","),
              format(sum(d$imputed_trt1_flag, na.rm = TRUE),
                     big.mark = ",")))
}

## ---- 2. Acute and gamma kernels for the primary event ------------------

# Acute (event-instant peak, monotone decay)
acute_kernel <- function(t, tau) {
  ifelse(is.na(t), 0, exp(-t / tau))
}
# Gamma (peaks at t = tau, asymmetric)
# Normalized so peak == 1.0 when t == tau:
#   peak value = exp(0) = 1
gamma_kernel <- function(t, tau) {
  ifelse(is.na(t), 0, (t / tau) * exp(1 - t / tau))
}

if ("years_since_dstrb1_imp" %in% names(d)) {
  d[, acute_kernel_5yr  := acute_kernel(years_since_dstrb1_imp, opts$tau_acute_short)]
  d[, acute_kernel_10yr := acute_kernel(years_since_dstrb1_imp, opts$tau_acute_long)]
  d[, gamma_kernel_3yr  := gamma_kernel(years_since_dstrb1_imp, opts$tau_gamma_short)]
  d[, gamma_kernel_5yr  := gamma_kernel(years_since_dstrb1_imp, opts$tau_gamma_long)]
}

if ("years_since_trt1_imp" %in% names(d)) {
  d[, trt_acute_kernel_5yr  := acute_kernel(years_since_trt1_imp, opts$tau_acute_short)]
  d[, trt_acute_kernel_10yr := acute_kernel(years_since_trt1_imp, opts$tau_acute_long)]
  d[, trt_gamma_kernel_3yr  := gamma_kernel(years_since_trt1_imp, opts$tau_gamma_short)]
  d[, trt_gamma_kernel_5yr  := gamma_kernel(years_since_trt1_imp, opts$tau_gamma_long)]
}

## ---- 3. Plantation age taper -------------------------------------------

# Sigmoid attenuation: plantation effect is full at young ages, drops to
# half at STDAGE = tau_plant, near-zero by tau_plant + 2*scale.
# Pure logistic.

plant_age_taper <- function(stdage, tau, scale) {
  ifelse(is.na(stdage), 1, 1 / (1 + exp((stdage - tau) / scale)))
}

if ("STDAGE_cond1" %in% names(d)) {
  d[, plant_age_taper_30 := plant_age_taper(STDAGE_cond1,
                                              opts$tau_plant_30,
                                              opts$scale_plant)]
  d[, plant_age_taper_50 := plant_age_taper(STDAGE_cond1,
                                              opts$tau_plant_50,
                                              opts$scale_plant)]
  cat(sprintf("Plantation age taper computed for %s rows with STDAGE\n",
              format(sum(!is.na(d$STDAGE_cond1)), big.mark = ",")))
}

## ---- 4. Stacked disturbances (event 2 and 3) ---------------------------

for (i in 2:3) {
  cd_col  <- sprintf("DSTRBCD%d_cond1", i)
  yr_col  <- sprintf("DSTRBYR%d_cond1", i)
  if (!cd_col %in% names(d) || !yr_col %in% names(d)) next

  imp_flag <- sprintf("imputed_dstrb%d_flag", i)
  yst_col  <- sprintf("years_since_dstrb%d_imp", i)
  ack_5    <- sprintf("acute_kernel_5yr_%d", i)
  ack_10   <- sprintf("acute_kernel_10yr_%d", i)
  gk_3     <- sprintf("gamma_kernel_3yr_%d", i)

  d[, (imp_flag) := as.integer(
    !is.na(get(cd_col)) & get(cd_col) > 0 &
    (is.na(get(yr_col)) | get(yr_col) <= 0)
  )]
  d[, (yst_col) := fifelse(
    get(imp_flag) == 1L,
    clip_yst(YEARS / 2),
    clip_yst(INVYR1 - get(yr_col))
  )]
  d[is.na(get(cd_col)) | get(cd_col) <= 0,
    (yst_col) := NA_real_]
  d[, (ack_5)  := acute_kernel(get(yst_col), opts$tau_acute_short)]
  d[, (ack_10) := acute_kernel(get(yst_col), opts$tau_acute_long)]
  d[, (gk_3)   := gamma_kernel(get(yst_col), opts$tau_gamma_short)]

  # Type indicators for stacked event
  d[, sprintf("had_fire_t1_%d",    i) :=
      as.integer(!is.na(get(cd_col)) & get(cd_col) == 30)]
  d[, sprintf("had_insect_t1_%d",  i) :=
      as.integer(!is.na(get(cd_col)) & get(cd_col) %in% 10:19)]
  d[, sprintf("had_disease_t1_%d", i) :=
      as.integer(!is.na(get(cd_col)) & get(cd_col) %in% 20:29)]
  d[, sprintf("had_wind_t1_%d",    i) :=
      as.integer(!is.na(get(cd_col)) & get(cd_col) == 50)]
  d[, sprintf("had_harvest_t1_%d", i) :=
      as.integer(!is.na(get(cd_col)) & get(cd_col) == 80)]

  cat(sprintf("Stacked event %d: %s rows have a code, %s of those imputed\n",
              i,
              format(sum(!is.na(d[[cd_col]]) & d[[cd_col]] > 0),
                     big.mark = ","),
              format(sum(d[[imp_flag]], na.rm = TRUE),
                     big.mark = ",")))
}

## Count of stacked events per row
d[, n_stacked_events := as.integer(!is.na(DSTRBCD1_cond1) & DSTRBCD1_cond1 > 0)]
if ("DSTRBCD2_cond1" %in% names(d)) {
  d[, n_stacked_events := n_stacked_events +
      as.integer(!is.na(DSTRBCD2_cond1) & DSTRBCD2_cond1 > 0)]
}
if ("DSTRBCD3_cond1" %in% names(d)) {
  d[, n_stacked_events := n_stacked_events +
      as.integer(!is.na(DSTRBCD3_cond1) & DSTRBCD3_cond1 > 0)]
}

cat("\nSummary of new continuous-time columns:\n")
for (v in c("imputed_dstrb1_flag", "imputed_trt1_flag",
            "years_since_dstrb1_imp", "years_since_trt1_imp",
            "acute_kernel_5yr", "acute_kernel_10yr",
            "gamma_kernel_3yr", "gamma_kernel_5yr",
            "trt_gamma_kernel_3yr",
            "plant_age_taper_30", "plant_age_taper_50",
            "n_stacked_events")) {
  if (v %in% names(d)) {
    x <- d[[v]]
    cat(sprintf("  %-26s n_nonzero=%10s  mean=%.4f  max=%.4f\n",
                v,
                format(sum(x > 0, na.rm = TRUE), big.mark = ","),
                mean(x, na.rm = TRUE),
                max(x, na.rm = TRUE)))
  }
}

cat(sprintf("\nWriting %s ...\n", opts$out_rds))
saveRDS(d, opts$out_rds, compress = "xz")
sz <- file.info(opts$out_rds)$size
cat(sprintf("Done. %s rows, %d cols, %.1f MB\n",
            format(nrow(d), big.mark = ","), ncol(d), sz / 1e6))
