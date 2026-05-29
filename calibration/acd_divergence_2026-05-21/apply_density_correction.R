# apply_density_correction.R  (v33, supersedes v31 coefficients)
# AcadianGY 12.3.9 bridge-level post-projection correction for the
# density-dependent residual identified in v30 and validated on a fresh
# out-of-sample run (v32). Coefficients are refit on the pooled n=184 sample.
# Cap is now asymmetric: the correction can only reduce BA_pred, never push it
# up. This protects high-density stands where the model is already near
# observed.
#
# 5-fold CV on the pooled n=184 sample, 50 random shuffles:
#   uncorrected:                BA bias = +14.46 percent, R^2 = 0.380
#   asymmetric (0, +25):        BA bias = +0.21 +/- 0.11, R^2 = 0.484 +/- 0.003
#   symmetric +/-25 (v31 ship): BA bias = +1.38 +/- 0.12, R^2 = 0.472 +/- 0.005
#
# Usage (downstream of AcadianGYOneStand):
#   source("apply_density_correction.R")
#   ba_corr <- apply_density_correction(BA_pred, BA_t1)
#
# Arguments:
#   BA_pred    numeric, predicted stand basal area at projection horizon (ft^2/ac)
#   BA_t1      numeric, observed stand basal area at start of projection (ft^2/ac)
#   upper_cap  numeric, max correction the function will subtract (default 25)
#   lower_cap  numeric, min correction (default 0; never push BA up)
#
# Returns a numeric vector of corrected BA predictions.

# Production coefficients, fit on pooled n=184 (v30 seed=42 + v32 seed=2027,
# both ME FIA, 10-yr remeasurement, AcadianGY 12.3.9 with MORTCAL=TRUE,
# CutPoint=0, CSI_SCALE=0.7).
ACD_DENSITY_CORRECTION <- list(
  a         = 36.9549,    # intercept (ft^2/ac)
  b         = -0.235987,  # slope on BA_t1 (per ft^2/ac)
  upper_cap = 25,         # max subtraction from BA_pred
  lower_cap = 0,          # min subtraction (0 = never increase BA_pred)
  crossover_BA = 156.6,   # BA_t1 where raw correction = 0
  n         = 184,        # n plots in calibration
  r2_fit    = 0.0956,     # in-sample fit R^2
  cv_bias_mean = 0.21,    # 5-fold CV mean BA bias percent
  cv_bias_sd   = 0.11,
  cv_r2_mean   = 0.484,
  cv_r2_sd     = 0.003,
  fitted_on  = "2026-05-28 v32+v30 pooled, n=184 ME FIA",
  applies_to = "ME FIA conditions, 10-yr remeasurement, 12.3.9 production posture",
  supersedes = "v31 coefficients (a=40.6, b=-0.334, symmetric +/-25 cap)"
)

apply_density_correction <- function(BA_pred, BA_t1,
                                     upper_cap = ACD_DENSITY_CORRECTION$upper_cap,
                                     lower_cap = ACD_DENSITY_CORRECTION$lower_cap) {
  stopifnot(length(BA_pred) == length(BA_t1))
  raw_correction <- ACD_DENSITY_CORRECTION$a + ACD_DENSITY_CORRECTION$b * BA_t1
  bounded <- pmax(lower_cap, pmin(upper_cap, raw_correction))
  BA_pred - bounded
}

apply_density_correction_verbose <- function(BA_pred, BA_t1,
                                             upper_cap = ACD_DENSITY_CORRECTION$upper_cap,
                                             lower_cap = ACD_DENSITY_CORRECTION$lower_cap) {
  raw <- ACD_DENSITY_CORRECTION$a + ACD_DENSITY_CORRECTION$b * BA_t1
  bounded <- pmax(lower_cap, pmin(upper_cap, raw))
  data.frame(
    BA_t1 = BA_t1,
    BA_pred_raw = BA_pred,
    raw_correction = raw,
    bounded_correction = bounded,
    BA_pred_corrected = BA_pred - bounded,
    hit_upper = raw > upper_cap,
    hit_lower = raw < lower_cap
  )
}

if (interactive() || sys.nframe() == 0) {
  cat("ACD_DENSITY_CORRECTION (v33):\n")
  for (k in names(ACD_DENSITY_CORRECTION))
    cat(sprintf("  %-18s %s\n", k, ACD_DENSITY_CORRECTION[[k]]))
  cat("\nSmoke test on representative BA_t1 values (default 0 to +25 cap):\n")
  test_BA_t1 <- c(10, 30, 60, 100, 156.6, 180, 220, 260)
  fake_BA_pred <- test_BA_t1 * 1.15
  d <- apply_density_correction_verbose(fake_BA_pred, test_BA_t1)
  print(format(d, digits = 4))
}
