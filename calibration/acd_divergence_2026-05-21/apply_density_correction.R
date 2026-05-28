# apply_density_correction.R
# AcadianGY 12.3.9 bridge-level post-projection correction for the
# density-dependent residual identified in v30 (n=93 ME FIA, R^2=0.188,
# p=1.4e-5). The model overshoots low-density stands and undershoots
# high-density ones; the post-hoc correction closes mean test bias from
# +10.9% to about +2.1% and lifts per-plot R^2 from 0.38 to 0.48 on a
# 200-iter 50/50 holdout.
#
# Usage (downstream of AcadianGYOneStand):
#   library(...)  # source this file
#   ba_corr <- apply_density_correction(BA_pred, BA_t1)
#
# Arguments:
#   BA_pred   numeric, predicted stand basal area at projection horizon (ft^2/ac)
#   BA_t1     numeric, observed stand basal area at start of projection (ft^2/ac)
#   cap       numeric, +/- clamp on the correction in ft^2/ac (default 25)
#
# Returns a numeric vector of corrected BA predictions.

# Production coefficients, fit on the full v30 sample (n=93, ME FIA, 10-yr
# remeasurement, AcadianGY 12.3.9 with MORTCAL=TRUE, CutPoint=0, CSI_SCALE=0.7)
ACD_DENSITY_CORRECTION <- list(
  a   = 40.6345,        # intercept (ft^2/ac)
  b   = -0.334383,      # slope on BA_t1 (per ft^2/ac)
  cap = 25,             # default +/- clamp
  crossover_BA = 121.5, # BA_t1 where predicted residual = 0
  n   = 93,             # n plots in calibration
  r2_fit = 0.1877,      # fit R^2
  r2_holdout_mean = 0.479,  # 200-iter 50/50 holdout mean test R^2
  r2_holdout_sd   = 0.115,
  fitted_on  = "2026-05-28 v30 cardinal_acadgy_residualcal_v30.R",
  applies_to = "ME FIA conditions, 10-yr remeasurement, 12.3.9 production posture"
)

apply_density_correction <- function(BA_pred, BA_t1, cap = ACD_DENSITY_CORRECTION$cap) {
  stopifnot(length(BA_pred) == length(BA_t1))
  raw_residual <- ACD_DENSITY_CORRECTION$a + ACD_DENSITY_CORRECTION$b * BA_t1
  if (!is.null(cap) && is.finite(cap) && cap > 0) {
    raw_residual <- pmax(-cap, pmin(cap, raw_residual))
  }
  BA_pred - raw_residual
}

# Convenience wrapper that also returns the raw vs corrected for diagnostics.
apply_density_correction_verbose <- function(BA_pred, BA_t1, cap = ACD_DENSITY_CORRECTION$cap) {
  raw_residual <- ACD_DENSITY_CORRECTION$a + ACD_DENSITY_CORRECTION$b * BA_t1
  capped <- if (!is.null(cap) && is.finite(cap) && cap > 0)
    pmax(-cap, pmin(cap, raw_residual)) else raw_residual
  data.frame(
    BA_t1 = BA_t1,
    BA_pred_raw = BA_pred,
    predicted_residual = raw_residual,
    capped_residual = capped,
    BA_pred_corrected = BA_pred - capped,
    was_capped = abs(raw_residual) > (cap %||% Inf)
  )
}

# Smoke test
if (interactive() || sys.nframe() == 0) {
  cat("ACD_DENSITY_CORRECTION coefficients:\n")
  for (k in names(ACD_DENSITY_CORRECTION)) cat(sprintf("  %-18s %s\n", k, ACD_DENSITY_CORRECTION[[k]]))
  cat("\nSmoke test on representative BA_t1 values (default cap=25):\n")
  test_BA_t1 <- c(20, 50, 80, 100, 121.5, 150, 180, 220)
  fake_BA_pred <- test_BA_t1 * 1.1  # pretend 10% overshoot
  d <- apply_density_correction_verbose(fake_BA_pred, test_BA_t1)
  print(format(d, digits = 4))
}
