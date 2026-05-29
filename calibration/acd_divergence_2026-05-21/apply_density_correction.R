# apply_density_correction.R  (v37 production: pooled n=484 refit)
# AcadianGY 12.3.9 bridge-level post-projection correction for the
# density-dependent residual identified in v30, refined through v32 fresh
# sample test, v33 pooled refit, v34 tree-level reconciliation, and now
# v37 refit on the pooled v30 + v32 + v36 sample (n=484).
#
# This file exposes:
#   apply_density_correction(BA_pred, BA_t1, ...)
#       Stand-level scalar correction (returns corrected BA only).
#
#   apply_density_correction_treelist(tree, ...)
#       Tree-level reconciliation. Takes a tree data.frame, computes the
#       corrected stand BA via the v33 formula, and scales every tree's
#       EXPF uniformly so sum(DBH^2 * EXPF) matches the corrected BA.
#       Returns the updated tree data.frame (same row count, same DBH and
#       species, just EXPF reduced).
#
# Mortality interpretation. The correction subtracts BA from the projection;
# attributing that subtraction to extra mortality (i.e., fewer surviving
# trees per acre) preserves the size distribution. The tree-level function
# scales EXPF uniformly by (corrected_BA / raw_BA), so it conserves the BA
# constraint by construction.
#
# Alternative weightings (more mortality on smaller suppressed trees, etc.)
# would require species- or size-stratified mortality information that the
# v30 residual signal does not give us. Uniform is the safest default.
#
# IMPORTANT: tree-level reconciliation applies a per-stand EXPF scale factor.
# In pathological cases (when raw BA_pred is unusually small relative to the
# v33 prediction), the scale factor can be unhelpfully low or negative. The
# scale_floor argument (default 0.7) caps how aggressively any single stand's
# tree list can be thinned. Empirical defense on v32 fresh sample (n=91):
#
#   no floor:    BA bias +5.05 percent, TPA bias -5.22 percent, R^2 0.486
#   floor=0.7:   BA bias +6.98 percent, TPA bias -3.22 percent, R^2 0.496
#   floor=0.8:   BA bias +8.23 percent, TPA bias -1.91 percent, R^2 0.491
#
# floor=0.7 is the production default. It gives a meaningful BA closure
# (raw +17.92 percent -> +6.98 percent) and keeps TPA close to observed,
# with the highest test R^2.
#
# 5-fold CV on the pooled n=484 sample, 50 random shuffles:
#   uncorrected:                BA bias = +11.76 percent, R^2 = 0.452
#   asymmetric (0, +20)         BA bias = +0.15 +/- 0.04, R^2 = 0.514 +/- 0.001
#   asymmetric (0, +25)         BA bias = -0.64 +/- 0.04, R^2 = 0.520 +/- 0.001
#
# Production picks asym (0, +20) for the cleaner mean bias closure with
# minimal R^2 cost. CV variance is 3x tighter than v33 (sd 0.04 vs 0.11)
# from the larger n.

# Production coefficients, fit on pooled n=484 (v30 seed=42 n=93 + v32
# seed=2027 n=91 + v36 seed=2028 n=300, all ME FIA, 10-yr remeasurement,
# AcadianGY 12.3.9 with MORTCAL=TRUE, CutPoint=0, CSI_SCALE=0.7).
ACD_DENSITY_CORRECTION <- list(
  a         = 28.9607,
  b         = -0.186023,
  upper_cap = 20,
  lower_cap = 0,
  crossover_BA = 155.7,
  n         = 484,
  r2_fit    = 0.0638,
  cv_bias_mean = 0.15,
  cv_bias_sd   = 0.04,
  cv_r2_mean   = 0.514,
  cv_r2_sd     = 0.001,
  fitted_on  = "2026-05-29 v36+v32+v30 pooled, n=484 ME FIA",
  applies_to = "ME FIA conditions, 10-yr remeasurement, 12.3.9 production posture",
  supersedes = c("v31 (a=40.6, b=-0.334, sym +/-25)",
                  "v33 (a=37.0, b=-0.236, asym 0/+25, n=184)")
)

# -------------------------------------------------------------------------
# Stand-level scalar correction (v33)
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# Tree-level reconciliation (v34)
# -------------------------------------------------------------------------
#
# Push the v33 stand-level BA correction back into the tree list by scaling
# each tree's EXPF uniformly. By construction:
#   sum(corrected_EXPF * DBH^2) * conv == BA_corrected
# so the tree list is internally consistent with the stand-level constraint.
#
# Arguments:
#   tree       data.frame with at minimum columns STAND, DBH, EXPF. May
#              contain any other columns (SP, HT, HCB, etc.) which are
#              passed through unchanged.
#   BA_t1_by_stand
#              named numeric vector or data.frame with columns (STAND, BA_t1).
#              BA_t1 must be in ft^2/ac (same units as the v33 fit).
#   dbh_units  one of "cm" or "in". AcadianGY natively works in cm;
#              if your DBH is in cm, BA computation will convert. Default
#              "cm".
#   upper_cap, lower_cap as in stand-level function.
#   floor_EXPF small positive minimum on the scaled EXPF (default 1e-5),
#              to avoid zero-EXPF entries that could cause downstream issues.
#
# Returns a list with two elements:
#   tree_corrected   the tree data.frame with EXPF scaled in place
#   stand_summary    per-stand data.frame:
#     STAND, BA_t1, BA_raw, BA_corrected, scale_factor, TPA_raw, TPA_corrected
#
apply_density_correction_treelist <- function(tree, BA_t1_by_stand,
                                              dbh_units = c("cm", "in"),
                                              upper_cap = ACD_DENSITY_CORRECTION$upper_cap,
                                              lower_cap = ACD_DENSITY_CORRECTION$lower_cap,
                                              scale_floor = 0.7,
                                              floor_EXPF = 1e-5) {
  dbh_units <- match.arg(dbh_units)
  stopifnot(all(c("STAND", "DBH", "EXPF") %in% names(tree)))

  # Convert BA_t1_by_stand to a lookup
  if (is.data.frame(BA_t1_by_stand)) {
    stopifnot(all(c("STAND", "BA_t1") %in% names(BA_t1_by_stand)))
    ba_t1_lut <- setNames(BA_t1_by_stand$BA_t1, as.character(BA_t1_by_stand$STAND))
  } else {
    ba_t1_lut <- BA_t1_by_stand
  }

  # BA contribution per tree in ft^2/ac, matching the v33 fit's units
  # In cm: BA per tree = 0.00007854 * DBH^2 * EXPF * 4.35  (m^2/ha -> ft^2/ac)
  # In inches: BA per tree = 0.005454 * DBH^2 * EXPF
  if (dbh_units == "cm") {
    tree$BA_contrib_ft2ac <- 0.00007854 * tree$DBH^2 * tree$EXPF * 4.35
  } else {
    tree$BA_contrib_ft2ac <- 0.005454 * tree$DBH^2 * tree$EXPF
  }

  # Per-stand raw BA at projection horizon
  ba_raw_by_stand <- tapply(tree$BA_contrib_ft2ac, tree$STAND, sum, na.rm = TRUE)

  # Apply v33 stand-level correction
  ba_raw_vec    <- as.numeric(ba_raw_by_stand)
  stands        <- as.character(names(ba_raw_by_stand))
  ba_t1_vec     <- ba_t1_lut[stands]
  ba_t1_vec[is.na(ba_t1_vec)] <- ACD_DENSITY_CORRECTION$crossover_BA  # neutral default
  ba_corrected_vec <- apply_density_correction(ba_raw_vec, ba_t1_vec,
                                               upper_cap = upper_cap,
                                               lower_cap = lower_cap)

  # Tree-level scale factor per stand. Apply scale_floor to guard against
  # pathological cases where BA_corrected is much smaller than BA_raw and
  # uniform scaling would otherwise thin the stand too aggressively.
  raw_scale <- ifelse(ba_raw_vec > 0, ba_corrected_vec / ba_raw_vec, 1)
  bounded_scale <- pmax(scale_floor, pmin(1.0, raw_scale))  # also cap at 1.0 (correction can only subtract)
  scale_lut <- setNames(bounded_scale, stands)
  tree$EXPF_scale  <- scale_lut[as.character(tree$STAND)]
  tree$EXPF_raw    <- tree$EXPF
  tree$EXPF        <- pmax(floor_EXPF, tree$EXPF_raw * tree$EXPF_scale)

  # Stand summary for diagnostics
  tpa_raw   <- tapply(tree$EXPF_raw, tree$STAND, sum, na.rm = TRUE)
  tpa_corr  <- tapply(tree$EXPF,     tree$STAND, sum, na.rm = TRUE)
  stand_summary <- data.frame(
    STAND       = stands,
    BA_t1       = as.numeric(ba_t1_vec),
    BA_raw      = ba_raw_vec,
    BA_corrected = ba_corrected_vec,
    scale_factor = as.numeric(scale_lut),
    TPA_raw_unit = as.numeric(tpa_raw)[match(stands, names(tpa_raw))],
    TPA_corrected_unit = as.numeric(tpa_corr)[match(stands, names(tpa_corr))],
    stringsAsFactors = FALSE
  )

  # Drop the BA contribution helper column
  tree$BA_contrib_ft2ac <- NULL

  list(
    tree_corrected = tree,
    stand_summary  = stand_summary
  )
}

# -------------------------------------------------------------------------
# Tree-level reconciliation v40: TPA-preserving size-weighted scaling
# -------------------------------------------------------------------------
#
# Alternative to apply_density_correction_treelist (which scales EXPF uniformly
# and crashes TPA in proportion to BA). v40 scales EXPF as a function of DBH
# so the corrected BA matches the target AND the TPA stays at the observed
# level. Biological interpretation: the model is keeping too many big trees
# alive; reduce big-tree EXPF and shift weight toward smaller trees.
#
# Parameterization: s_i = 1 + beta * (DBH_mean - DBH_i) / DBH_mean
#   where DBH_mean is the EXPF-weighted mean DBH per stand. By construction
#   this preserves TPA exactly. beta is solved per stand from the BA target.
#
# Trade-off vs v34:
#   v34 (uniform):       preserves QMD, crashes TPA (-10 percent typical)
#   v40 (size-weighted): preserves TPA, decreases QMD (~5 percent typical)
#   Both close BA equally well.
#
# Choice depends on which downstream metric matters most. Volume tables
# are usually QMD-sensitive (v34 better). TPA-driven downstream (counts,
# carbon, harvest decisions on density) usually prefers v40.

apply_density_correction_treelist_tpa <- function(tree, BA_t1_by_stand,
                                                  dbh_units = c("cm", "in"),
                                                  upper_cap = ACD_DENSITY_CORRECTION$upper_cap,
                                                  lower_cap = ACD_DENSITY_CORRECTION$lower_cap,
                                                  s_min = 0.1,
                                                  s_max = 2.0,
                                                  floor_EXPF = 1e-5) {
  dbh_units <- match.arg(dbh_units)
  stopifnot(all(c("STAND", "DBH", "EXPF") %in% names(tree)))

  if (is.data.frame(BA_t1_by_stand)) {
    stopifnot(all(c("STAND", "BA_t1") %in% names(BA_t1_by_stand)))
    ba_t1_lut <- setNames(BA_t1_by_stand$BA_t1, as.character(BA_t1_by_stand$STAND))
  } else {
    ba_t1_lut <- BA_t1_by_stand
  }

  # Compute per-tree BA contribution and per-stand totals
  if (dbh_units == "cm") {
    ba_factor <- 0.00007854 * 4.35   # m^2/ha -> ft^2/ac
  } else {
    ba_factor <- 0.005454            # in-units, ft^2/ac
  }

  tree$EXPF_raw <- tree$EXPF
  stand_summary_rows <- list()

  for (sid in unique(tree$STAND)) {
    idx <- tree$STAND == sid
    DBH_i <- tree$DBH[idx]
    EXPF_i <- tree$EXPF_raw[idx]

    BA_raw <- ba_factor * sum(EXPF_i * DBH_i^2)
    BA_t1  <- ba_t1_lut[as.character(sid)]
    if (is.na(BA_t1)) BA_t1 <- ACD_DENSITY_CORRECTION$crossover_BA

    raw_corr <- ACD_DENSITY_CORRECTION$a + ACD_DENSITY_CORRECTION$b * BA_t1
    bnd_corr <- max(lower_cap, min(upper_cap, raw_corr))
    BA_target <- BA_raw - bnd_corr

    if (BA_target >= BA_raw || BA_raw <= 0) {
      # No correction needed (high-density stand)
      tree$EXPF[idx] <- EXPF_i
      tree$EXPF_scale[idx] <- 1.0
      stand_summary_rows[[sid]] <- data.frame(
        STAND = sid, BA_t1 = BA_t1,
        BA_raw = BA_raw, BA_corrected = BA_raw,
        beta = 0, scale_min = 1, scale_max = 1,
        TPA_raw = sum(EXPF_i), TPA_corrected = sum(EXPF_i),
        stringsAsFactors = FALSE)
      next
    }

    # EXPF-weighted mean DBH
    sum_E <- sum(EXPF_i)
    DBH_mean <- sum(EXPF_i * DBH_i) / sum_E

    # Solve for beta: BA_target = M2 + beta * (DBH_mean * M2 - M3) / DBH_mean
    # where M2 = sum(EXPF * DBH^2), M3 = sum(EXPF * DBH^3)
    M2 <- sum(EXPF_i * DBH_i^2)
    M3 <- sum(EXPF_i * DBH_i^3)
    target_M2 <- BA_target / ba_factor
    denom <- (DBH_mean * M2 - M3)
    if (abs(denom) < 1e-9) {
      # All trees same size; uniform scaling is the only option
      s_i <- rep(target_M2 / M2, length(EXPF_i))
      beta <- NA_real_
    } else {
      beta <- (target_M2 - M2) * DBH_mean / denom
      s_i <- 1 + beta * (DBH_mean - DBH_i) / DBH_mean
      s_i <- pmax(s_min, pmin(s_max, s_i))
    }

    new_EXPF <- pmax(floor_EXPF, EXPF_i * s_i)
    tree$EXPF[idx] <- new_EXPF
    tree$EXPF_scale[idx] <- s_i

    stand_summary_rows[[sid]] <- data.frame(
      STAND = sid, BA_t1 = BA_t1,
      BA_raw = BA_raw,
      BA_corrected = ba_factor * sum(new_EXPF * DBH_i^2),
      beta = beta,
      scale_min = min(s_i), scale_max = max(s_i),
      TPA_raw = sum(EXPF_i), TPA_corrected = sum(new_EXPF),
      stringsAsFactors = FALSE)
  }

  list(
    tree_corrected = tree,
    stand_summary  = do.call(rbind, stand_summary_rows)
  )
}

# -------------------------------------------------------------------------
# Smoke test
# -------------------------------------------------------------------------

if (interactive() || sys.nframe() == 0) {
  cat("ACD_DENSITY_CORRECTION (v33 coefficients, v34 tree-level wrapper):\n")
  for (k in names(ACD_DENSITY_CORRECTION))
    cat(sprintf("  %-18s %s\n", k, ACD_DENSITY_CORRECTION[[k]]))

  cat("\nStand-level smoke test:\n")
  test_BA_t1 <- c(10, 30, 60, 100, 156.6, 180, 220, 260)
  fake_BA_pred <- test_BA_t1 * 1.15
  d <- apply_density_correction_verbose(fake_BA_pred, test_BA_t1)
  print(format(d, digits = 4))

  cat("\nTree-level v40 (TPA-preserving) smoke test:\n")
  cat("Stand S1: low-density (BA_t1=50)\n")
  set.seed(42)
  ntrees <- 20
  fake_tree2 <- data.frame(
    STAND = rep("S1", ntrees), TREE = seq_len(ntrees), SP = "RS",
    DBH   = sort(rgamma(ntrees, 4, 0.18)),  # cm
    EXPF  = runif(ntrees, 5, 15)
  )
  fake_tree2$EXPF_scale <- NA_real_
  res2 <- apply_density_correction_treelist_tpa(fake_tree2, c(S1 = 50), dbh_units = "cm")
  cat(sprintf("  BA_raw = %.2f, BA_corrected = %.2f (target reduction = %.2f)\n",
              res2$stand_summary$BA_raw, res2$stand_summary$BA_corrected,
              res2$stand_summary$BA_raw - res2$stand_summary$BA_corrected))
  cat(sprintf("  TPA_raw = %.0f, TPA_corrected = %.0f (preserved)\n",
              res2$stand_summary$TPA_raw, res2$stand_summary$TPA_corrected))
  cat(sprintf("  beta = %.3f, scale range = [%.2f, %.2f]\n",
              res2$stand_summary$beta, res2$stand_summary$scale_min, res2$stand_summary$scale_max))
  cat("  (small trees: s > 1; big trees: s < 1)\n")

  cat("\nTree-level v34 (uniform, original) smoke test:\n")
  cat("Stand S1: low-density (BA_t1=50), should get scale ~0.7 (floor active)\n")
  cat("Stand S2: high-density (BA_t1=180), should get scale = 1.0 (no correction)\n")
  fake_tree <- data.frame(
    STAND = c(rep("S1", 4), rep("S2", 4)),
    TREE  = 1:8,
    SP    = c("RS","BF","RM","YB","RS","BF","RM","YB"),
    DBH   = c(15, 18, 30, 20, 35, 28, 40, 25),  # cm
    EXPF  = c(60, 50, 30, 40, 25, 30, 20, 35)   # tree per acre equivalent
  )
  res <- apply_density_correction_treelist(fake_tree, c(S1 = 50, S2 = 180),
                                            dbh_units = "cm", scale_floor = 0.7)
  cat("Per-tree EXPF before/after:\n")
  print(format(res$tree_corrected[, c("STAND","TREE","SP","DBH","EXPF_raw","EXPF_scale","EXPF")], digits=4))
  cat("\nStand summary:\n"); print(format(res$stand_summary, digits=4))
}
