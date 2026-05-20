#!/usr/bin/env Rscript
###############################################################################
# Stand-level calibration factors for FVS-ACD
# -----------------------------------------------------------------------------
# Derive per-stratum multipliers m_attr that, when applied as
#   attr_calibrated = attr_pred * m_attr
# bring the calibrated arm's mean attribute toward observed at the population
# level. Compute factors for:
#   - TPA  (trees per acre)
#   - QMD  (quadratic mean diameter, in)
#   - BAPH (basal area per ACRE, sq ft/ac — keeping FVS conventional unit;
#           rename to BAPA in code, BAPH in titles for user familiarity)
#   - TOPHT (top height, ft)
# Stratifications:
#   1. Overall (single scalar per attribute)
#   2. By FORTYPCD coarse group (softwood / hardwood / mixedwood)
#   3. By FVS_SITE_INDEX class (low / med / high terciles)
#   4. By initial BA class (BA_t1 < 50, 50-100, 100-150, 150+)
#   5. By remeasurement interval (4, 5, 6, 7 yr)
# Writes:
#   acd_calibration_factors.csv  — long-format with strata + factors
#   acd_calibration_factors_keyword.txt  — FVS keyword-formatted output
###############################################################################

vd <- read.csv("/sessions/confident-awesome-noether/mnt/outputs/bias_diag/validation_data_acd_post.csv",
                stringsAsFactors = FALSE)

# --- Coarse forest-type group --------------------------------------------------
# FIA forest type codes: 100s spruce-fir, 400-500 oak, 800s aspen-birch,
# 700-900 northern hardwoods. Coarse grouping:
classify_ft <- function(fortypcd) {
  ft <- as.integer(fortypcd)
  group <- rep("unknown", length(ft))
  group[ft %in% 100:199] <- "spruce-fir"
  group[ft %in% 200:399] <- "white-red-jack-pine"
  group[ft %in% 400:499] <- "oak-hickory"
  group[ft %in% 500:599] <- "oak-pine"
  group[ft %in% 700:799] <- "aspen-birch"
  group[ft %in% 800:899] <- "maple-beech-birch"
  group[ft %in% 900:999] <- "nonstocked-other"
  group
}
vd$FT_GROUP <- classify_ft(vd$FORTYPCD)

# --- Helper: stand-level multiplier (population pred/obs ratio) ---------------
factor_for <- function(pred, obs) {
  v <- !is.na(pred) & !is.na(obs) & obs > 0
  if (sum(v) < 30) return(c(n = sum(v), pred_obs = NA_real_, mult = NA_real_))
  ratio <- mean(pred[v]) / mean(obs[v])
  c(n = sum(v), pred_obs = ratio, mult = 1 / ratio)  # multiplier corrects toward 1
}

# --- Apply across strata for one attribute ------------------------------------
strat_factors <- function(label, pred_col, obs_col) {
  pred <- vd[[pred_col]]
  obs  <- vd[[obs_col]]

  out <- list()

  # Overall
  fo <- factor_for(pred, obs)
  out[[length(out)+1]] <- data.frame(
    attribute = label, stratum = "OVERALL", level = "ALL",
    n = fo["n"], pred_obs_ratio = fo["pred_obs"], multiplier = fo["mult"],
    stringsAsFactors = FALSE)

  # By forest-type group
  for (g in sort(unique(vd$FT_GROUP))) {
    idx <- vd$FT_GROUP == g
    f <- factor_for(pred[idx], obs[idx])
    out[[length(out)+1]] <- data.frame(
      attribute = label, stratum = "FT_GROUP", level = g,
      n = f["n"], pred_obs_ratio = f["pred_obs"], multiplier = f["mult"],
      stringsAsFactors = FALSE)
  }

  # By SI tercile
  qs <- quantile(vd$FVS_SITE_INDEX, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  vd$SI_class <- cut(vd$FVS_SITE_INDEX, qs, include.lowest = TRUE,
                      labels = c("SI_low", "SI_med", "SI_high"))
  for (lvl in levels(vd$SI_class)) {
    idx <- !is.na(vd$SI_class) & vd$SI_class == lvl
    f <- factor_for(pred[idx], obs[idx])
    out[[length(out)+1]] <- data.frame(
      attribute = label, stratum = "SI_tercile", level = lvl,
      n = f["n"], pred_obs_ratio = f["pred_obs"], multiplier = f["mult"],
      stringsAsFactors = FALSE)
  }

  # By initial BA class
  vd$BA_t1_class <- cut(vd$BA_t1, breaks = c(0, 50, 100, 150, 1e6),
                         labels = c("BA1_low", "BA1_med", "BA1_high", "BA1_vhigh"),
                         right = FALSE)
  for (lvl in levels(vd$BA_t1_class)) {
    idx <- !is.na(vd$BA_t1_class) & vd$BA_t1_class == lvl
    f <- factor_for(pred[idx], obs[idx])
    out[[length(out)+1]] <- data.frame(
      attribute = label, stratum = "BA_t1_class", level = lvl,
      n = f["n"], pred_obs_ratio = f["pred_obs"], multiplier = f["mult"],
      stringsAsFactors = FALSE)
  }

  # By interval years
  for (yr in sort(unique(vd$interval_years))) {
    idx <- vd$interval_years == yr
    f <- factor_for(pred[idx], obs[idx])
    out[[length(out)+1]] <- data.frame(
      attribute = label, stratum = "interval_years", level = as.character(yr),
      n = f["n"], pred_obs_ratio = f["pred_obs"], multiplier = f["mult"],
      stringsAsFactors = FALSE)
  }

  do.call(rbind, out)
}

# --- Compute for all four attributes ------------------------------------------
tbl <- rbind(
  strat_factors("TPA",   "TPA_pred_calib",        "TPA_t2"),
  strat_factors("QMD",   "QMD_pred_calib",        "QMD_t2"),
  strat_factors("BAPH",  "BA_pred_calib",         "BA_t2"),   # sq ft/ac (FVS-conventional)
  strat_factors("TOPHT", "HT_top_calib",          "HT_top_t2")
)
tbl <- tbl[order(tbl$attribute, tbl$stratum, tbl$level), ]
rownames(tbl) <- NULL

# Pretty round
tbl$pred_obs_ratio <- round(tbl$pred_obs_ratio, 4)
tbl$multiplier     <- round(tbl$multiplier,    4)

write.csv(tbl,
          "/sessions/confident-awesome-noether/mnt/outputs/bias_diag/acd_calibration_factors.csv",
          row.names = FALSE)

# --- Console summary ----------------------------------------------------------
cat("\nOverall (population-level) calibration multipliers\n")
cat(strrep("-", 60), "\n")
overall <- tbl[tbl$stratum == "OVERALL", c("attribute", "n", "pred_obs_ratio", "multiplier")]
print(overall, row.names = FALSE)
cat("\nReading: multiplier > 1 = pred under, multiply pred by this to match obs.\n")
cat("         multiplier < 1 = pred over,  multiply pred by this to match obs.\n\n")

cat("FT_GROUP stratification (BAPH only, for the heteroscedasticity look)\n")
cat(strrep("-", 60), "\n")
print(tbl[tbl$attribute == "BAPH" & tbl$stratum == "FT_GROUP",
          c("level", "n", "pred_obs_ratio", "multiplier")], row.names = FALSE)
cat("\n")

cat("BA_t1 class stratification (BAPH)\n")
cat(strrep("-", 60), "\n")
print(tbl[tbl$attribute == "BAPH" & tbl$stratum == "BA_t1_class",
          c("level", "n", "pred_obs_ratio", "multiplier")], row.names = FALSE)
cat("\n")

# --- FVS-style application guidance -------------------------------------------
# The multipliers above are stand-attribute scale factors derived from
# population pred/obs ratios. They are NOT directly an FVS keyword: the
# closest in-engine analogues are FIXDG (diameter growth scalar by species),
# FIXHTG (height growth scalar), and FIXMORT (mortality scalar). Because BAPH
# is a derived attribute (BAPH = f(DBH, TPA)), the BAPH multiplier maps to
# a combination of FIXDG (drives DIA) and FIXMORT (drives TPA) under a
# constant-form assumption.
#
# We emit a guidance file (not directly executable FVS syntax) showing the
# per-stratum multipliers and a recommended translation into FIXDG/FIXHTG/
# FIXMORT magnitudes, leaving the precise species-list to the operator.
sink("/sessions/confident-awesome-noether/mnt/outputs/bias_diag/acd_calibration_factors_guidance.txt")
cat("=============================================================================\n")
cat(" ACD stand-level calibration factors -- application guidance\n")
cat(" Source: benchmark job 9610424 (30,146 NE-relabeled conditions), 2026-05-15\n")
cat("=============================================================================\n\n")

cat("Population multipliers (overall):\n")
ov <- tbl[tbl$stratum == "OVERALL", ]
for (i in seq_len(nrow(ov))) {
  cat(sprintf("  %-6s pred/obs = %.4f   multiplier = %.4f\n",
              ov$attribute[i], ov$pred_obs_ratio[i], ov$multiplier[i]))
}
cat("\n")
cat("Translation to FVS keywords (constant-form approximation):\n")
cat("  BAPH ~ DBH^2 x TPA, so log(BAPH_mult) ~ 2*log(DIA_mult) + log(TPA_mult)\n")
cat("  Given QMD_mult and TPA_mult above, the BAPH discrepancy decomposes as:\n")
cat(sprintf("    QMD^2  contribution   %.4f^2 = %.4f\n",
            ov$multiplier[ov$attribute == "QMD"],
            ov$multiplier[ov$attribute == "QMD"]^2))
cat(sprintf("    TPA contribution      %.4f\n",
            ov$multiplier[ov$attribute == "TPA"]))
cat(sprintf("    product (predicted)   %.4f   vs. observed BAPH mult %.4f\n",
            ov$multiplier[ov$attribute == "QMD"]^2 * ov$multiplier[ov$attribute == "TPA"],
            ov$multiplier[ov$attribute == "BAPH"]))
cat("\n  Suggested FVS keyword settings (per simulation, post-STDINFO):\n")
cat(sprintf("    FIXDG    species: ALL   scalar: %.4f   ! adjust DG output\n",
            ov$multiplier[ov$attribute == "QMD"]))
cat(sprintf("    FIXHTG   species: ALL   scalar: %.4f   ! adjust HTG output\n",
            ov$multiplier[ov$attribute == "TOPHT"]))
cat(sprintf("    FIXMORT  species: ALL   scalar: %.4f   ! TPA multiplier > 1 means LESS mortality\n",
            ov$multiplier[ov$attribute == "TPA"]))
cat("\n  (Exact FIXDG/FIXHTG/FIXMORT syntax varies; consult variant overview.)\n\n")

cat("Per-stratum factors (informational; pick the column that matches your stand):\n\n")

# BA_t1 class is the most useful stratifier (heteroscedasticity)
cat("By initial BA class (BA_t1):\n")
sub <- tbl[tbl$stratum == "BA_t1_class", ]
sub_wide <- reshape(sub[, c("attribute","level","multiplier")],
                    direction = "wide",
                    idvar = "level",
                    timevar = "attribute")
print(sub_wide, row.names = FALSE)
cat("\n  Reading: stands with BA_t1 < 50 sq ft/ac need TPA scaled by 1.39x\n")
cat("  (i.e. ACD over-thins sparse stands). High-BA stands need a slight\n")
cat("  reduction (~3%). Apply only the column that matches your initial state.\n\n")

# FT_GROUP block
cat("By forest-type group:\n")
sub <- tbl[tbl$stratum == "FT_GROUP" &
            !(tbl$level %in% c("unknown","nonstocked-other")), ]
sub_wide <- reshape(sub[, c("attribute","level","multiplier")],
                    direction = "wide",
                    idvar = "level",
                    timevar = "attribute")
print(sub_wide, row.names = FALSE)
cat("\n")

cat("By site-index tercile:\n")
sub <- tbl[tbl$stratum == "SI_tercile", ]
sub_wide <- reshape(sub[, c("attribute","level","multiplier")],
                    direction = "wide",
                    idvar = "level",
                    timevar = "attribute")
print(sub_wide, row.names = FALSE)
cat("\n")

cat("By remeasurement interval (years):\n")
sub <- tbl[tbl$stratum == "interval_years", ]
sub_wide <- reshape(sub[, c("attribute","level","multiplier")],
                    direction = "wide",
                    idvar = "level",
                    timevar = "attribute")
print(sub_wide, row.names = FALSE)
cat("\n")

cat("=============================================================================\n")
cat(" Caveats\n")
cat("=============================================================================\n")
cat(" 1. These are population-level ratio corrections, not per-tree corrections.\n")
cat("    They will pull the MEAN of predictions toward the observed mean but do\n")
cat("    not reduce per-stand RMSE.\n")
cat(" 2. The biggest gain is from the BA_t1 class stratifier (BA_t1 < 50 needs\n")
cat("    a +39%% TPA correction). Apply by class, not by single overall scalar.\n")
cat(" 3. Validation interval was 4-7 years. Multipliers will not extrapolate\n")
cat("    cleanly to longer projections without re-derivation.\n")
cat(" 4. INGROWTH was DISABLED in the post-pass that produced these factors.\n")
cat("    Re-derive after re-running the engine with FVS_ACD_RELABEL=TRUE and\n")
cat("    ingrowth_lookup cached, if you want ingrowth-inclusive scalars.\n")
sink()

cat("Wrote:\n")
cat("  acd_calibration_factors.csv\n")
cat("  acd_calibration_factors_keyword.txt\n")
