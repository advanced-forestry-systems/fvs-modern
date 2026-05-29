#!/usr/bin/env Rscript
# silc_cfi_pair_list.R
# =====================================================================
# Build the 24 prediction tasks for the SILC CFI four-model benchmark.
# For each STAND_PAI record, capture:
#   * plot, year_prev, year_curr, period_yr
#   * observed stand-level metrics at year_prev (from STAND_METRICS)
#   * observed stand-level metrics at year_curr (from STAND_METRICS)
#   * observed PAI components from STAND_PAI
#   * year-T_prev tree list as a CSV per pair, ready for AcadianGY
#
# Output:
#   silc_cfi_pair_summary.csv   one row per pair, used by the scorecard
#   pair_input/pair_<P>_<Yprev>_<Ycurr>_tree.csv  per-pair tree list
#   silc_cfi_naive_scorecard.csv  zero-growth + linear-PAI baselines
# =====================================================================
od  <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
in_dir <- file.path(od, "pair_input")
dir.create(in_dir, showWarnings = FALSE)

sm  <- read.csv(file.path(od, "STAND_METRICS.csv"))
pai <- read.csv(file.path(od, "STAND_PAI.csv"))
tr  <- read.csv(file.path(od, "TREE.csv"))

pair_ok <- pai[!is.na(pai$PAI_BA_NET_FT2ACY), ]
cat(sprintf("Building pair list: %d records\n", nrow(pair_ok)))

EXPF <- 5.0  # plot expansion factor, 1/5 acre fixed plot

# Tree-list export per pair
emit_tree_list <- function(plot_id, year_prev) {
  t <- tr[tr$PLOT == plot_id & tr$MEASYEAR == year_prev &
          tr$STATUSCD == 1 & is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  if (nrow(t) == 0) return(0)
  out <- data.frame(
    STAND = sprintf("CFI_%04d", plot_id),
    YEAR  = year_prev,
    PLOT  = plot_id,
    TREE  = t$TREE,
    SP    = t$COMMON_NAME,
    SPCD  = t$SPCD,
    DBH   = t$DIA_IN,
    HT    = t$HT_FT,
    HCB   = NA,
    EXPF  = EXPF,
    Form  = 1,
    Risk  = 1
  )
  f <- file.path(in_dir,
                 sprintf("pair_%04d_%d_tree.csv", plot_id, year_prev))
  write.csv(out, f, row.names = FALSE)
  nrow(out)
}

# Build summary
summ <- pair_ok[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR",
                    "BA_PREV_FT2AC","BA_CURR_FT2AC",
                    "PAI_BA_SURV_FT2ACY","PAI_BA_INGR_FT2ACY",
                    "PAI_BA_MORT_FT2ACY","PAI_BA_NET_FT2ACY",
                    "MEAN_ANN_DIA_IN_YR","N_GROWTH_RECS","N_MORT_RECS")]
# Add TPA / QMD prev and curr from STAND_METRICS
get_sm <- function(p, y, col) {
  v <- sm[sm$PLOT==p & sm$MEASYEAR==y, col]
  if (length(v) == 0) NA else v
}
summ$TPA_PREV  <- mapply(get_sm, summ$PLOT, summ$YEAR_PREV, "TPA")
summ$TPA_CURR  <- mapply(get_sm, summ$PLOT, summ$YEAR_CURR, "TPA")
summ$QMD_PREV  <- mapply(get_sm, summ$PLOT, summ$YEAR_PREV, "QMD_IN")
summ$QMD_CURR  <- mapply(get_sm, summ$PLOT, summ$YEAR_CURR, "QMD_IN")
summ$N_TREES_PREV_LIVE <- mapply(function(p, y) {
  sum(tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 & tr$DIA_IN >= 4.5)
}, summ$PLOT, summ$YEAR_PREV)

# Write per-pair tree lists
summ$tree_list_file <- mapply(function(p, y) {
  n <- emit_tree_list(p, y)
  if (n == 0) "" else
    sprintf("pair_input/pair_%04d_%d_tree.csv", p, y)
}, summ$PLOT, summ$YEAR_PREV)
summ$tree_count_in_list <- mapply(function(p, y) {
  sum(tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 & tr$DIA_IN >= 4.5)
}, summ$PLOT, summ$YEAR_PREV)

# Baseline predictors:
# (1) zero-growth: BA_pred = BA_prev
# (2) linear PAI: assumes observed net PAI continues at same rate; trivial
# (3) FIA mean net PAI: use the FVS-NE FIA-fit average net BA growth rate
#     for spruce-fir as a generic prior. We use 0.93 ft^2/ac/yr (Westfall+CFRU)
#     for spruce-fir / mixedwood; 0.65 for hardwood; 0.40 for cedar.

# Map plot -> forest type from the strata mapping
strata_map <- read.csv(file.path(od, "silc_cfi_plot_strata_map.csv"))
ft_lookup <- setNames(strata_map$forest_type, strata_map$PLOT)
summ$forest_type <- ft_lookup[as.character(summ$PLOT)]

prior_PAI <- c("Cedar"               = 0.40,
               "Hardwood"            = 0.65,
               "Mixedwood"           = 0.85,
               "Commercial Softwood" = 0.93,
               "Other Softwood"      = 0.70,
               "Unclassifiable"      = NA_real_)
summ$prior_PAI_BA_ft2ac_yr <- prior_PAI[summ$forest_type]
summ$BA_pred_zero_growth   <- summ$BA_PREV_FT2AC
summ$BA_pred_FIA_prior     <- summ$BA_PREV_FT2AC +
                                summ$prior_PAI_BA_ft2ac_yr *
                                summ$PERIOD_YR

bias <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)
rmse <- function(p, o) sqrt(mean((p - o)^2, na.rm=TRUE))

cat("\n=== Naive baseline scorecard (BA at year_curr, ft^2/ac) ===\n")
sb <- data.frame(
  predictor = c("zero growth",  "FIA prior PAI"),
  bias_pct  = c(bias(summ$BA_pred_zero_growth, summ$BA_CURR_FT2AC),
                bias(summ$BA_pred_FIA_prior,   summ$BA_CURR_FT2AC)),
  RMSE      = c(rmse(summ$BA_pred_zero_growth, summ$BA_CURR_FT2AC),
                rmse(summ$BA_pred_FIA_prior,   summ$BA_CURR_FT2AC)),
  n = sum(is.finite(summ$BA_CURR_FT2AC))
)
print(sb, row.names=FALSE, digits=3)
write.csv(sb, file.path(od, "silc_cfi_naive_scorecard.csv"), row.names=FALSE)

write.csv(summ, file.path(od, "silc_cfi_pair_summary.csv"),
          row.names = FALSE)
cat(sprintf("\nWrote pair summary (%d pairs) and %d per-pair tree-list CSVs.\n",
            nrow(summ), sum(summ$tree_count_in_list > 0)))
