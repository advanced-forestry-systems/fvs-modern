#!/usr/bin/env Rscript
# silc_cfi_longhorizon_pairs.R
# =====================================================================
# Build the long-horizon prediction pair list: for each CFI plot, pair
# its earliest reliable measurement with its latest. This is the
# operationally relevant test for SILC since they project 10-50 years
# out, not 5.
#
# Output:
#   silc_cfi_longhorizon_pairs.csv  one row per plot (10 rows expected)
#   pair_input_long/pair_<PLOT>_<Yprev>_tree.csv  tree list at year_prev
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
in_dir <- file.path(od, "pair_input_long")
dir.create(in_dir, showWarnings = FALSE)

sm <- read.csv(file.path(od, "STAND_METRICS.csv"))
tr <- read.csv(file.path(od, "TREE.csv"))
EXPF <- 5.0

# Reliable plot-years (METRICS_RELIABLE == "Y") with BA > 10 to avoid
# establishment-only intervals
sm_ok <- sm[sm$METRICS_RELIABLE == "Y" & sm$BA_FT2_AC >= 10, ]

# Per plot: earliest year, latest year
plot_span <- aggregate(MEASYEAR ~ PLOT, data = sm_ok,
                       FUN = function(x) c(min = min(x), max = max(x)))
plot_span <- do.call(data.frame,
                     list(PLOT = plot_span$PLOT,
                          YEAR_PREV = plot_span$MEASYEAR[, "min"],
                          YEAR_CURR = plot_span$MEASYEAR[, "max"]))
plot_span$PERIOD_YR <- plot_span$YEAR_CURR - plot_span$YEAR_PREV
plot_span <- plot_span[plot_span$PERIOD_YR >= 5, ]

# Pull stand-level metrics at year_prev and year_curr
sm_lookup <- function(p, y, col) {
  v <- sm[sm$PLOT == p & sm$MEASYEAR == y, col]
  if (length(v) == 0) NA else v
}
plot_span$BA_PREV_FT2AC <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_PREV, "BA_FT2_AC")
plot_span$BA_CURR_FT2AC <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "BA_FT2_AC")
plot_span$TPA_PREV     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_PREV, "TPA")
plot_span$TPA_CURR     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "TPA")
plot_span$QMD_PREV     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_PREV, "QMD_IN")
plot_span$QMD_CURR     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "QMD_IN")
plot_span$SDI_obs_curr <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "SDI")
plot_span$CurtisRD_obs_curr <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "CURTIS_RD")

# Per-pair tree list
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
    Form  = 1, Risk = 1
  )
  f <- file.path(in_dir,
                 sprintf("pair_%04d_%d_tree.csv", plot_id, year_prev))
  write.csv(out, f, row.names = FALSE)
  nrow(out)
}
plot_span$tree_list_file <- mapply(function(p, y) {
  n <- emit_tree_list(p, y)
  if (n == 0) "" else
    sprintf("pair_input_long/pair_%04d_%d_tree.csv", p, y)
}, plot_span$PLOT, plot_span$YEAR_PREV)
plot_span$tree_count_prev <- mapply(function(p, y) {
  sum(tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 &
      is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5)
}, plot_span$PLOT, plot_span$YEAR_PREV)

write.csv(plot_span,
          file.path(od, "silc_cfi_longhorizon_pairs.csv"),
          row.names = FALSE)
cat("=== Long-horizon CFI pairs ===\n")
print(plot_span[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR",
                     "BA_PREV_FT2AC","BA_CURR_FT2AC","tree_count_prev")],
      row.names = FALSE)
cat(sprintf("\nMean horizon: %.1f years\n", mean(plot_span$PERIOD_YR)))
cat(sprintf("Wrote %d pair tree-list files to %s/\n",
            sum(plot_span$tree_count_prev > 0), in_dir))
