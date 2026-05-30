#!/usr/bin/env Rscript
# silc_cfi_longhorizon_pairs_v16.R
# v16 of SILC CFI: 154 plots (was 10 in v3). Build long horizon pairs
# (earliest reliable to latest reliable per plot) and emit tree lists.
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
in_dir <- file.path(od, "pair_input_long_v16")
dir.create(in_dir, showWarnings = FALSE)

sm <- read.csv(file.path(od, "v16/STAND_METRICS.csv"), stringsAsFactors = FALSE)
tr <- read.csv(file.path(od, "v16/TREE.csv"), stringsAsFactors = FALSE)
EXPF <- 5.0

sm_ok <- sm[!is.na(sm$METRICS_RELIABLE) & sm$METRICS_RELIABLE == "Y" &
            !is.na(sm$BA_FT2_AC) & sm$BA_FT2_AC >= 10, ]

plot_span <- aggregate(MEASYEAR ~ PLOT, data = sm_ok,
                       FUN = function(x) c(min = min(x), max = max(x)))
plot_span <- do.call(data.frame,
                     list(PLOT = plot_span$PLOT,
                          YEAR_PREV = plot_span$MEASYEAR[, "min"],
                          YEAR_CURR = plot_span$MEASYEAR[, "max"]))
plot_span$PERIOD_YR <- plot_span$YEAR_CURR - plot_span$YEAR_PREV
plot_span <- plot_span[plot_span$PERIOD_YR >= 5, ]

sm_lookup <- function(p, y, col) {
  v <- sm[sm$PLOT == p & sm$MEASYEAR == y, col]
  if (length(v) == 0) NA else v[1]
}
plot_span$BA_PREV_FT2AC <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_PREV, "BA_FT2_AC")
plot_span$BA_CURR_FT2AC <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "BA_FT2_AC")
plot_span$TPA_PREV     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_PREV, "TPA")
plot_span$TPA_CURR     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "TPA")
plot_span$QMD_PREV     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_PREV, "QMD_IN")
plot_span$QMD_CURR     <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "QMD_IN")
plot_span$SDI_obs_curr <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "SDI")
plot_span$CurtisRD_obs_curr <- mapply(sm_lookup, plot_span$PLOT, plot_span$YEAR_CURR, "CURTIS_RD")

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
  f <- file.path(in_dir, sprintf("pair_%04d_%d_tree.csv", plot_id, year_prev))
  write.csv(out, f, row.names = FALSE)
  nrow(out)
}
plot_span$tree_list_file <- mapply(function(p, y) {
  n <- emit_tree_list(p, y)
  if (n == 0) "" else
    sprintf("pair_input_long_v16/pair_%04d_%d_tree.csv", p, y)
}, plot_span$PLOT, plot_span$YEAR_PREV)
plot_span$tree_count_prev <- mapply(function(p, y) {
  sum(tr$PLOT==p & tr$MEASYEAR==y & tr$STATUSCD==1 &
      is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5)
}, plot_span$PLOT, plot_span$YEAR_PREV)

# Filter to plots with non-empty tree list
plot_span <- plot_span[plot_span$tree_count_prev > 0, ]

write.csv(plot_span, file.path(od, "silc_cfi_longhorizon_pairs_v16.csv"), row.names = FALSE)
cat(sprintf("=== Long-horizon CFI pairs v16: n=%d plots ===\n", nrow(plot_span)))
cat(sprintf("Mean horizon: %.1f yr (range %d-%d)\n",
            mean(plot_span$PERIOD_YR), min(plot_span$PERIOD_YR), max(plot_span$PERIOD_YR)))
cat(sprintf("Mean trees per plot at year_prev: %.0f\n", mean(plot_span$tree_count_prev)))
cat(sprintf("BA_PREV range: %.1f-%.1f ft²/ac\n", min(plot_span$BA_PREV_FT2AC), max(plot_span$BA_PREV_FT2AC)))
cat(sprintf("BA_CURR range: %.1f-%.1f ft²/ac\n", min(plot_span$BA_CURR_FT2AC), max(plot_span$BA_CURR_FT2AC)))
print(head(plot_span[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","BA_PREV_FT2AC","BA_CURR_FT2AC","tree_count_prev")], 10))
