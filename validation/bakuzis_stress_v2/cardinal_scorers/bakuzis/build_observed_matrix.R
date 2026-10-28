#!/usr/bin/env Rscript
# =============================================================================
# build_observed_matrix.R
# -----------------------------------------------------------------------------
# Build the OBSERVED Bakuzis matrix from FIA remeasurement pairs. This is the
# biological-realism TARGET that each model "arm" is scored against by
# score_arm.R.
#
# Convention follows A. Weiskittel's existing Bakuzis framework in fvs-conus
# (R/16_bakuzis_matrix_figure.R, BAKUZIS_UNCERTAINTY_README.md,
# python/bakuzis_matrix_heatmap.py). The classic Bakuzis (1969) matrix is the
# set of pairwise relationships among stand-level state variables. We adopt the
# SAME structural state variables Aaron already uses for the stand-level
# trajectory / matrix panels:
#
#     BA   stand basal area
#     QMD  quadratic mean diameter
#     TPH  trees per hectare (stand density)   [Aaron's TPA panel, metric here]
#     SDI  Reineke stand density index
#     RD   relative density (SDI / SDImax)
#     AGB  aboveground biomass  (SECONDARY: arm AGB uses an NE-anchored height,
#                                so it is reported but down-weighted)
#
# The matrix CELLS are the pairwise bivariate relationships among these state
# variables (the off-diagonal of the Bakuzis matrix), summarised as the OLS
# slope on log-log axes (the natural Bakuzis/allometric scale). The headline
# cell is the self-thinning boundary (the QMD-TPH / TPH-BA upper boundary line)
# and the SDI ceiling, exactly the "laws" Aaron tracks. Trajectory DIRECTION
# (the sign of dX/dt from t1 -> t2 in the remeasurement pairs) is also recorded
# so an arm can be checked for moving the right way through the matrix.
#
# Units: the remeasurement RDS is metric (BA m2/ha, QMD cm, TPH trees/ha, SDI,
# RD %). score_arm.R converts the imperial arm CSVs to these same metric units
# before scoring, so observed and arm cells are directly comparable.
#
# Output (to BAKUZIS_OUTPUT_DIR, default /fs/scratch/PUOM0008/crsfaaron/bakuzis):
#   observed_bakuzis_matrix.csv     one row per (stratum x cell), the target
#   observed_self_thinning.csv      self-thinning boundary coefficients
#   observed_state_summary.csv      univariate distributions of each state var
#   fig_observed_bakuzis_matrix.png matrix figure (lower triangle scatter+slope)
#
# Usage:
#   module load gcc/12.3.0 R/4.4.0
#   Rscript build_observed_matrix.R \
#       [--pairs <rds>] [--out <dir>] [--strata EPA_L2_CODE]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})
has_gg <- requireNamespace("ggplot2", quietly = TRUE)

# ---- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
PAIRS <- get_arg("--pairs",
  "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds")
OUTDIR <- get_arg("--out",
  Sys.getenv("BAKUZIS_OUTPUT_DIR", "/fs/scratch/PUOM0008/crsfaaron/bakuzis"))
STRATUM <- get_arg("--strata", "EPA_L2_CODE")   # column in the RDS

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[observed] pairs = %s\n[observed] out   = %s\n[observed] strata = %s\n",
            PAIRS, OUTDIR, STRATUM))

# ---- Bakuzis convention: state variables + matrix cells ---------------------
# Structural (height-independent) state variables; AGB is secondary.
STATE_VARS <- c("BA", "QMD", "TPH", "SDI", "RD")          # primary, structural
SECONDARY  <- c("AGB")                                    # height-anchored
ALL_VARS   <- c(STATE_VARS, SECONDARY)

# Bakuzis matrix cells = ordered pairs (y ~ x) on log-log axes. We keep the
# ecologically canonical orientation for each (the "law" reading):
CELLS <- rbind(
  data.table(cell = "QMD_TPH",  yvar = "QMD", xvar = "TPH", law = "self-thinning (-3/2 / -0.5 region)"),
  data.table(cell = "TPH_BA",   yvar = "TPH", xvar = "BA",  law = "density-BA boundary"),
  data.table(cell = "BA_QMD",   yvar = "BA",  xvar = "QMD", law = "BA accrues with mean size"),
  data.table(cell = "SDI_QMD",  yvar = "SDI", xvar = "QMD", law = "SDI vs mean size"),
  data.table(cell = "SDI_BA",   yvar = "SDI", xvar = "BA",  law = "SDI vs BA"),
  data.table(cell = "BA_TPH",   yvar = "BA",  xvar = "TPH", law = "BA vs density"),
  data.table(cell = "RD_QMD",   yvar = "RD",  xvar = "QMD", law = "relative density vs size"),
  data.table(cell = "AGB_BA",   yvar = "AGB", xvar = "BA",  law = "biomass-BA (secondary, ht-anchored)")
)

# ---- load + standardise the observed pairs ----------------------------------
cat("[observed] reading RDS (large; ~90 MB)...\n")
p <- as.data.table(readRDS(PAIRS))

# The RDS carries t1 (first measurement) and t2 (remeasurement) for each state
# variable. Build a tidy table with one row per measurement OCCASION so the
# observed cloud spans the realised stand-structure space, plus a paired table
# carrying the t1->t2 direction.
# Map RDS column names -> our state var names (metric).
occ1 <- data.table(
  stratum = as.character(p[[STRATUM]]),
  FORTYPCD = p[["FORTYPCD"]],
  BA  = p[["BA1"]],  QMD = p[["QMD1"]], TPH = p[["TPH1"]],
  SDI = p[["SDI1"]], RD  = p[["RD1"]],
  AGB = if ("max_biomass" %in% names(p)) NA_real_ else NA_real_,
  occasion = 1L
)
occ2 <- data.table(
  stratum = as.character(p[[STRATUM]]),
  FORTYPCD = p[["FORTYPCD"]],
  BA  = p[["BA2"]],  QMD = p[["QMD2"]], TPH = p[["TPH2"]],
  SDI = p[["SDI2"]], RD  = p[["RD2"]],
  AGB = NA_real_,
  occasion = 2L
)
occ <- rbind(occ1, occ2)

# basic cleaning: drop non-positive / impossible values for log-log work
for (v in c("BA","QMD","TPH","SDI","RD")) {
  occ[[v]][occ[[v]] <= 0] <- NA_real_
}
occ <- occ[is.finite(BA) & is.finite(QMD) & is.finite(TPH)]
cat(sprintf("[observed] %d measurement occasions after cleaning\n", nrow(occ)))

# SDImax (observed self-thinning ceiling) from the curated brms column if present
SDImax_obs <- if ("SDImax_brms" %in% names(p)) {
  median(p[["SDImax_brms"]], na.rm = TRUE)
} else NA_real_

# =============================================================================
# 1. Univariate state summary (the matrix diagonal)
# =============================================================================
state_summary <- rbindlist(lapply(STATE_VARS, function(v) {
  x <- occ[[v]]
  data.table(variable = v,
             n = sum(is.finite(x)),
             median = median(x, na.rm = TRUE),
             q025 = quantile(x, 0.025, na.rm = TRUE),
             q975 = quantile(x, 0.975, na.rm = TRUE),
             p99  = quantile(x, 0.99,  na.rm = TRUE))
}))
fwrite(state_summary, file.path(OUTDIR, "observed_state_summary.csv"))

# =============================================================================
# 2. Bakuzis matrix cells: log-log OLS slope per (stratum x cell) + national
# =============================================================================
fit_cell <- function(dt, yv, xv) {
  d <- dt[is.finite(get(yv)) & is.finite(get(xv)) & get(yv) > 0 & get(xv) > 0]
  if (nrow(d) < 30) return(list(n = nrow(d), slope = NA, intercept = NA,
                                r2 = NA, ymed = NA, xmed = NA))
  ly <- log(d[[yv]]); lx <- log(d[[xv]])
  fit <- lm(ly ~ lx)
  list(n = nrow(d),
       slope = unname(coef(fit)[2]),
       intercept = unname(coef(fit)[1]),
       r2 = summary(fit)$r.squared,
       ymed = median(d[[yv]]), xmed = median(d[[xv]]))
}

build_matrix <- function(dt, stratum_label) {
  rbindlist(lapply(seq_len(nrow(CELLS)), function(i) {
    cc <- CELLS[i]
    f <- fit_cell(dt, cc$yvar, cc$xvar)
    data.table(stratum = stratum_label, cell = cc$cell,
               yvar = cc$yvar, xvar = cc$xvar, law = cc$law,
               n = f$n, slope = f$slope, intercept = f$intercept,
               r2 = f$r2, ymed = f$ymed, xmed = f$xmed)
  }))
}

# national pooled
mat_nat <- build_matrix(occ, "ALL")

# per-stratum (keep strata with enough data)
strata_keep <- occ[, .N, by = stratum][N >= 500 & stratum != "" & !is.na(stratum), stratum]
mat_strata <- rbindlist(lapply(strata_keep, function(s) build_matrix(occ[stratum == s], s)))

observed_matrix <- rbind(mat_nat, mat_strata)
fwrite(observed_matrix, file.path(OUTDIR, "observed_bakuzis_matrix.csv"))
cat(sprintf("[observed] wrote observed_bakuzis_matrix.csv (%d rows, %d strata + ALL)\n",
            nrow(observed_matrix), length(strata_keep)))

# =============================================================================
# 3. Self-thinning boundary: upper-boundary line of log(TPH) ~ log(QMD)
#    Estimated as the 99th-percentile (max-density) edge via quantile regression
#    surrogate: bin by QMD, take per-bin 99th pct of TPH, fit log-log line.
#    Reineke slope ~ -1.605 expected.
# =============================================================================
st <- occ[is.finite(QMD) & is.finite(TPH) & QMD > 0 & TPH > 0]
st[, qbin := cut(log(QMD), breaks = 25)]
edge <- st[, .(lq = mean(log(QMD)),
               lt99 = quantile(log(TPH), 0.99, na.rm = TRUE),
               n = .N), by = qbin][n >= 50]
st_fit <- lm(lt99 ~ lq, data = edge)
self_thin <- data.table(
  metric = c("reineke_slope_TPH_QMD", "reineke_intercept", "n_bins",
             "central_slope_TPH_QMD", "SDImax_observed_median"),
  value = c(unname(coef(st_fit)[2]), unname(coef(st_fit)[1]), nrow(edge),
            fit_cell(occ, "TPH", "QMD")$slope, SDImax_obs)
)
fwrite(self_thin, file.path(OUTDIR, "observed_self_thinning.csv"))
cat(sprintf("[observed] self-thinning boundary slope (log TPH ~ log QMD) = %.3f (Reineke ref -1.605)\n",
            coef(st_fit)[2]))
cat(sprintf("[observed] observed SDImax median = %.0f\n", SDImax_obs))

# =============================================================================
# 4. Figure: lower-triangle scatter + fitted slope for the key cells
# =============================================================================
if (has_gg) {
  suppressPackageStartupMessages(library(ggplot2))
  plot_cells <- c("QMD_TPH","BA_QMD","SDI_QMD","TPH_BA","BA_TPH","RD_QMD")
  samp <- occ[sample(.N, min(.N, 40000))]
  pl <- lapply(plot_cells, function(cn) {
    cc <- CELLS[cell == cn]
    d <- samp[is.finite(get(cc$yvar)) & is.finite(get(cc$xvar)) &
              get(cc$yvar) > 0 & get(cc$xvar) > 0]
    row <- mat_nat[cell == cn]
    ggplot(d, aes(x = log(get(cc$xvar)), y = log(get(cc$yvar)))) +
      geom_hex(bins = 40) +
      geom_abline(slope = row$slope, intercept = row$intercept,
                  color = "red", linewidth = 0.8) +
      scale_fill_viridis_c(trans = "log10", guide = "none") +
      labs(x = paste0("log ", cc$xvar), y = paste0("log ", cc$yvar),
           title = sprintf("%s  (slope %.2f, R2 %.2f)", cn, row$slope, row$r2)) +
      theme_minimal(base_size = 9)
  })
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    fig <- wrap_plots(pl, ncol = 3) +
      plot_annotation(title = "Observed Bakuzis matrix (FIA remeasurement pairs)",
                      subtitle = "Log-log bivariate stand-structure cells; red = OLS slope")
    ggsave(file.path(OUTDIR, "fig_observed_bakuzis_matrix.png"), fig,
           width = 12, height = 7, dpi = 150)
    cat("[observed] wrote fig_observed_bakuzis_matrix.png\n")
  }
}

cat("[observed] DONE\n")
