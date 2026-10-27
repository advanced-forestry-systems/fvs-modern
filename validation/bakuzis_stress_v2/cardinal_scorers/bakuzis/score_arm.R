#!/usr/bin/env Rscript
# =============================================================================
# score_arm.R
# -----------------------------------------------------------------------------
# Compute one model ARM's Bakuzis matrix from its per-cycle stand-metric CSVs
# and score each matrix cell against the OBSERVED reference built by
# build_observed_matrix.R.
#
# An "arm" is one model build: a (VARIANT-set x CONFIG) combination, e.g.
#   - engine default      CONFIG = default
#   - engine calibrated   CONFIG = calibrated
#   - conus species-dep   CONFIG = conus_b2
#   - conus species-free  CONFIG = conus_b1
#
# Arm inputs are per-cycle stand metrics (one row per STAND_CN x PROJ_YEAR):
#   STAND_CN, VARIANT, CONFIG, YEAR, PROJ_YEAR, AGB_TONS_AC, BA_FT2AC, QMD_IN, TPH
# (engine arm v2 CSVs use the same schema; if column names differ slightly they
#  are auto-mapped below.)
#
# Bakuzis convention adopted = the same as build_observed_matrix.R:
#   structural state vars BA, QMD, TPH, SDI, RD (AGB secondary, ht-anchored),
#   matrix cells = log-log OLS slope of canonical bivariate relationships,
#   headline = self-thinning boundary slope + SDI ceiling.
#
# SCORING (per cell):
#   slope_dev   = arm_slope - obs_slope               (signed)
#   abs_dev     = |slope_dev|
#   pct_dev     = 100 * slope_dev / |obs_slope|
#   sign_match  = sign(arm_slope) == sign(obs_slope)  (right direction?)
#   cell_score  = exp(-|slope_dev| / 0.5)             (1 = perfect, ->0 worse)
#                 0.5 = one "law-scale" slope unit; structural cells score on
#                 this. AGB cell is flagged secondary and excluded from the
#                 headline realism index.
# Arm realism index = mean(cell_score) over the structural cells, plus a
# self-thinning sub-score = exp(-|arm_reineke - obs_reineke| / 0.3).
#
# Units: arm CSV is imperial; converted to the observed metric units before
# any cell is computed:
#   BA  ft2/ac -> m2/ha   (* 0.2295684)
#   QMD inches -> cm      (* 2.54)
#   TPH already trees/ha  (no change)
#   SDI computed metric (Reineke): SDI = TPH * (QMD_cm / 25)^1.605
#   RD  = 100 * SDI / SDImax     (SDImax from observed reference)
#   AGB tons/ac -> Mg/ha  (* 2.2417)   [secondary]
#
# Usage:
#   module load gcc/12.3.0 R/4.4.0
#   Rscript score_arm.R --arm-name <label> \
#       --glob '<dir>/conus_eq_*_b1_metrics.csv' \
#       [--config conus_b1] \
#       [--observed <dir>/observed_bakuzis_matrix.csv] \
#       [--self-thin <dir>/observed_self_thinning.csv] \
#       [--out <dir>]
#
# You may pass --files 'f1,f2,...' instead of --glob.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ---- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
ARM      <- get_arg("--arm-name", "arm")
GLOB     <- get_arg("--glob", NULL)
FILES    <- get_arg("--files", NULL)
CONFIG   <- get_arg("--config", NULL)             # optional filter on CONFIG col
OUTDIR   <- get_arg("--out", Sys.getenv("BAKUZIS_OUTPUT_DIR",
                    "/fs/scratch/PUOM0008/crsfaaron/bakuzis"))
OBS      <- get_arg("--observed",  file.path(OUTDIR, "observed_bakuzis_matrix.csv"))
OBS_ST   <- get_arg("--self-thin", file.path(OUTDIR, "observed_self_thinning.csv"))
YEAR_MIN <- as.numeric(get_arg("--proj-year-min", "0"))   # keep cycles >= this

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- conversions ------------------------------------------------------------
FT2AC_TO_M2HA <- 0.2295684
IN_TO_CM      <- 2.54
TONSAC_TO_MGHA<- 2.2417

# ---- Bakuzis cell definitions (MUST match build_observed_matrix.R) ----------
CELLS <- rbind(
  data.table(cell="QMD_TPH", yvar="QMD", xvar="TPH", secondary=FALSE),
  data.table(cell="TPH_BA",  yvar="TPH", xvar="BA",  secondary=FALSE),
  data.table(cell="BA_QMD",  yvar="BA",  xvar="QMD", secondary=FALSE),
  data.table(cell="SDI_QMD", yvar="SDI", xvar="QMD", secondary=FALSE),
  data.table(cell="SDI_BA",  yvar="SDI", xvar="BA",  secondary=FALSE),
  data.table(cell="BA_TPH",  yvar="BA",  xvar="TPH", secondary=FALSE),
  data.table(cell="RD_QMD",  yvar="RD",  xvar="QMD", secondary=FALSE),
  data.table(cell="AGB_BA",  yvar="AGB", xvar="BA",  secondary=TRUE)
)

# ---- load arm CSVs ----------------------------------------------------------
if (!is.null(GLOB)) {
  arm_files <- Sys.glob(GLOB)
} else if (!is.null(FILES)) {
  arm_files <- strsplit(FILES, ",")[[1]]
} else {
  stop("provide --glob or --files")
}
arm_files <- arm_files[file.exists(arm_files)]
if (length(arm_files) == 0) stop("no arm files matched")
cat(sprintf("[score:%s] %d arm files\n", ARM, length(arm_files)))

read_arm <- function(f) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  nm <- toupper(names(d)); setnames(d, names(d), nm)
  # auto-map common column variants
  ren <- c(BA_FT2AC="BA_FT2AC", QMD_IN="QMD_IN", TPH="TPH",
           AGB_TONS_AC="AGB_TONS_AC", PROJ_YEAR="PROJ_YEAR", CONFIG="CONFIG",
           VARIANT="VARIANT", STAND_CN="STAND_CN")
  for (a in c("BA","BA_FT2_AC")) if (a %in% nm && !"BA_FT2AC" %in% nm) setnames(d, a, "BA_FT2AC")
  for (a in c("QMD","QMD_INCHES")) if (a %in% nm && !"QMD_IN" %in% nm) setnames(d, a, "QMD_IN")
  d
}
arm <- rbindlist(lapply(arm_files, read_arm), use.names = TRUE, fill = TRUE)
if (!is.null(CONFIG) && "CONFIG" %in% names(arm)) {
  .cfg <- CONFIG; arm <- arm[CONFIG == .cfg]
}
if ("PROJ_YEAR" %in% names(arm)) arm <- arm[is.finite(PROJ_YEAR) & PROJ_YEAR >= YEAR_MIN]
cat(sprintf("[score:%s] %d arm stand-cycle rows\n", ARM, nrow(arm)))

# ---- compute metric structural state variables ------------------------------
arm[, BA  := BA_FT2AC * FT2AC_TO_M2HA]
arm[, QMD := QMD_IN * IN_TO_CM]
arm[, TPH := as.numeric(TPH)]
arm[, AGB := if ("AGB_TONS_AC" %in% names(arm)) AGB_TONS_AC * TONSAC_TO_MGHA else NA_real_]
# Reineke SDI (metric, ref dia 25 cm), and RD against observed SDImax
arm[, SDI := TPH * (QMD / 25)^1.605]

# ---- observed reference -----------------------------------------------------
obs_mat <- fread(OBS)[stratum == "ALL"]
obs_st  <- if (file.exists(OBS_ST)) fread(OBS_ST) else NULL
SDImax  <- if (!is.null(obs_st)) {
  v <- obs_st[metric == "SDImax_observed_median", value]; if (length(v)) v[1] else 1000
} else 1000
arm[, RD := 100 * SDI / SDImax]
obs_reineke <- if (!is.null(obs_st)) obs_st[metric == "reineke_slope_TPH_QMD", value][1] else -1.605

# ---- fit arm cells ----------------------------------------------------------
fit_cell <- function(dt, yv, xv) {
  d <- dt[is.finite(get(yv)) & is.finite(get(xv)) & get(yv) > 0 & get(xv) > 0]
  if (nrow(d) < 30) return(list(n=nrow(d), slope=NA, intercept=NA, r2=NA))
  f <- lm(log(d[[yv]]) ~ log(d[[xv]]))
  list(n=nrow(d), slope=unname(coef(f)[2]), intercept=unname(coef(f)[1]),
       r2=summary(f)$r.squared)
}

scored <- rbindlist(lapply(seq_len(nrow(CELLS)), function(i) {
  cc <- CELLS[i]
  f  <- fit_cell(arm, cc$yvar, cc$xvar)
  o  <- obs_mat[cell == cc$cell]
  obs_slope <- if (nrow(o)) o$slope else NA_real_
  slope_dev <- f$slope - obs_slope
  data.table(
    arm = ARM, cell = cc$cell, yvar = cc$yvar, xvar = cc$xvar,
    secondary = cc$secondary,
    n = f$n, arm_slope = f$slope, arm_r2 = f$r2,
    obs_slope = obs_slope,
    slope_dev = slope_dev,
    abs_dev = abs(slope_dev),
    pct_dev = 100 * slope_dev / abs(obs_slope),
    sign_match = is.finite(f$slope) & is.finite(obs_slope) &
                 (sign(f$slope) == sign(obs_slope)),
    cell_score = ifelse(is.finite(slope_dev), exp(-abs(slope_dev) / 0.5), NA_real_)
  )
}))

# ---- arm self-thinning + SDI ceiling ---------------------------------------
st <- arm[is.finite(QMD) & is.finite(TPH) & QMD > 0 & TPH > 0]
arm_reineke <- NA_real_
if (nrow(st) > 200) {
  st[, qbin := cut(log(QMD), breaks = 25)]
  edge <- st[, .(lq = mean(log(QMD)), lt99 = quantile(log(TPH), 0.99, na.rm=TRUE),
                 n = .N), by = qbin][n >= 30]
  if (nrow(edge) >= 5) arm_reineke <- coef(lm(lt99 ~ lq, data = edge))[2]
}
arm_sdi_ceiling <- quantile(arm$SDI, 0.99, na.rm = TRUE)
obs_sdimax <- if (!is.null(obs_st)) obs_st[metric=="SDImax_observed_median", value][1] else NA

# ---- realism index ----------------------------------------------------------
struct <- scored[secondary == FALSE]
realism_index <- mean(struct$cell_score, na.rm = TRUE)
sign_rate     <- mean(struct$sign_match, na.rm = TRUE)
self_thin_score <- exp(-abs(arm_reineke - obs_reineke) / 0.3)

summary_dt <- data.table(
  arm = ARM,
  n_stand_cycle_rows = nrow(arm),
  realism_index = realism_index,          # 0..1, structural cells only
  sign_match_rate = sign_rate,
  self_thin_score = unname(self_thin_score),
  arm_reineke_slope = unname(arm_reineke),
  obs_reineke_slope = obs_reineke,
  arm_sdi_ceiling = unname(arm_sdi_ceiling),
  obs_sdimax = obs_sdimax,
  agb_cell_score = scored[cell=="AGB_BA", cell_score]
)

# ---- write ------------------------------------------------------------------
safe <- gsub("[^A-Za-z0-9_.-]", "_", ARM)
fwrite(scored,     file.path(OUTDIR, sprintf("arm_%s_bakuzis_cells.csv", safe)))
fwrite(summary_dt, file.path(OUTDIR, sprintf("arm_%s_bakuzis_summary.csv", safe)))

cat("\n================ ARM SCORE:", ARM, "================\n")
print(scored[, .(cell, arm_slope = round(arm_slope,3), obs_slope = round(obs_slope,3),
                 dev = round(slope_dev,3), sign_match, score = round(cell_score,3),
                 secondary)])
cat(sprintf("\nrealism_index (structural)= %.3f   sign_match_rate= %.2f\n",
            realism_index, sign_rate))
cat(sprintf("self-thinning: arm Reineke= %.3f  obs= %.3f  score= %.3f\n",
            arm_reineke, obs_reineke, self_thin_score))
cat(sprintf("SDI ceiling:   arm= %.0f  obs SDImax= %.0f\n",
            arm_sdi_ceiling, obs_sdimax))
cat(sprintf("[score:%s] wrote arm_%s_bakuzis_{cells,summary}.csv\n", ARM, safe))
