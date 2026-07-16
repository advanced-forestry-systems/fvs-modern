##=============================================================================
## 50m_extract_mortality_plotlevel.R
##
## Mortality residual extractor at PLOT level.
## Groups trees by PLOT_CN, computes observed proportion mortality and
## mean predicted p_mort. Residual on cloglog scale:
##
##   p_obs_plot = sum(died) / sum(trees_plot)         (weighted)
##   p_pred_plot = weighted mean of p_mort_T_hat over plot
##   residual = cloglog(p_obs_clamped) - cloglog(p_pred_clamped)
##   weight = sqrt(n_trees_plot)
##
## CLI: --fit, --meta, --pairs, --traits, --out, --subsample
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

FIT_FILE  <- get_arg("fit")
META_FILE <- get_arg("meta", sub("_fit\\.rds$", "_meta.rds", FIT_FILE %||% ""))
PAIRS_FILE  <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_FILE  <- get_arg("out")
SUBSAMPLE <- as.integer(get_arg("subsample", "0"))

stopifnot(!is.null(FIT_FILE), file.exists(FIT_FILE))
stopifnot(file.exists(META_FILE))
stopifnot(!is.null(OUT_FILE))

cat("== 50m_extract_mortality_plotlevel.R ==\n")
cat("  fit:", FIT_FILE, "\n")
cat("  out:", OUT_FILE, "\n\n")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

fit   <- readRDS(FIT_FILE)
meta  <- readRDS(META_FILE)
pairs <- as.data.table(readRDS(PAIRS_FILE))

sp_levels <- meta$sp_levels
L1_levels <- meta$L1_levels
L2_levels <- meta$L2_levels
L3_levels <- meta$L3_levels
FT_levels <- meta$FT_levels

pairs[, alive := as.integer(TREESTATUS2 == 1)]
pairs[, rd_ratio := sdi_additive1 / SDImax_brms]
pairs[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]
if ("climate_si" %in% names(pairs)) {
  med <- median(pairs$climate_si, na.rm = FALSE)
  pairs[!is.finite(climate_si), climate_si := med]
  pairs[, ln_csi := log(pmax(climate_si, 0.1))]
} else { pairs[, ln_csi := 0] }
pairs[!is.finite(ln_csi), ln_csi := 0]

# Identify plot column (prefer unique plot key, fall back to PLOT)
plot_col <- intersect(c("plot_key","PLT_CN_cond1","PLOT_CN","PLT_CN","plot_id","PLOT"),
                      names(pairs))[1]
if (is.na(plot_col)) stop("No plot column found")
cat("Plot column:", plot_col, "(unique values =", uniqueN(pairs[[plot_col]]), ")\n")

filt <- with(pairs,
  TREESTATUS1 == 1 & !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1,2) &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
  is.finite(BA1) & BA1 >= 0 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 & is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(rd_ratio) & rd_ratio >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0)
d <- pairs[filt]

# Match species (need to skip rows whose species not in base)
d <- d[SPCD %in% sp_levels]
d[, sp_idx := match(SPCD, sp_levels)]
d[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
d[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
d[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
d[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]
d <- d[!is.na(sp_idx) & !is.na(L1_idx) & !is.na(L2_idx) &
        !is.na(L3_idx) & !is.na(FT_idx)]
cat("Tree-level rows:", nrow(d), "\n")

# Coerce NAs on disturbance/treatment indicators to 0 BEFORE aggregation
ind_cols <- c("is_plantation","had_fire_t1","had_insect_t1","had_disease_t1",
              "had_wind_t1","had_harvest_t1","had_cutting_t1","had_site_prep_t1")
for (c in intersect(ind_cols, names(d))) {
  v <- d[[c]]; v[is.na(v)] <- 0L; d[[c]] <- as.integer(v)
}
for (c in grep("_decay_", names(d), value = TRUE)) {
  v <- d[[c]]; v[is.na(v)] <- 0.0; d[[c]] <- as.numeric(v)
}

# Posterior-mean tree-level p_mort_T
get_mean <- function(v) {
  out <- tryCatch(fit$draws(variables = v, format = "draws_matrix"),
                   error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  mean(as.numeric(out))
}
vec_mean <- function(base, n) {
  out <- tryCatch(fit$draws(variables = base, format = "draws_matrix"),
                   error = function(e) NULL)
  if (is.null(out)) return(rep(0, n))
  colMeans(as.matrix(out))
}

trait_effect <- vec_mean("trait_effect", length(sp_levels))
z_L1 <- vec_mean("z_L1", length(L1_levels))
z_L2 <- vec_mean("z_L2", length(L2_levels))
z_L3 <- vec_mean("z_L3", length(L3_levels))
z_FT <- vec_mean("z_FT", length(FT_levels))
b0 <- get_mean("b0"); b1 <- get_mean("b1"); b2 <- get_mean("b2")
b3 <- get_mean("b3"); b4 <- get_mean("b4"); b5 <- get_mean("b5"); b6 <- get_mean("b6")

eta <- b0 + trait_effect[d$sp_idx] +
  z_L1[d$L1_idx] + z_L2[d$L2_idx] + z_L3[d$L3_idx] + z_FT[d$FT_idx] +
  b1 * d$DBH1 + b2 * d$DBH1^2 +
  b3 * d$CR1 + b4 * d$ln_csi +
  b5 * (d$BAL_SW1 + d$BAL_HW1) + b6 * d$sqrt_ba_rd
eta_safe <- pmin(pmax(eta, -20), 5)
d[, p_mort_T := 1 - exp(-exp(eta_safe) * YEARS)]
d[, dead := 1L - alive]

# Aggregate to plot level: take FIRST disturbance/treatment within plot
# (assume one event per plot interval)
plot_id_sym <- as.name(plot_col)
agg <- d[, .(
  n_trees     = .N,
  n_dead      = sum(dead),
  p_obs       = sum(dead) / .N,
  p_pred      = mean(p_mort_T),
  EPA_L1_CODE = first(EPA_L1_CODE),
  L1_idx      = first(L1_idx),
  YEARS       = mean(YEARS),
  is_plantation    = max(is_plantation, na.rm = FALSE),
  had_fire_t1      = max(had_fire_t1, na.rm = FALSE),
  had_insect_t1    = max(had_insect_t1, na.rm = FALSE),
  had_disease_t1   = max(had_disease_t1, na.rm = FALSE),
  had_wind_t1      = max(had_wind_t1, na.rm = FALSE),
  had_harvest_t1   = max(had_harvest_t1, na.rm = FALSE),
  had_cutting_t1   = max(had_cutting_t1, na.rm = FALSE),
  had_site_prep_t1 = max(had_site_prep_t1, na.rm = FALSE),
  dstrb_decay_5yr  = first(dstrb_decay_5yr),
  dstrb_decay_10yr = first(dstrb_decay_10yr),
  dstrb_decay_20yr = first(dstrb_decay_20yr),
  trt_decay_5yr    = first(trt_decay_5yr),
  trt_decay_10yr   = first(trt_decay_10yr),
  trt_decay_20yr   = first(trt_decay_20yr)
), by = eval(plot_col)]

# Filter plots with too few trees
agg <- agg[n_trees >= 3]
cat("Plot-level rows:", nrow(agg), "\n")

# Cloglog residual: cloglog(p) = log(-log(1-p))
# Clamp p away from 0/1: for p_safe in (eps, 1-eps), cloglog is finite
clog <- function(p) {
  p_safe <- pmin(pmax(p, 1e-3), 1 - 1e-3)
  log(-log(1 - p_safe))
}
agg[, residual := clog(p_obs) - clog(p_pred)]
agg[, weight := sqrt(n_trees)]

if (SUBSAMPLE > 0 && SUBSAMPLE < nrow(agg)) {
  set.seed(42)
  agg <- agg[sample(.N, SUBSAMPLE)]
  cat("Subsampled to:", nrow(agg), "plots\n")
}

# Coerce NAs on indicators
for (c in grep("^(is_plantation|had_)", names(agg), value = TRUE)) {
  v <- agg[[c]]; v[is.na(v) | is.infinite(v)] <- 0L; agg[[c]] <- as.integer(v)
}
for (c in grep("_decay_", names(agg), value = TRUE)) {
  v <- agg[[c]]; v[is.na(v)] <- 0.0; agg[[c]] <- as.numeric(v)
}

cat(sprintf("\nPlot-level residual: n=%d mean=%.4f sd=%.4f p01=%.3f p99=%.3f\n",
            sum(is.finite(agg$residual)),
            mean(agg$residual, na.rm = FALSE),
            sd(agg$residual, na.rm = FALSE),
            quantile(agg$residual, 0.01, na.rm = FALSE),
            quantile(agg$residual, 0.99, na.rm = FALSE)))

bundle <- list(
  model = "mort_speciesfree_plotlevel",
  family = "cloglog_plot",
  fit_path = FIT_FILE,
  data = agg,
  sp_levels = sp_levels,
  L1_levels = L1_levels,
  n_rows = nrow(agg),
  resid_sd = sd(agg$residual, na.rm = FALSE)
)

saveRDS(bundle, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")
cat("sigma_resid (plot-level cloglog):", round(bundle$resid_sd, 3), "\n")
cat("Done.\n")
