##=============================================================================
## 50i_extract_ingrowth_residuals.R
##
## Modifier-residual extractor for ingrowth v4 NB fit. Plot-level residuals
## on log(rate) scale (Anscombe-style transform for NB):
##
##   eta_full = eta_base + log(years)
##   residual = log(n_recruits + 0.5) - eta_full
##   weight   = sqrt(years)
##
## Disturbance/treatment indicators aggregated to plot level by MAX across
## trees. PLT_CN kept as character to preserve 14-digit precision.
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(bit64)     # for integer64 character conversion
  library(cmdstanr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

FIT_FILE   <- get_arg("fit")
META_FILE  <- get_arg("meta", sub("_fit\\.rds$", "_meta.rds", FIT_FILE %||% ""))
PAIRS_FILE <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
GRM_DIR    <- get_arg("grm_dir", "calibration/data/raw_fia")
OUT_FILE   <- get_arg("out")
SUBSAMPLE  <- as.integer(get_arg("subsample", "0"))

stopifnot(!is.null(FIT_FILE), file.exists(FIT_FILE))
stopifnot(file.exists(META_FILE))
stopifnot(!is.null(OUT_FILE))

cat("== 50i_extract_ingrowth_residuals.R ==\n")
cat("  fit:", FIT_FILE, "\n  out:", OUT_FILE, "\n\n")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

fit  <- readRDS(FIT_FILE)
meta <- readRDS(META_FILE)

cat("Loading pairs ..."); flush.console()
trees <- as.data.table(readRDS(PAIRS_FILE))
cat(" done\n")

# Coerce indicator columns to integer (no NAs) BEFORE aggregation
ind_cols <- c("is_plantation", "had_fire_t1", "had_insect_t1", "had_disease_t1",
              "had_wind_t1", "had_harvest_t1", "had_cutting_t1",
              "had_site_prep_t1")
for (col in intersect(ind_cols, names(trees))) {
  v <- trees[[col]]; v[is.na(v)] <- 0L; trees[[col]] <- as.integer(v)
}
decay_cols <- grep("_decay_", names(trees), value = TRUE)
for (col in decay_cols) {
  v <- trees[[col]]; v[is.na(v)] <- 0.0; trees[[col]] <- as.numeric(v)
}

# PLT_CN_cond1 may be integer64; convert AFTER aggregation (matching v4)

plot_covs <- trees[, .(
  ba_t1       = first(BA1) * 0.2296,
  bal_mean_t1 = mean(BAL1, na.rm = TRUE) * 0.2296,
  ht40_t1     = first(HT40_1) * 0.3048,
  htlorey_t1  = first(HTlorey1) * 0.3048,
  rd_sdimax   = first(rd_sdimax),
  years       = first(YEARS),
  climate_si  = first(climate_si),
  clim_pca1   = first(clim_pca1),
  EPA_L1_CODE = first(EPA_L1_CODE),
  EPA_L2_CODE = first(EPA_L2_CODE),
  EPA_L3_CODE = first(EPA_L3_CODE),
  is_plantation    = as.integer(max(is_plantation)),
  had_fire_t1      = as.integer(max(had_fire_t1)),
  had_insect_t1    = as.integer(max(had_insect_t1)),
  had_disease_t1   = as.integer(max(had_disease_t1)),
  had_wind_t1      = as.integer(max(had_wind_t1)),
  had_harvest_t1   = as.integer(max(had_harvest_t1)),
  had_cutting_t1   = as.integer(max(had_cutting_t1)),
  had_site_prep_t1 = as.integer(max(had_site_prep_t1)),
  dstrb_decay_5yr  = as.numeric(first(dstrb_decay_5yr)),
  dstrb_decay_10yr = as.numeric(first(dstrb_decay_10yr)),
  dstrb_decay_20yr = as.numeric(first(dstrb_decay_20yr)),
  trt_decay_5yr    = as.numeric(first(trt_decay_5yr)),
  trt_decay_10yr   = as.numeric(first(trt_decay_10yr)),
  trt_decay_20yr   = as.numeric(first(trt_decay_20yr))
), by = PLT_CN_cond1]
setnames(plot_covs, "PLT_CN_cond1", "PLT_CN")
# Use bit64::as.character.integer64 to get actual integer string, then numeric
plot_covs[, PLT_CN := as.numeric(bit64::as.character.integer64(PLT_CN))]
cat("Plot covariates aggregated:", nrow(plot_covs), "plots\n")
cat("  plot_covs PLT_CN sample:", head(plot_covs$PLT_CN, 3), "\n")

# Load GRM ingrowth counts (PLT_CN as character)
cat("Loading GRM ingrowth counts ..."); flush.console()
grm_files <- list.files(GRM_DIR, pattern = "_TREE_GRM_COMPONENT\\.csv$",
                        full.names = TRUE)
# Skip known-corrupt files (HTML placeholders)
grm_files <- grm_files[file.info(grm_files)$size > 10000]
cat(" reading", length(grm_files), "files ...")
gr_cols <- c("TRE_CN","PLT_CN","STATECD","SUBP_COMPONENT_AL_FOREST")
read_safe <- function(f) {
  tryCatch(fread(f, select = gr_cols, showProgress = FALSE),
           error = function(e) NULL)
}
grm <- rbindlist(lapply(grm_files, read_safe))
ingrowth <- grm[SUBP_COMPONENT_AL_FOREST == "INGROWTH"]
plot_counts <- ingrowth[, .(n_recruits = .N), by = PLT_CN]
all_plots <- unique(grm[, .(PLT_CN)])
all_plots <- merge(all_plots, plot_counts, by = "PLT_CN", all.x = TRUE)
all_plots[is.na(n_recruits), n_recruits := 0L]
all_plots[, PLT_CN := as.numeric(bit64::as.character.integer64(PLT_CN))]
rm(grm, ingrowth); gc()
cat(" done. plots with ingrowth:", sum(all_plots$n_recruits > 0), "\n")
cat("  all_plots PLT_CN sample:", head(all_plots$PLT_CN, 3), "\n")
cat("  intersect first 5000:",
    length(intersect(head(plot_covs$PLT_CN, 5000), head(all_plots$PLT_CN, 5000))),
    "\n")

dat <- merge(all_plots, plot_covs, by = "PLT_CN")
cat("After join:", nrow(dat), "plots\n")

L1_levels <- meta$L1_levels
L2_levels <- meta$L2_levels
L3_levels <- meta$L3_levels

dat <- dat[
  is.finite(years) & years >= 1 & years <= 20 &
  is.finite(ba_t1) & ba_t1 >= 0 &
  is.finite(rd_sdimax) & rd_sdimax > 0 & rd_sdimax < 3.0 &
  is.finite(ht40_t1) & ht40_t1 > 0 &
  is.finite(climate_si) & climate_si > 0 &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE != ""
]
cat("After filters:", nrow(dat), "plots\n")

dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
dat <- dat[!is.na(L1_idx) & !is.na(L2_idx) & !is.na(L3_idx)]
cat("After level match:", nrow(dat), "plots\n")

if (SUBSAMPLE > 0 && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  dat <- dat[sample(.N, SUBSAMPLE)]
  cat("Subsampled to:", nrow(dat), "plots\n")
}

dat[, ln_ba := log(ba_t1 + 1)]
dat[, ln_bal := log(bal_mean_t1 + 1)]
dat[, ln_csi := log(pmax(climate_si, 0.1))]
dat[, ln_ht40 := log(ht40_t1 + 1)]
dat[, log_years := log(years)]
dat[!is.finite(clim_pca1), clim_pca1 := 0]

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

b0 <- get_mean("b0"); b1 <- get_mean("b1"); b2 <- get_mean("b2")
b3 <- get_mean("b3"); b4 <- get_mean("b4"); b5 <- get_mean("b5")
b6 <- get_mean("b6")
z_L1 <- vec_mean("z_L1", length(L1_levels))
z_L2 <- vec_mean("z_L2", length(L2_levels))
z_L3 <- vec_mean("z_L3", length(L3_levels))

eta_base <- b0 +
  z_L1[dat$L1_idx] + z_L2[dat$L2_idx] + z_L3[dat$L3_idx] +
  b1 * dat$ln_ba +
  b2 * dat$ln_bal +
  b3 * dat$rd_sdimax +
  b4 * dat$ln_csi +
  b5 * dat$ln_ht40 +
  b6 * dat$clim_pca1

eta_full <- eta_base + dat$log_years
dat[, eta_full := eta_full]
dat[, residual := log(n_recruits + 0.5) - eta_full]
dat[, weight := sqrt(years)]

resid <- dat$residual
cat(sprintf("\nResidual: n=%d mean=%.4f sd=%.4f p01=%.3f p99=%.3f\n",
            sum(is.finite(resid)),
            mean(resid, na.rm = TRUE),
            sd(resid, na.rm = TRUE),
            quantile(resid, 0.01, na.rm = TRUE),
            quantile(resid, 0.99, na.rm = TRUE)))

keep <- intersect(c("PLT_CN", "L1_idx", "years", "n_recruits",
                     "is_plantation", "had_fire_t1", "had_insect_t1",
                     "had_disease_t1", "had_wind_t1", "had_harvest_t1",
                     "had_cutting_t1", "had_site_prep_t1",
                     "dstrb_decay_5yr", "dstrb_decay_10yr", "dstrb_decay_20yr",
                     "trt_decay_5yr", "trt_decay_10yr", "trt_decay_20yr",
                     "eta_full", "residual", "weight"),
                   names(dat))
out <- dat[, ..keep]

bundle <- list(
  model     = "ingrowth_v4_plotlevel",
  family    = "nb_log_link",
  fit_path  = FIT_FILE,
  data      = out,
  L1_levels = L1_levels,
  n_rows    = nrow(out),
  resid_sd  = sd(resid, na.rm = TRUE)
)
saveRDS(bundle, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\nsigma_resid:", round(bundle$resid_sd, 3), "\n")

cat("\n=== Plot-level disturbance/treatment prevalence ===\n")
for (col in grep("^(is_plantation|had_)", names(out), value = TRUE)) {
  cat(sprintf("  %-20s : %.4f\n", col, mean(out[[col]])))
}
