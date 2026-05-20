##=============================================================================
## 35_fit_ingrowth_negbinom.R
## Plot-level ingrowth (recruitment) model for FVS-CONUS
##
## Counts the number of new trees (DBH >= 2.54 cm) that enter the inventory
## between plot t1 and t2. Fits a negative binomial GLM with hierarchical
## ecoregion random effects and a dominant-species trait composition fixed
## effect, with log(years) as an offset for annualization.
##
## Usage:
##   Rscript --vanilla 35_fit_ingrowth_negbinom.R \
##     --stan_file=calibration/stan/ingrowth_negbinom.stan \
##     --outdir=calibration/output/conus/ingrowth \
##     --outname=ingrowth_negbinom_v1
##
## CLI flags:
##   --stan_file=PATH    override default
##   --outdir=PATH       override default
##   --outname=NAME      override default
##   --subsample=N       fit on a random N-plot subsample (default: full)
##   --smoke             50+50 iter, 2 chains, tiny test
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-13
##=============================================================================

library(data.table)
library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/ingrowth_negbinom.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/ingrowth")
OUT_NAME  <- get_arg("outname",   "ingrowth_negbinom_v1")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 35_fit_ingrowth_negbinom.R ==\n")
cat("Stan: ", STAN_FILE, "\n")
cat("Out:  ", OUT_DIR, "/", OUT_NAME, "_*\n", sep = "")
cat("Smoke:", SMOKE, "\n\n")

## 1. Load data ---------------------------------------------------------------

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- "calibration/traits/species_traits.rds"
stopifnot(file.exists(DATA_FILE), file.exists(STAN_FILE), file.exists(TRAITS_FILE))

cat("Loading data ..."); flush.console()
trees <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done\n\n")

## 2. Aggregate to plot level + identify ingrowth ----------------------------
##
## Ingrowth definition: new trees in the inventory at t2 that were NOT live
## or measured at t1. FIA marks these via TREESTATUS1 (0 or NA) combined with
## TREESTATUS2 == 1 and DBH2 above the measurement threshold. We approximate
## by counting trees where t1 status indicates "not yet present" and t2 has
## a valid DBH above ~2.54 cm.

cat("Identifying ingrowth events ...\n"); flush.console()

# Ingrowth flag at the tree level
trees[, is_ingrowth := (is.na(TREESTATUS1) | TREESTATUS1 == 0) &
                       TREESTATUS2 == 1 &
                       !is.na(DBH2) & DBH2 >= 2.54]

# Aggregate to plot level. Key columns to retain: plot_key, EPA codes,
# overstory state at t1, climate.
plot_dat <- trees[, .(
  n_recruits      = sum(is_ingrowth, na.rm = TRUE),
  ba_t1           = sum(BA1[TREESTATUS1 == 1], na.rm = TRUE),
  bal_mean_t1     = mean(BAL1[TREESTATUS1 == 1], na.rm = TRUE),
  qmd_t1          = first(QMD1),
  sdi_t1          = first(SDI1),
  rd_t1           = first(RD1),
  stand_age       = first(STDAGE),
  years           = first(YEARS),
  cspi            = first(cspi),
  bgi             = first(bgi),
  clim_pca1       = first(clim_pca1),
  clim_pca2       = first(clim_pca2),
  EPA_L1_CODE     = first(EPA_L1_CODE),
  EPA_L2_CODE     = first(EPA_L2_CODE),
  EPA_L3_CODE     = first(EPA_L3_CODE),
  dom_spcd        = SPCD[which.max(DBH1 * (TREESTATUS1 == 1))]
), by = plot_key]

cat("  total plots:", nrow(plot_dat), "\n")
cat("  plots with >= 1 ingrowth:", sum(plot_dat$n_recruits > 0), "\n")
cat("  mean ingrowth per plot:", round(mean(plot_dat$n_recruits), 2), "\n")
cat("  max ingrowth per plot:", max(plot_dat$n_recruits), "\n\n")

## 3. Filter and transform ---------------------------------------------------

CSPI_SHIFT <- 1.0

plot_dat <- plot_dat[
  is.finite(years) & years >= 1 & years <= 20 &
  is.finite(ba_t1) & ba_t1 >= 0 &
  is.finite(rd_t1) & rd_t1 > 0 &
  is.finite(cspi) &
  !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
  !is.na(EPA_L2_CODE) & EPA_L2_CODE != "" &
  !is.na(EPA_L3_CODE) & EPA_L3_CODE != "" &
  !is.na(dom_spcd)
]
plot_dat[is.na(bal_mean_t1), bal_mean_t1 := 0]
plot_dat[is.na(stand_age), stand_age := 0]
plot_dat[is.na(clim_pca1), clim_pca1 := 0]

plot_dat[, ln_ba     := log(ba_t1 * 0.2296 + 1.0)]
plot_dat[, ln_bal    := log(bal_mean_t1 * 0.2296 + 1.0)]
plot_dat[, ln_csi    := log(pmax(cspi + CSPI_SHIFT, 0.01))]
plot_dat[, log_years := log(years)]

cat("  plots after filtering:", nrow(plot_dat), "\n")

## 4. Build dominant-species trait matrix W_dom ------------------------------
trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
traits_sub <- traits[, c("SPCD", trait_cols), with = FALSE]
for (col in trait_cols) {
  med <- median(traits_sub[[col]], na.rm = TRUE)
  traits_sub[is.na(get(col)), (col) := med]
}
trait_means <- sapply(traits_sub[, trait_cols, with = FALSE], mean)
trait_sds   <- sapply(traits_sub[, trait_cols, with = FALSE], sd)

# Match each plot's dominant species to its trait row, then standardize
plot_dat <- merge(plot_dat, traits_sub, by.x = "dom_spcd", by.y = "SPCD",
                  all.x = TRUE)
for (col in trait_cols) {
  med <- median(plot_dat[[col]], na.rm = TRUE)
  plot_dat[is.na(get(col)), (col) := med]
}
W_dom <- as.matrix(plot_dat[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W_dom))) {
  W_dom[, j] <- (W_dom[, j] - trait_means[j]) / trait_sds[j]
}

## 5. Indices ----------------------------------------------------------------
L1_levels <- sort(unique(as.character(plot_dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(plot_dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(plot_dat$EPA_L3_CODE)))
plot_dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
plot_dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
plot_dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]

## 6. Subsample if requested -------------------------------------------------
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(plot_dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(plot_dat), SUBSAMPLE))
  plot_dat <- plot_dat[idx]
  W_dom <- W_dom[idx, , drop = FALSE]
  cat("  subsampled to:", nrow(plot_dat), "plots\n")
}

## 7. Stan data --------------------------------------------------------------
stan_data <- list(
  N_plots = nrow(plot_dat),
  N_L1    = length(L1_levels),
  N_L2    = length(L2_levels),
  N_L3    = length(L3_levels),
  P_trait = ncol(W_dom),

  n_recruits = as.integer(plot_dat$n_recruits),
  log_years  = plot_dat$log_years,

  ln_ba     = plot_dat$ln_ba,
  ln_bal    = plot_dat$ln_bal,
  rd        = plot_dat$rd_t1,
  ln_csi    = plot_dat$ln_csi,
  stand_age = plot_dat$stand_age,
  clim_pca1 = plot_dat$clim_pca1,

  L1_idx = plot_dat$L1_idx,
  L2_idx = plot_dat$L2_idx,
  L3_idx = plot_dat$L3_idx,

  W_dom = W_dom
)

cat("\n=== Stan data ready ===\n")
cat("N_plots =", stan_data$N_plots, "\n")
cat("N_L1    =", stan_data$N_L1, "\n")
cat("N_L2    =", stan_data$N_L2, "\n")
cat("N_L3    =", stan_data$N_L3, "\n")
cat("P_trait =", stan_data$P_trait, "\n")
cat("mean ingrowth:", round(mean(stan_data$n_recruits), 2),
    "max:", max(stan_data$n_recruits), "\n\n")

## 8. Compile + sample -------------------------------------------------------
cat("Compiling Stan model:", STAN_FILE, "\n"); flush.console()
mod <- cmdstan_model(STAN_FILE)

if (SMOKE) {
  iter_warmup   <- 50
  iter_sampling <- 50
  chains        <- 2
} else {
  iter_warmup   <- 1000
  iter_sampling <- 1000
  chains        <- 4
}

t_start <- Sys.time()
fit <- mod$sample(
  data            = stan_data,
  chains          = chains,
  parallel_chains = chains,
  iter_warmup     = iter_warmup,
  iter_sampling   = iter_sampling,
  seed            = 42,
  adapt_delta     = 0.9,
  max_treedepth   = 10,
  refresh         = 100
)
wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("\nSampling wall time:", round(wall_min, 1), "min\n\n")

## 9. Save -------------------------------------------------------------------
fit_path  <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
meta_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

cat("Saving fit ..."); flush.console()
fit$save_object(fit_path)

vars <- c(paste0("b", 0:6), paste0("gamma[", seq_len(ncol(W_dom)), "]"),
          "sigma_L1", "sigma_L2", "sigma_L3", "phi")
summary_df <- fit$summary(variables = vars,
                          "mean", "median", "sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat", "ess_bulk", "ess_tail")
names(summary_df)[names(summary_df) %in% c("5%", "95%")] <- c("q5", "q95")
data.table::fwrite(summary_df, summ_path)

meta <- list(
  form       = "negbinom",
  variant    = "plot_level_ingrowth",
  site_var   = "cspi",
  trait_cols = trait_cols,
  stan_file  = STAN_FILE,
  L1_levels  = L1_levels,
  L2_levels  = L2_levels,
  L3_levels  = L3_levels,
  summary    = summary_df,
  n_plots    = stan_data$N_plots,
  n_zero     = sum(stan_data$n_recruits == 0),
  mean_ingrowth = mean(stan_data$n_recruits),
  wall_min   = wall_min
)
saveRDS(meta, meta_path)

cat("Saved meta and summary.\n\n")

cat("=== b coefficient posterior summary ===\n")
print(summary_df[grepl("^b[0-9]", summary_df$variable), ])
cat("\n=== gamma trait coefficients ===\n")
print(summary_df[grepl("^gamma", summary_df$variable), ])
cat("\n=== sigma scales + phi ===\n")
print(summary_df[grepl("^sigma|^phi", summary_df$variable), ])

cat("\nDone.\n")
