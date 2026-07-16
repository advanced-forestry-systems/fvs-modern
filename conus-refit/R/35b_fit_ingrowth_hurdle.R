##=============================================================================
## 35b_fit_ingrowth_hurdle.R
##
## Stage 1 + Stage 2 of the four-stage ingrowth model.
## Hurdle negative binomial: occurrence (logistic) + zero-truncated count
## (NegBin2). Same data inputs as 35_fit_ingrowth_negbinom_v4.R; different
## Stan model (ingrowth_hurdle_v1.stan).
##
## CLI:
##   --stan_file=PATH (default ingrowth_hurdle_v1.stan)
##   --outdir=PATH
##   --outname=NAME
##   --subsample=N
##   --traits=PATH (default species_traits_v2.rds when available)
##   --smoke
##=============================================================================

suppressPackageStartupMessages({
  library(data.table); library(cmdstanr); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/ingrowth_hurdle_v1.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/ingrowth/hurdle")
OUT_NAME  <- get_arg("outname",   "ingrowth_hurdle_v1")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
TRAITS_FILE <- get_arg("traits",
                       if (file.exists("calibration/traits/species_traits_v2.rds"))
                         "calibration/traits/species_traits_v2.rds"
                       else "calibration/traits/species_traits.rds")
SMOKE     <- has_flag("smoke")
MAX_TREEDEPTH <- as.integer(get_arg("max_treedepth", "10"))
ADAPT_DELTA   <- as.numeric(get_arg("adapt_delta", "0.9"))
ITER_WARMUP   <- as.integer(get_arg("iter_warmup",   if (has_flag("smoke")) "50"   else "1000"))
ITER_SAMPLING <- as.integer(get_arg("iter_sampling", if (has_flag("smoke")) "50"   else "1000"))
CHAINS        <- as.integer(get_arg("chains",        if (has_flag("smoke")) "2"    else "4"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 35b_fit_ingrowth_hurdle.R ==\n")
cat("Stan:  ", STAN_FILE, "\n")
cat("Traits:", TRAITS_FILE, "\n\n")

## 1. Load GRM + plot data -------------------------------------------------
GRM_DIR <- "calibration/data/raw_fia"
grm_files <- list.files(GRM_DIR, pattern = "_TREE_GRM_COMPONENT\\.csv$",
                        full.names = TRUE)
cat("GRM files:", length(grm_files), "\n")
gr_cols <- c("TRE_CN","PLT_CN","STATECD","SUBP_COMPONENT_AL_FOREST")
cat("Loading GRM..."); flush.console()
grm <- rbindlist(lapply(grm_files, fread, select=gr_cols, showProgress=FALSE))
cat(" done. Rows:", nrow(grm), "\n")

plot_recruit <- grm[, .(n_recruits = sum(SUBP_COMPONENT_AL_FOREST %in%
                                          c("INGROWTH","REVERSION1","REVERSION2"))),
                    by = PLT_CN]
total_plots <- nrow(plot_recruit)
n_with <- sum(plot_recruit$n_recruits > 0)
cat(sprintf("Plot-level: total = %d  with ingrowth = %d (%.1f%%)\n",
            total_plots, n_with, 100 * n_with / total_plots))

## 2. Load remeasurement pairs to get plot covariates ---------------------
PAIRS_FILE <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
cat("Loading pairs..."); flush.console()
pairs <- as.data.table(readRDS(PAIRS_FILE))
# pairs uses PLT_CN_cond1; rename to PLT_CN for join with GRM data
if (!"PLT_CN" %in% names(pairs) && "PLT_CN_cond1" %in% names(pairs)) {
  setnames(pairs, "PLT_CN_cond1", "PLT_CN")
}
cat(" done\n")

# Plot-level aggregation from tree-level pairs
plot_cov <- pairs[, .(
    EPA_L1_CODE = first(EPA_L1_CODE),
    EPA_L2_CODE = first(EPA_L2_CODE),
    EPA_L3_CODE = first(EPA_L3_CODE),
    years    = first(YEARS),
    cspi     = first(cspi),
    clim_pca1 = first(clim_pca1),
    BA       = sum(BA1[TREESTATUS1 == 1], na.rm = TRUE),
    BAL_mean = mean(BAL1[TREESTATUS1 == 1], na.rm = TRUE),
    RD       = first(SDI1) / first(SDImax_brms),
    ht40     = if ("HT40_DOM_t1" %in% names(.SD)) first(HT40_DOM_t1) else NA_real_,
    dom_spcd = SPCD[which.max(DBH1 * (TREESTATUS1 == 1))]
  ), by = PLT_CN]

dat <- plot_recruit[plot_cov, on = "PLT_CN", nomatch = 0]
cat("Inner-join plots:", nrow(dat), "\n")

## 3. Filter + transformations --------------------------------------------
dat <- dat[is.finite(years) & years >= 1 & years <= 20 &
           is.finite(BA) & BA >= 0 &
           is.finite(RD) & RD > 0 & RD < 3 &
           !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
           !is.na(dom_spcd) & is.finite(cspi)]
dat[is.na(BAL_mean), BAL_mean := 0]
dat[is.na(ht40), ht40 := 5.0]
dat[is.na(clim_pca1), clim_pca1 := 0]

dat[, ln_ba     := log(BA * 0.2296 + 1.0)]
dat[, ln_bal    := log(BAL_mean * 0.2296 + 1.0)]
dat[, ln_csi    := log(pmax(cspi + 1.0, 0.01))]
dat[, log_years := log(years)]
cat("After filter:", nrow(dat), "plots\n")
cat(sprintf("  ingrowth events range: %d-%d, mean = %.2f, %% zero = %.1f%%\n",
            min(dat$n_recruits), max(dat$n_recruits), mean(dat$n_recruits),
            100 * mean(dat$n_recruits == 0)))

## 4. DF + LP variety split (use new SPCDs in trait file) -----------------
traits <- as.data.table(readRDS(TRAITS_FILE))
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  pre <- sum(dat$dom_spcd == 202L)
  dat[dom_spcd == 202L & as.character(EPA_L1_CODE) == "7", dom_spcd := 2020L]
  dat[dom_spcd == 202L, dom_spcd := 2021L]
  cat(sprintf("DF variety split (dom_spcd): %d -> %d coastal + %d rocky\n",
              pre, sum(dat$dom_spcd == 2020L), sum(dat$dom_spcd == 2021L)))
}
if (any(traits$SPCD == 1080L) && any(traits$SPCD == 1081L)) {
  pre <- sum(dat$dom_spcd == 108L)
  dat[dom_spcd == 108L & as.character(EPA_L1_CODE) == "7", dom_spcd := 1080L]
  dat[dom_spcd == 108L, dom_spcd := 1081L]
  cat(sprintf("LP variety split (dom_spcd): %d -> %d shore + %d rocky\n",
              pre, sum(dat$dom_spcd == 1080L), sum(dat$dom_spcd == 1081L)))
}

## 5. Subsample ------------------------------------------------------------
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  dat <- dat[sample(.N, SUBSAMPLE)]
  cat("Subsampled to:", nrow(dat), "plots\n")
}

## 6. Build trait matrix W_dom --------------------------------------------
# Auto-detect traits_v3 (decomposed Potter VCC: CE + S + LAC) vs v2 (composite vuln_score + S)
use_v3_traits <- all(c("climate_exposure", "low_adaptive_cap") %in% names(traits))
if (use_v3_traits) {
  cat("[traits] detected v3 layout: using decomposed Potter components (CE+S+LAC)\n")
  trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                  "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                  "climate_exposure", "sensitivity", "low_adaptive_cap")
} else {
  trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                  "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                  "vulnerability_score", "sensitivity")
}
sp_levels <- sort(unique(dat$dom_spcd))
W_dom <- as.matrix(traits[match(sp_levels, SPCD), trait_cols, with = FALSE])
for (j in seq_len(ncol(W_dom))) {
  na <- is.na(W_dom[, j])
  if (any(na)) W_dom[na, j] <- median(W_dom[!na, j], na.rm = TRUE)
  W_dom[, j] <- (W_dom[, j] - mean(W_dom[, j])) / sd(W_dom[, j])
}

## 7. Indices --------------------------------------------------------------
dat[, dom_sp_idx := match(dom_spcd, sp_levels)]
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
dat <- dat[!is.na(L1_idx) & !is.na(L2_idx) & !is.na(L3_idx) & !is.na(dom_sp_idx)]
cat("\nFinal N_plots:", nrow(dat), "\n")

## 8. Stan data ------------------------------------------------------------
stan_data <- list(
  N_plots = nrow(dat),
  N_sp    = length(sp_levels),
  N_L1    = length(L1_levels),
  N_L2    = length(L2_levels),
  N_L3    = length(L3_levels),
  P_trait = ncol(W_dom),
  n_recruits = as.integer(dat$n_recruits),
  log_years  = dat$log_years,
  ln_ba      = dat$ln_ba,
  ln_bal     = dat$ln_bal,
  rd         = dat$RD,
  ht40       = dat$ht40,
  ln_csi     = dat$ln_csi,
  clim_pca1  = dat$clim_pca1,
  dom_sp_idx = dat$dom_sp_idx,
  L1_idx     = dat$L1_idx,
  L2_idx     = dat$L2_idx,
  L3_idx     = dat$L3_idx,
  W_dom      = W_dom
)
cat("Stan data: N_plots =", stan_data$N_plots,
    ", % zero =", round(100 * mean(stan_data$n_recruits == 0), 1), "%\n\n")

## 9. Compile + sample -----------------------------------------------------
cat("Compiling Stan model:", STAN_FILE, "\n"); flush.console()
mod <- cmdstan_model(STAN_FILE)
chains <- CHAINS
iw     <- ITER_WARMUP
is_    <- ITER_SAMPLING

t0 <- Sys.time()
fit <- mod$sample(data = stan_data, chains = chains, parallel_chains = chains,
                  iter_warmup = iw, iter_sampling = is_,
                  seed = 42, adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH, refresh = 100)
wm <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat("\nWall:", round(wm, 1), "min\n\n")

## 10. Save ----------------------------------------------------------------
fit$save_object(file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds")))
vars <- c("a0_p","a0_c",
          paste0("b", 1:6, "_p"), paste0("b", 1:6, "_c"),
          paste0("gamma_p[", seq_len(ncol(W_dom)), "]"),
          paste0("gamma_c[", seq_len(ncol(W_dom)), "]"),
          "sigma_L1_p","sigma_L2_p","sigma_L3_p",
          "sigma_L1_c","sigma_L2_c","sigma_L3_c","phi")
summ <- fit$summary(variables = vars, "mean","median","sd",
                    ~quantile(.x, c(0.05, 0.95), na.rm = TRUE),
                    "rhat","ess_bulk","ess_tail")
names(summ)[names(summ) %in% c("5%","95%")] <- c("q5","q95")
data.table::fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv")))
saveRDS(list(form="hurdle_v1", trait_cols=trait_cols, stan_file=STAN_FILE,
             sp_levels=sp_levels, L1_levels=L1_levels, L2_levels=L2_levels,
             L3_levels=L3_levels, summary=summ, n_plots=stan_data$N_plots,
             wall_min=wm),
        file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds")))
cat("Done.\n")
