##=============================================================================
## 34_fit_mortality_speciesfree.R
## Gompit (cloglog) mortality species-free B1 with exposure offset.
## Architecture: trait_effect + L1/L2/L3 ecoregion REs + forest type RE.
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

STAN_FILE <- get_arg("stan_file", "calibration/stan/gompit_mortality_speciesfree.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/mort/speciesfree")
OUT_NAME  <- get_arg("outname",   "mort_speciesfree")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 34_fit_mortality_speciesfree.R ==\n")

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- "calibration/traits/species_traits.rds"

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done. Rows:", nrow(dat), "\n\n")

MIN_OBS_SPECIES <- 5000

# Mortality: alive = (TREESTATUS2 == 1)
# Include trees that started alive (TREESTATUS1==1) and have valid status2
dat[, alive := as.integer(TREESTATUS2 == 1)]

# Climate variable - coalesce NA to median
if ("climate_si" %in% names(dat)) {
  med <- median(dat$climate_si, na.rm = TRUE)
  dat[!is.finite(climate_si), climate_si := med]
  dat[, ln_csi := log(pmax(climate_si, 0.1))]
} else {
  dat[, ln_csi := 0]
}
dat[!is.finite(ln_csi), ln_csi := 0]

# SDI ratio for density-mortality interaction
dat[, rd_ratio := sdi_additive1 / SDImax_brms]
dat[, sqrt_ba_rd := sqrt(pmax(BA1 * 0.2296, 0) * pmax(rd_ratio, 0))]

dat <- dat[
  TREESTATUS1 == 1 &
  !is.na(TREESTATUS2) & TREESTATUS2 %in% c(1, 2) &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  is.finite(BA1) & BA1 >= 0 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 &
  is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(rd_ratio) & rd_ratio >= 0 &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
]
cat("After column filters:", nrow(dat), "rows\n")
cat("Mortality rate:", round(1 - mean(dat$alive), 4), "\n")

sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("After species filter:", nrow(dat), "rows;", nrow(sp_counts), "species\n\n")

sp_levels <- sort(unique(dat$SPCD))
L1_levels <- sort(unique(as.character(dat$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(dat$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(dat$EPA_L3_CODE)))
FT_levels <- sort(unique(as.integer(dat$FORTYPCD_cond1)))

dat[, sp_idx := match(SPCD, sp_levels)]
dat[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
dat[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
dat[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
dat[, FT_idx := match(as.integer(FORTYPCD_cond1), FT_levels)]

trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j])
  if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42)
  idx <- sort(sample.int(nrow(dat), SUBSAMPLE))
  dat <- dat[idx]
  cat("Subsampled to:", nrow(dat), "rows\n\n")
}

stan_data <- list(
  N_obs = nrow(dat),
  N_sp = length(sp_levels),
  N_L1 = length(L1_levels),
  N_L2 = length(L2_levels),
  N_L3 = length(L3_levels),
  N_FT = length(FT_levels),
  P_trait = ncol(W),

  alive = dat$alive,
  T_years = dat$YEARS,
  dbh = dat$DBH1,
  dbh_sq = dat$DBH1^2,
  cr_init = dat$CR1,
  ln_csi = dat$ln_csi,
  bal_metric = (dat$BAL_SW1 + dat$BAL_HW1),
  sqrt_ba_rd = dat$sqrt_ba_rd,

  sp_idx = dat$sp_idx,
  L1_idx = dat$L1_idx,
  L2_idx = dat$L2_idx,
  L3_idx = dat$L3_idx,
  FT_idx = dat$FT_idx,
  W = W
)

cat("=== Stan data ready ===\n")
cat("N_obs   =", stan_data$N_obs, "\n")
cat("Survival rate:", round(mean(stan_data$alive), 4), "\n\n")

mod <- cmdstan_model(STAN_FILE)

if (SMOKE) {
  iter_warmup <- 50; iter_sampling <- 50; chains <- 2
} else {
  iter_warmup <- 1000; iter_sampling <- 1000; chains <- 4
}

t_start <- Sys.time()
fit <- mod$sample(
  data = stan_data, chains = chains, parallel_chains = chains,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  seed = 42, adapt_delta = 0.9, max_treedepth = 10, refresh = 100
)
wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("\nWall:", round(wall_min, 1), "min\n\n")

fit_path  <- file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds"))
meta_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds"))
summ_path <- file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv"))

fit$save_object(fit_path)
vars <- c("b0", paste0("b", 1:6),
          paste0("gamma[", seq_len(ncol(W)), "]"),
          "sigma_L1","sigma_L2","sigma_L3","sigma_FT")
summary_df <- fit$summary(variables = vars, "mean","median","sd",
                          ~quantile(.x, c(0.05, 0.95)),
                          "rhat","ess_bulk","ess_tail")
names(summary_df)[names(summary_df) %in% c("5%","95%")] <- c("q5","q95")
data.table::fwrite(summary_df, summ_path)
saveRDS(list(form = "mort_speciesfree", trait_cols = trait_cols,
             stan_file = STAN_FILE,
             sp_levels = sp_levels, L1_levels = L1_levels,
             L2_levels = L2_levels, L3_levels = L3_levels,
             FT_levels = FT_levels,
             summary = summary_df, n_obs = stan_data$N_obs,
             wall_min = wall_min), meta_path)

cat("=== b coefs ===\n")
print(summary_df[grepl("^b[0-9]", summary_df$variable), ])
cat("\n=== sigmas ===\n")
print(summary_df[grepl("^sigma", summary_df$variable), ])
cat("\nDone.\n")
