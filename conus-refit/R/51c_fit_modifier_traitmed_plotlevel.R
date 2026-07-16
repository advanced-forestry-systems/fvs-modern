##=============================================================================
## 51c_fit_modifier_traitmed_plotlevel.R
##
## Fit modifier_traitmed_plotlevel.stan on plot-level residual bundles
## (mortality, ingrowth, etc). Builds a plot-level trait matrix W_plot
## by aggregating each plot's per-tree trait values (tree-count weighted
## or just first(SPCD)-based dominant-species traits).
##
## CLI:
##   --residuals=PATH (plot-level residual bundle RDS)
##   --pairs=PATH (matched-pairs for plot composition lookup)
##   --traits=PATH (species traits RDS)
##   --stan_file=PATH (default modifier_traitmed_plotlevel.stan)
##   --lambda=N (5/10/20)
##   --subsample=N
##   --out=PATH
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(bit64)
  library(cmdstanr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

RESID_FILE <- get_arg("residuals")
PAIRS_FILE <- get_arg("pairs", "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits.rds")
STAN_FILE  <- get_arg("stan_file", "calibration/stan/modifier_traitmed_plotlevel.stan")
LAMBDA     <- as.integer(get_arg("lambda", "10"))
SUBSAMPLE  <- as.integer(get_arg("subsample", "30000"))
OUT_DIR    <- get_arg("out", "calibration/output/conus/mort/modifier_traitmed_plotlevel")

stopifnot(!is.null(RESID_FILE), file.exists(RESID_FILE))
stopifnot(file.exists(STAN_FILE))
stopifnot(LAMBDA %in% c(5L, 10L, 20L))

cat("== 51c_fit_modifier_traitmed_plotlevel.R ==\n")
cat("  residuals:", RESID_FILE, "\n")
cat("  stan:", STAN_FILE, "\n")
cat("  lambda:", LAMBDA, "yr\n")
cat("  out:", OUT_DIR, "\n\n")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(RESID_FILE)
d <- bundle$data
setDT(d)
cat("Residual bundle:", bundle$model, " family=", bundle$family,
    " rows=", nrow(d), "\n", sep = "")

# Drop NA residuals
d <- d[is.finite(residual)]
cat("After NA filter:", nrow(d), "rows\n")

# Identify plot key column - varies by bundle
plot_key <- if ("PLT_CN" %in% names(d)) "PLT_CN" else {
  intersect(c("plot_key","PLT_CN_cond1"), names(d))[1]
}
if (is.na(plot_key)) stop("No plot key found in bundle")
cat("Plot key column:", plot_key, "\n")

# Subsample early
if (SUBSAMPLE > 0 && SUBSAMPLE < nrow(d)) {
  set.seed(42)
  d <- d[sample(.N, SUBSAMPLE)]
  cat("Subsampled to:", nrow(d), "plots\n")
}

# Build plot-level trait matrix W_plot by joining matched-pairs trees
# back, aggregating traits weighted by tree count per plot
cat("Loading pairs for plot composition ..."); flush.console()
trees <- as.data.table(readRDS(PAIRS_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done\n")

# Match plot keys (data.table needs explicit column access)
if (plot_key == "PLT_CN") {
  trees[, plot_match := as.numeric(bit64::as.character.integer64(PLT_CN_cond1))]
  d[, plot_match := PLT_CN]
} else if (plot_key == "plot_key") {
  trees[, plot_match := plot_key]
  d[, plot_match := plot_key]
} else if (plot_key == "PLT_CN_cond1") {
  trees[, plot_match := as.numeric(bit64::as.character.integer64(PLT_CN_cond1))]
  d[, plot_match := as.numeric(bit64::as.character.integer64(PLT_CN_cond1))]
} else {
  stop("Unsupported plot key: ", plot_key)
}

# Get species composition per plot
plot_spcd <- trees[, .(SPCD, plot_match)]
plot_spcd <- plot_spcd[plot_match %in% d$plot_match]
cat("Tree-rows in target plots:", nrow(plot_spcd), "\n")

# Trait columns (same set as base models)
trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
trait_cols <- intersect(trait_cols, names(traits))
traits_sub <- traits[, c("SPCD", trait_cols), with = FALSE]

# Z-score traits
for (col in trait_cols) {
  v <- traits_sub[[col]]
  med <- median(v, na.rm = TRUE); if (is.na(med)) med <- 0
  v[is.na(v)] <- med
  s <- sd(v); if (!is.finite(s) || s < 1e-12) s <- 1
  traits_sub[, (col) := (v - mean(v)) / s]
}

# Join traits to per-tree composition, then aggregate per plot
plot_spcd <- merge(plot_spcd, traits_sub, by = "SPCD", all.x = TRUE)
# Plot-mean trait values
W_plot_dt <- plot_spcd[, lapply(.SD, mean, na.rm = TRUE),
                       by = plot_match, .SDcols = trait_cols]
cat("W_plot rows:", nrow(W_plot_dt), "\n")

# Align W_plot to d by plot_match
d_w <- merge(d, W_plot_dt, by = "plot_match", all.x = TRUE)
# Fill any plots missing trait values with zeros
for (col in trait_cols) {
  v <- d_w[[col]]; v[!is.finite(v)] <- 0; d_w[[col]] <- v
}
d_w <- d_w[is.finite(residual)]
cat("Final rows for fit:", nrow(d_w), "\n")

W <- as.matrix(d_w[, trait_cols, with = FALSE])
P_trait <- ncol(W)

# Decay envelope
dstrb_col <- sprintf("dstrb_decay_%dyr", LAMBDA)
trt_col   <- sprintf("trt_decay_%dyr",   LAMBDA)

stan_data <- list(
  N_obs       = nrow(d_w),
  N_L1        = length(bundle$L1_levels),
  P_trait     = P_trait,
  residual    = d_w$residual,
  weight      = d_w$weight,
  is_plantation = as.numeric(d_w$is_plantation),
  d_fire      = as.numeric(d_w$had_fire_t1),
  d_insect    = as.numeric(d_w$had_insect_t1),
  d_disease   = as.numeric(d_w$had_disease_t1),
  d_wind      = as.numeric(d_w$had_wind_t1),
  d_harvest   = as.numeric(d_w$had_harvest_t1),
  dstrb_decay = d_w[[dstrb_col]],
  t_cutting   = as.numeric(d_w$had_cutting_t1),
  t_site_prep = as.numeric(d_w$had_site_prep_t1),
  trt_decay   = d_w[[trt_col]],
  L1_idx      = as.integer(d_w$L1_idx),
  W_plot      = W
)

cat("\nStan data: N_obs=", stan_data$N_obs, " N_L1=", stan_data$N_L1,
    " P_trait=", P_trait, "\n", sep = "")

mod <- cmdstan_model(STAN_FILE)
t_start <- Sys.time()
fit <- mod$sample(
  data = stan_data, chains = 2, parallel_chains = 2,
  iter_warmup = 300, iter_sampling = 300,
  seed = 42, adapt_delta = 0.9, refresh = 100
)
wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("\nWall:", round(wall_min, 1), "min\n\n")

out_name <- sprintf("%s_traitmed_plotlevel_lambda%d",
                    sub("_plotlevel$", "", bundle$model), LAMBDA)
fit$save_object(file.path(OUT_DIR, paste0(out_name, "_fit.rds")))

vars <- c("alpha_0",
          paste0("alpha_", c("plant","fire","insect","disease","wind",
                              "harvest","cutting","siteprep")),
          paste0("gamma_alpha_plant[",   seq_len(P_trait), "]"),
          paste0("gamma_alpha_fire[",    seq_len(P_trait), "]"),
          paste0("gamma_alpha_insect[",  seq_len(P_trait), "]"),
          paste0("gamma_alpha_disease[", seq_len(P_trait), "]"),
          paste0("gamma_alpha_wind[",    seq_len(P_trait), "]"),
          paste0("gamma_alpha_harvest[", seq_len(P_trait), "]"),
          paste0("gamma_alpha_cutting[", seq_len(P_trait), "]"),
          paste0("gamma_alpha_siteprep[",seq_len(P_trait), "]"),
          "sigma_L1", "sigma_resid")

summ <- fit$summary(variables = vars, "mean","median","sd","mad",
                    ~quantile(.x, c(0.05, 0.95)), "rhat","ess_bulk","ess_tail")
names(summ)[names(summ) %in% c("5%","95%")] <- c("q5","q95")
# split global and gamma
gamma_rows <- grepl("gamma_alpha_", summ$variable)
glob <- summ[!gamma_rows, ]
gam  <- summ[gamma_rows, ]
gam[, trait := rep(trait_cols, 8)]

data.table::fwrite(glob, file.path(OUT_DIR, paste0(out_name, "_global_summary.csv")))
data.table::fwrite(gam,  file.path(OUT_DIR, paste0(out_name, "_gamma_summary.csv")))

cat("=== Global alphas ===\n")
print(glob[grepl("^alpha", glob$variable), ])
cat("\n=== Significant gammas (90% CI excludes 0) ===\n")
sig <- gam[(q5 > 0) | (q95 < 0)]
print(sig)
cat("\nDone. Output in:", OUT_DIR, "\n")
