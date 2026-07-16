##=============================================================================
## 36_fit_ingrowth_species_composition.R
##
## Stage 4 of the four-stage ingrowth model: species composition.
## Multinomial with trait-driven covariate effects. Reads the per-recruit
## dataset built by 01c_build_ingrowth_per_recruit.R.
##
## CLI:
##   --stan_file=PATH (default ingrowth_species_composition_v1.stan)
##   --recruit_data=PATH (default data/ingrowth_per_recruit.rds)
##   --traits=PATH (default species_traits_v2.rds)
##   --outdir=PATH
##   --outname=NAME
##   --subsample=N (plot subsample)
##   --top_n_sp=N (default 50; rest lumped as OTHER)
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

STAN_FILE     <- get_arg("stan_file", "calibration/stan/ingrowth_species_composition_v1.stan")
RECRUIT_DATA  <- get_arg("recruit_data", "calibration/data/ingrowth_per_recruit.rds")
TRAITS_FILE   <- get_arg("traits",
                         if (file.exists("calibration/traits/species_traits_v2.rds"))
                           "calibration/traits/species_traits_v2.rds"
                         else "calibration/traits/species_traits.rds")
OUT_DIR       <- get_arg("outdir",  "calibration/output/conus/ingrowth/composition")
OUT_NAME      <- get_arg("outname", "ingrowth_compos_v1")
SUBSAMPLE     <- as.integer(get_arg("subsample", NA_character_))
TOP_N_SP      <- as.integer(get_arg("top_n_sp", "50"))
SMOKE         <- has_flag("smoke")
MAX_TREEDEPTH <- as.integer(get_arg("max_treedepth", "10"))
ADAPT_DELTA   <- as.numeric(get_arg("adapt_delta", "0.9"))
ITER_WARMUP   <- as.integer(get_arg("iter_warmup",   if (has_flag("smoke")) "50"   else "1000"))
ITER_SAMPLING <- as.integer(get_arg("iter_sampling", if (has_flag("smoke")) "50"   else "1000"))
CHAINS        <- as.integer(get_arg("chains",        if (has_flag("smoke")) "2"    else "4"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 36_fit_ingrowth_species_composition.R ==\n")
cat("Stan:    ", STAN_FILE, "\n")
cat("Recruits:", RECRUIT_DATA, "\n")
cat("Traits:  ", TRAITS_FILE, "\n")
cat("Top N sp:", TOP_N_SP, "\n\n")

## 1. Load per-recruit data ------------------------------------------------
recruits <- as.data.table(readRDS(RECRUIT_DATA))
recruits[, PLT_CN := as.character(PLT_CN)]
cat(sprintf("Per-recruit rows: %s (plots = %s, species = %s)\n",
            format(nrow(recruits), big.mark = ","),
            format(uniqueN(recruits$PLT_CN), big.mark = ","),
            format(uniqueN(recruits$SPCD), big.mark = ",")))

## 2. Load plot covariates (from remeasurement pairs) ---------------------
PAIRS_FILE <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
pairs <- as.data.table(readRDS(PAIRS_FILE))
if (!"PLT_CN" %in% names(pairs) && "PLT_CN_cond1" %in% names(pairs)) {
  setnames(pairs, "PLT_CN_cond1", "PLT_CN")
}
pairs[, PLT_CN := as.character(PLT_CN)]
plot_cov <- pairs[, .(
    EPA_L1_CODE = first(EPA_L1_CODE),
    cspi      = first(cspi),
    clim_pca1 = first(clim_pca1),
    BA        = sum(BA1[TREESTATUS1 == 1], na.rm = TRUE),
    BAL_mean  = mean(BAL1[TREESTATUS1 == 1], na.rm = TRUE),
    RD        = first(SDI1) / first(SDImax_brms),
    ht40      = if ("HT40_DOM_t1" %in% names(.SD)) first(HT40_DOM_t1) else NA_real_
  ), by = PLT_CN]

## 3. Select top-N species + OTHER ----------------------------------------
sp_freq <- recruits[, .(n = .N), by = SPCD][order(-n)]
top_spcd <- sp_freq[1:min(TOP_N_SP, nrow(sp_freq)), SPCD]
recruits[, sp_label := ifelse(SPCD %in% top_spcd, SPCD, -1L)]   # -1 = OTHER
sp_levels <- c(sort(top_spcd), -1L)
n_sp <- length(sp_levels)
cat(sprintf("Species pool: top %d species + OTHER (n_sp = %d)\n",
            length(top_spcd), n_sp))

## 4. Apply DF + LP variety split to RECRUIT species too -----------------
traits <- as.data.table(readRDS(TRAITS_FILE))
recruits <- recruits[plot_cov[, .(PLT_CN, EPA_L1_CODE)], on = "PLT_CN", nomatch = 0]
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  pre <- sum(recruits$sp_label == 202L)
  recruits[sp_label == 202L & as.character(EPA_L1_CODE) == "7", sp_label := 2020L]
  recruits[sp_label == 202L, sp_label := 2021L]
  if (2020L %in% recruits$sp_label || 2021L %in% recruits$sp_label) {
    sp_levels <- unique(c(sp_levels, intersect(c(2020L, 2021L), recruits$sp_label)))
    sp_levels <- sp_levels[sp_levels != 202L]
  }
  cat(sprintf("DF recruit split: %d -> coastal + rocky\n", pre))
}
if (any(traits$SPCD == 1080L) && any(traits$SPCD == 1081L)) {
  pre <- sum(recruits$sp_label == 108L)
  recruits[sp_label == 108L & as.character(EPA_L1_CODE) == "7", sp_label := 1080L]
  recruits[sp_label == 108L, sp_label := 1081L]
  if (1080L %in% recruits$sp_label || 1081L %in% recruits$sp_label) {
    sp_levels <- unique(c(sp_levels, intersect(c(1080L, 1081L), recruits$sp_label)))
    sp_levels <- sp_levels[sp_levels != 108L]
  }
  cat(sprintf("LP recruit split: %d -> shore + rocky\n", pre))
}
sp_levels <- sort(sp_levels)
n_sp <- length(sp_levels)
recruits[, sp_idx := match(sp_label, sp_levels)]

## 5. Build (plot, species) count matrix -----------------------------------
plot_sp <- recruits[, .N, by = .(PLT_CN, sp_idx)]
y_long <- plot_sp[CJ(PLT_CN = unique(PLT_CN), sp_idx = seq_len(n_sp), unique = TRUE),
                  on = c("PLT_CN", "sp_idx")]
y_long[is.na(N), N := 0]
y_wide <- dcast(y_long, PLT_CN ~ sp_idx, value.var = "N", fill = 0)
plot_ids <- y_wide$PLT_CN
y_count <- as.matrix(y_wide[, -1, with = FALSE])
total_recruits <- rowSums(y_count)
cat(sprintf("Plot x species count matrix: %d plots x %d species\n",
            nrow(y_count), ncol(y_count)))

## 6. Join with plot covariates, filter ----------------------------------
plot_data <- data.table(PLT_CN = plot_ids, total_recruits = total_recruits)
plot_data <- plot_cov[plot_data, on = "PLT_CN", nomatch = 0]
plot_data <- plot_data[is.finite(cspi) & is.finite(BA) & BA >= 0 &
                       is.finite(RD) & RD > 0 & RD < 3 &
                       !is.na(EPA_L1_CODE) & EPA_L1_CODE != "" &
                       total_recruits > 0]
plot_data[is.na(BAL_mean), BAL_mean := 0]
plot_data[is.na(ht40), ht40 := 5.0]
plot_data[is.na(clim_pca1), clim_pca1 := 0]
keep <- match(plot_data$PLT_CN, plot_ids)
y_count <- y_count[keep, , drop = FALSE]
cat("After filter:", nrow(plot_data), "plots\n")

## 7. Subsample plots ------------------------------------------------------
if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(plot_data)) {
  set.seed(42)
  idx <- sample(nrow(plot_data), SUBSAMPLE)
  plot_data <- plot_data[idx]
  y_count <- y_count[idx, , drop = FALSE]
  cat("Subsampled to:", nrow(plot_data), "plots\n")
}

## 8. Build standardized covariate matrix + trait matrix ------------------
plot_data[, ln_ba   := log(BA * 0.2296 + 1.0)]
plot_data[, ln_bal  := log(BAL_mean * 0.2296 + 1.0)]
plot_data[, ln_csi  := log(pmax(cspi + 1.0, 0.01))]
cov_cols <- c("ln_ba", "ln_bal", "RD", "ht40", "ln_csi", "clim_pca1")
X <- as.matrix(plot_data[, ..cov_cols])
# Drop columns with zero variance (e.g., ht40 when HT40_DOM_t1 missing)
keep_cols <- vapply(seq_len(ncol(X)),
                    function(j) is.finite(sd(X[, j], na.rm = TRUE)) && sd(X[, j], na.rm = TRUE) > 0,
                    logical(1))
if (!all(keep_cols)) {
  cat("Dropping zero-variance covariates:",
      paste(cov_cols[!keep_cols], collapse = ", "), "\n")
  X <- X[, keep_cols, drop = FALSE]
  cov_cols <- cov_cols[keep_cols]
}
scale_mean <- numeric(ncol(X)); scale_sd <- numeric(ncol(X))
for (j in seq_len(ncol(X))) {
  m <- mean(X[, j], na.rm = TRUE)
  s_ <- sd(X[, j], na.rm = TRUE)
  scale_mean[j] <- m; scale_sd[j] <- s_
  X[, j] <- (X[, j] - m) / s_
  X[is.na(X[, j]), j] <- 0  # any remaining NAs replaced with 0 (column mean post-standardization)
}
names(scale_mean) <- cov_cols; names(scale_sd) <- cov_cols  # saved in meta for exact out-of-sample prediction

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
# Trait for each recruit species (use medians for OTHER = -1)
W_rows <- list()
for (s in sp_levels) {
  if (s == -1L) {
    tr <- traits[, lapply(.SD, median, na.rm = TRUE), .SDcols = trait_cols]
  } else {
    tr <- traits[SPCD == s, ..trait_cols]
    if (nrow(tr) == 0) tr <- traits[, lapply(.SD, median, na.rm = TRUE), .SDcols = trait_cols]
  }
  W_rows[[length(W_rows) + 1]] <- as.numeric(tr[1, ])
}
W_sp <- do.call(rbind, W_rows)
# Impute and standardize
for (j in seq_len(ncol(W_sp))) {
  na <- is.na(W_sp[, j])
  if (any(na)) W_sp[na, j] <- median(W_sp[!na, j], na.rm = TRUE)
  W_sp[, j] <- (W_sp[, j] - mean(W_sp[, j])) / sd(W_sp[, j])
}
cat("X dim:", paste(dim(X), collapse = "x"),
    " W_sp dim:", paste(dim(W_sp), collapse = "x"), "\n")

## 9. Stan data ------------------------------------------------------------
stan_data <- list(
  N_plots = nrow(plot_data),
  N_sp    = n_sp,
  P_cov   = ncol(X),
  P_trait = ncol(W_sp),
  total_recruits = as.integer(plot_data$total_recruits),
  y_count = y_count,
  X = X,
  W_sp = W_sp
)

## 10. Compile + sample ---------------------------------------------------
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

fit$save_object(file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds")))
vars <- c(paste0("alpha_0[", 1:n_sp, "]"),
          paste0("b[", 1:ncol(X), "]"),
          paste0("gamma_int[", 1:ncol(W_sp), "]"),
          ## gamma_cov is the trait x covariate matrix that drives the composition gradient;
          ## it MUST be in the saved summary or the fit cannot predict out of sample.
          paste0("gamma_cov[", rep(1:ncol(W_sp), times=ncol(X)), ",",
                                rep(1:ncol(X), each=ncol(W_sp)), "]"),
          "sigma_alpha")
vars <- intersect(vars, variables(fit$draws()))
summ <- fit$summary(variables = vars, "mean","median","sd",
                    ~quantile(.x, c(0.05, 0.95), na.rm = TRUE),
                    "rhat","ess_bulk","ess_tail")
names(summ)[names(summ) %in% c("5%","95%")] <- c("q5","q95")
data.table::fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv")))
saveRDS(list(form="multinomial_v1", sp_levels=sp_levels, trait_cols=trait_cols,
             cov_cols=cov_cols, scale_mean=scale_mean, scale_sd=scale_sd,
             stan_file=STAN_FILE, summary=summ,
             n_plots=stan_data$N_plots, n_sp=n_sp, wall_min=wm),
        file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds")))
cat("Done.\n")
