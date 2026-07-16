##=============================================================================
## 36_fit_htdbh_chapman.R
##
## Chapman-Richards height-diameter fit driver (Aaron's preferred H-D form).
## Supports both:
##   - species-specific (trait-informed species RE): stan/ht_dbh_chapman.stan
##   - species-free      (trait projection only):    stan/ht_dbh_chapman_speciesfree.stan
## via --stan_file. Same stan_data convention as 36_fit_htdbh_speciesfree.R
## (v2 pairs, BAL_SW1+BAL_HW1, sdi_additive1/SDImax_brms RD, climate_si,
## DF/lodgepole variety split). Adds a held-out 10% validation reporting RMSE
## and bias overall and BY DBH CLASS (the metric that exposed the Wykoff
## large-tree over-tall bias), plus a monotonicity check of the fitted curve.
##
## CLI: --stan_file --traits --outdir --outname --subsample=N --smoke
##      --min_sp=N --holdout_frac=0.1
##
## Author: A. Weiskittel + Claude   Date: 2026-06-27
##=============================================================================
suppressMessages({
  library(data.table); library(cmdstanr); library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

STAN_FILE <- get_arg("stan_file", "calibration/stan/ht_dbh_chapman_speciesfree.stan")
OUT_DIR   <- get_arg("outdir",    "calibration/output/conus/ht_dbh/chapman_sf")
OUT_NAME  <- get_arg("outname",   "htdbh_chapman")
SUBSAMPLE <- as.integer(get_arg("subsample", NA_character_))
SMOKE     <- has_flag("smoke")
HOLDOUT   <- as.numeric(get_arg("holdout_frac", "0.1"))
MIN_OBS_SPECIES <- as.integer(get_arg("min_sp", "5000"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== 36_fit_htdbh_chapman.R ==\nStan:", STAN_FILE, "\n\n")

DATA_FILE   <- "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
TRAITS_FILE <- get_arg("traits", "calibration/traits/species_traits_v2.rds")

cat("Loading data ..."); flush.console()
dat    <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done. Rows:", nrow(dat), "\n\n")

## ln(climate site index), shifted, coalesce NA -> median
if ("climate_si" %in% names(dat)) {
  med <- median(dat$climate_si, na.rm = TRUE)
  dat[!is.finite(climate_si), climate_si := med]
  dat[, ln_cspi_shift := log(pmax(climate_si, 0.1))]
} else { dat[, ln_cspi_shift := 0] }
dat[!is.finite(ln_cspi_shift), ln_cspi_shift := 0]

dat[, rd_ratio := sdi_additive1 / SDImax_brms]
dat[, bal_metric := (BAL_SW1 + BAL_HW1)]
dat[, ba_metric  := BA1 * 0.2296]
dat[, sqrt_ba    := sqrt(pmax(ba_metric, 0))]
dat[, ba_x_rd    := ba_metric  * rd_ratio]
dat[, bal_x_rd   := bal_metric * rd_ratio]

dat <- dat[
  TREESTATUS1 == 1 &
  is.finite(HT1) & HT1 >= 1.37 & HT1 < 85 &
  is.finite(DBH1) & DBH1 >= 2.54 & DBH1 < 250 &
  is.finite(BA1) & BA1 > 0 &
  is.finite(bal_metric) & bal_metric >= 0 &
  is.finite(rd_ratio) & rd_ratio > 0 & rd_ratio < 2 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  !is.na(FORTYPCD_cond1) & FORTYPCD_cond1 > 0
]
cat("After column filters:", nrow(dat), "rows\n")

## DF / lodgepole variety splits (only when v2 traits in use)
if (any(traits$SPCD == 2020L) && any(traits$SPCD == 2021L)) {
  dat[SPCD == 202L & as.character(EPA_L1_CODE) == "7", SPCD := 2020L]
  dat[SPCD == 202L, SPCD := 2021L]
}
if (any(traits$SPCD == 1080L) && any(traits$SPCD == 1081L)) {
  dat[SPCD == 108L & as.character(EPA_L1_CODE) == "7", SPCD := 1080L]
  dat[SPCD == 108L, SPCD := 1081L]
}

sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
dat <- dat[SPCD %in% sp_counts$SPCD]
cat("After species filter:", nrow(dat), "rows;", nrow(sp_counts), "species\n")

if (!is.na(SUBSAMPLE) && SUBSAMPLE < nrow(dat)) {
  set.seed(42); dat <- dat[sort(sample.int(nrow(dat), SUBSAMPLE))]
  cat("Subsampled to:", nrow(dat), "rows\n")
}

## held-out split (random; trees are the unit). seed fixed for reproducibility.
set.seed(7)
dat[, .ho := runif(.N) < HOLDOUT]
train <- dat[.ho == FALSE]; test <- dat[.ho == TRUE]
cat(sprintf("Train %d / Holdout %d\n\n", nrow(train), nrow(test)))

sp_levels <- sort(unique(train$SPCD))
L1_levels <- sort(unique(as.character(train$EPA_L1_CODE)))
L2_levels <- sort(unique(as.character(train$EPA_L2_CODE)))
L3_levels <- sort(unique(as.character(train$EPA_L3_CODE)))

mk_idx <- function(d) {
  d[, sp_idx := match(SPCD, sp_levels)]
  d[, L1_idx := match(as.character(EPA_L1_CODE), L1_levels)]
  d[, L2_idx := match(as.character(EPA_L2_CODE), L2_levels)]
  d[, L3_idx := match(as.character(EPA_L3_CODE), L3_levels)]
  d
}
train <- mk_idx(train)
test  <- mk_idx(test)
## keep only holdout rows whose species/ecoregion appear in train
test <- test[!is.na(sp_idx) & !is.na(L1_idx) & !is.na(L2_idx) & !is.na(L3_idx)]

trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
trait_cols <- intersect(trait_cols, names(traits))
traits_sub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
W <- as.matrix(traits_sub[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W))) {
  na <- is.na(W[, j]); if (any(na)) W[na, j] <- median(W[!na, j], na.rm = TRUE)
  W[, j] <- (W[, j] - mean(W[, j])) / sd(W[, j])
}

stan_data <- list(
  N_obs = nrow(train), N_sp = length(sp_levels),
  N_L1 = length(L1_levels), N_L2 = length(L2_levels), N_L3 = length(L3_levels),
  P_trait = ncol(W),
  ht_obs = train$HT1, dbh = train$DBH1,
  bal = train$bal_metric, sqrt_ba = train$sqrt_ba,
  ln_cspi_shift = train$ln_cspi_shift,
  ba_x_rd = train$ba_x_rd, bal_x_rd = train$bal_x_rd,
  sp_idx = train$sp_idx, L1_idx = train$L1_idx,
  L2_idx = train$L2_idx, L3_idx = train$L3_idx, W = W
)
cat("=== Stan data ready ===  N_obs =", stan_data$N_obs,
    " N_sp =", stan_data$N_sp, "\n\n")

mod <- cmdstan_model(STAN_FILE)
if (SMOKE) { iw <- 75; is_ <- 75; ch <- 2 } else { iw <- 1000; is_ <- 1000; ch <- 4 }

t0 <- Sys.time()
fit <- mod$sample(data = stan_data, chains = ch, parallel_chains = ch,
                  iter_warmup = max(iw, 1500L), iter_sampling = is_, seed = 42,
                  adapt_delta = 0.90, max_treedepth = 12, refresh = 100,
                  init = function() list(a0 = 3.2, b_rate = 0.04, c_shape = 1.0,
                    a_bal = 0, a_ba = 0, a_cspi = 0.2, a_bard = 0, a_blrd = 0,
                    gamma = rep(0, ncol(W)), s0 = 0.7, s1 = 0.5,
                    sigma_L1 = 0.1, sigma_L2 = 0.1, sigma_L3 = 0.1))
wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat("\nWall:", round(wall_min, 1), "min\n\n")

svars <- fit$metadata()$stan_variables
want <- c("a0","b_rate","c_shape","a_bal","a_ba","a_cspi","a_bard","a_blrd",
          paste0("gamma[", seq_len(ncol(W)), "]"),
          "sigma_sp","sigma_L1","sigma_L2","sigma_L3","s0","s1")
want <- want[ sub("\\[.*$","",want) %in% svars ]
summ <- fit$summary(variables = want, "mean","median","sd",
                    ~quantile(.x, c(0.05,0.95)), "rhat","ess_bulk","ess_tail")
names(summ)[names(summ) %in% c("5%","95%")] <- c("q5","q95")
fwrite(summ, file.path(OUT_DIR, paste0(OUT_NAME, "_summary.csv")))
cat("=== key coefficients ===\n"); print(summ[grepl("^a0|^b_rate|^c_shape|^a_|^s0|^s1", summ$variable), ])

## ---- posterior-mean point predictions for held-out validation ----
pm <- function(v) { if (v %in% want || v %in% svars) mean(fit$draws(v, format="draws_matrix")) else NA_real_ }
a0 <- pm("a0"); b_rate <- pm("b_rate"); c_shape <- pm("c_shape")
a_bal <- pm("a_bal"); a_ba <- pm("a_ba"); a_cspi <- pm("a_cspi")
a_bard <- pm("a_bard"); a_blrd <- pm("a_blrd")
gam <- colMeans(fit$draws(paste0("gamma[", seq_len(ncol(W)), "]"), format="draws_matrix"))
spvec_name <- if ("trait_effect" %in% svars) "trait_effect" else "z_sp"
spvec <- colMeans(fit$draws(spvec_name, format="draws_matrix"))   # length N_sp
zL1 <- colMeans(fit$draws("z_L1", format="draws_matrix"))
zL2 <- colMeans(fit$draws("z_L2", format="draws_matrix"))
zL3 <- colMeans(fit$draws("z_L3", format="draws_matrix"))

pred_ht <- function(d) {
  logA <- a0 + spvec[d$sp_idx] + zL1[d$L1_idx] + zL2[d$L2_idx] + zL3[d$L3_idx] +
          a_bal*d$bal_metric + a_ba*d$sqrt_ba + a_cspi*d$ln_cspi_shift +
          a_bard*d$ba_x_rd + a_blrd*d$bal_x_rd
  A <- exp(logA)
  term <- pmax(1 - exp(-b_rate * d$DBH1), 1e-8)
  1.37 + A * term^c_shape
}
test[, ht_hat := pred_ht(test)]
test[, resid := ht_hat - HT1]
val <- function(o,p){ r<-p-o; data.table(n=length(r), bias=mean(r),
  bias_pct=100*mean(r)/mean(o), rmse=sqrt(mean(r^2)),
  rmse_pct=100*sqrt(mean(r^2))/mean(o)) }
overall <- val(test$HT1, test$ht_hat)
cat("\n=== HELD-OUT validation (overall) ===\n"); print(overall)

## bias by DBH class (cm bins mirroring the Jun-17 in/inch classes)
test[, dbh_class := cut(DBH1, breaks=c(0,7.6,12.7,22.9,33,48.3,1000),
     labels=c("1-3in","3-5in","5-9in","9-13in","13-19in","19-40in"),
     right=FALSE)]
by_class <- test[, val(HT1, ht_hat), by = dbh_class][order(dbh_class)]
cat("\n=== HELD-OUT bias by DBH class (pred - obs) ===\n"); print(by_class)
fwrite(by_class, file.path(OUT_DIR, paste0(OUT_NAME, "_holdout_bias_by_dbhclass.csv")))
fwrite(overall,  file.path(OUT_DIR, paste0(OUT_NAME, "_holdout_overall.csv")))

## monotonicity check: predicted curve over DBH grid for the median species/eco
grid <- data.table(DBH1 = seq(2.5, 120, by = 2.5))
grid[, `:=`(sp_idx = which.max(tabulate(test$sp_idx)),
            L1_idx = which.max(tabulate(test$L1_idx)),
            L2_idx = which.max(tabulate(test$L2_idx)),
            L3_idx = which.max(tabulate(test$L3_idx)),
            bal_metric = median(test$bal_metric), sqrt_ba = median(test$sqrt_ba),
            ln_cspi_shift = median(test$ln_cspi_shift),
            ba_x_rd = median(test$ba_x_rd), bal_x_rd = median(test$bal_x_rd))]
grid[, ht_hat := pred_ht(grid)]
mono <- all(diff(grid$ht_hat) >= -1e-6)
cat("\nMonotone non-decreasing in DBH over [2.5,120] cm:", mono,
    "  asymptote ~", round(max(grid$ht_hat),1), "m\n")
fwrite(grid[, .(DBH1, ht_hat)], file.path(OUT_DIR, paste0(OUT_NAME, "_curve_grid.csv")))

saveRDS(list(form="chapman_richards", stan_file=STAN_FILE,
             species_free=("trait_effect" %in% svars),
             trait_cols=trait_cols, sp_levels=sp_levels,
             L1_levels=L1_levels, L2_levels=L2_levels, L3_levels=L3_levels,
             summary=summ, overall=overall, by_class=by_class,
             monotone=mono, n_train=stan_data$N_obs, n_holdout=nrow(test),
             wall_min=wall_min),
        file.path(OUT_DIR, paste0(OUT_NAME, "_meta.rds")))
if (!SMOKE && !has_flag("nofit")) fit$save_object(file.path(OUT_DIR, paste0(OUT_NAME, "_fit.rds")))
cat("\nDone.\n")
