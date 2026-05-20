# =============================================================================
# Title: Score predictions on held-out species using only their trait vector
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: After 80_holdout_species_fit.R runs the fit on N-k species,
#              this script predicts annual diameter growth for the k held-out
#              species using ONLY their trait vector W (no species-level RE
#              available, because species-free architecture removed it).
#              Compares predicted distribution to the held-out observations.
#
# This is the strongest test of the trait-driven generalization claim. If
# prediction error on held-out species is comparable to in-sample error,
# the trait-substitution-for-species-RE story holds operationally.
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_holdout_predict.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/80b_holdout_species_predict.R \
#     #   --outdir=calibration/output/conus/dg/holdout_pilot_TIMESTAMP
#
# Outputs:
#   {outdir}/heldout_predictions.csv          per-tree predictions
#   {outdir}/heldout_species_summary.csv      RMSE / bias by held-out species
#   {outdir}/heldout_validation_panel.png/pdf observed vs predicted by species
# =============================================================================

library(tidyverse)
library(data.table)
library(posterior)
library(ggsci)
library(patchwork)

# --- CLI parsing -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
OUTDIR_REL <- get_arg("outdir", "calibration/output/conus/dg/holdout_pilot")
N_DRAWS    <- as.integer(get_arg("n_draws", "200"))

PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR <- file.path(PROJ_ROOT, OUTDIR_REL)
stopifnot(dir.exists(OUT_DIR))

# --- Load held-out metadata --------------------------------------------------
holdout_meta_path <- file.path(OUT_DIR, "holdout_species.csv")
stopifnot(file.exists(holdout_meta_path))
holdout <- read_csv(holdout_meta_path, show_col_types = FALSE)
cat("Held-out species:", paste(holdout$SPCD, collapse = ", "), "\n")

# --- Load fit and extract posterior draws of fixed-effect coefficients ------
FIT_PATH <- file.path(OUT_DIR, "dg_kuehne_b1_holdout_fit.rds")
stopifnot(file.exists(FIT_PATH))

cat("Loading fit ..."); flush.console()
fit <- readRDS(FIT_PATH)
cat(" done\n")

coef_names <- c(paste0("b", 0:8), paste0("gamma[", 1:8, "]"), "sigma")
all_draws <- fit$draws(variables = coef_names, format = "draws_matrix")
sel <- sample.int(nrow(all_draws), min(N_DRAWS, nrow(all_draws)))
draws <- all_draws[sel, , drop = FALSE]
rm(fit, all_draws); gc()

# --- Load held-out species observations -------------------------------------
DATA_FILE <- file.path(PROJ_ROOT,
  "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- file.path(PROJ_ROOT, "calibration/traits/species_traits.rds")

cat("Loading data ..."); flush.console()
dat <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done\n")

# Apply the same filters as the production driver, then keep only held-out
heldout_dat <- dat[
  SPCD %in% holdout$SPCD &
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(DBH2) &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 &
  is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(cspi) &
  is.finite(ba_x_rd) & is.finite(bal_x_rd)
]
cat("Held-out observations after filtering:", nrow(heldout_dat), "\n")

# Compute the same transformations used in the fit
CSPI_SHIFT <- 1.0
heldout_dat[, dg_obs_a := (DBH2 - DBH1) / YEARS]
heldout_dat[, sqrt_years := sqrt(YEARS)]
heldout_dat[, ln_dbh := log(DBH1)]
heldout_dat[, ln_cr_adj := log((CR1 + 0.2) / 1.2)]
heldout_dat[, ln_bal_sw_adj := log(BAL_SW1 + 0.01)]
heldout_dat[, ln_csi := log(pmax(cspi + CSPI_SHIFT, 0.01))]
heldout_dat <- heldout_dat[dg_obs_a > -0.5 & dg_obs_a < 5.0]
cat("After dg_obs_a range filter:", nrow(heldout_dat), "\n")

# Optional subsample for speed (keep all species represented)
set.seed(42)
SUBSAMPLE_PER_SPECIES <- 1000
heldout_dat <- heldout_dat[, .SD[sample(.N, min(.N, SUBSAMPLE_PER_SPECIES))],
                           by = SPCD]
cat("After per-species subsample:", nrow(heldout_dat), "\n")

# --- Build trait matrix for held-out species --------------------------------
trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
# Standardize using the same scaling that the fit used: load training species
# means and SDs from the meta file
META_PATH <- file.path(OUT_DIR, "dg_kuehne_b1_holdout_meta.rds")
meta <- readRDS(META_PATH)
training_traits <- traits[SPCD %in% meta$prep_meta$sp,
                          c("SPCD", trait_cols), with = FALSE]
# Same imputation strategy used in the fit driver
for (col in trait_cols) {
  med <- median(training_traits[[col]], na.rm = TRUE)
  training_traits[is.na(get(col)), (col) := med]
}
trait_means <- sapply(training_traits[, trait_cols, with = FALSE], mean)
trait_sds   <- sapply(training_traits[, trait_cols, with = FALSE], sd)

heldout_traits <- traits[SPCD %in% holdout$SPCD,
                         c("SPCD", trait_cols), with = FALSE]
# Impute with training medians
for (col in trait_cols) {
  med <- median(training_traits[[col]], na.rm = TRUE)
  heldout_traits[is.na(get(col)), (col) := med]
}
W_heldout <- as.matrix(heldout_traits[, trait_cols, with = FALSE])
for (j in seq_len(ncol(W_heldout))) {
  W_heldout[, j] <- (W_heldout[, j] - trait_means[j]) / trait_sds[j]
}
rownames(W_heldout) <- heldout_traits$SPCD

cat("Built W_heldout dim:", paste(dim(W_heldout), collapse = " x "), "\n")

# --- Predict ----------------------------------------------------------------
# eta = b0 + W*gamma + b1 ln(dbh) + b2 dbh + b3 ln_cr_adj
#     + b4 ln_bal_sw_adj + b5 bal_hw + b6 ln_csi + b7 ba_x_rd + b8 bal_x_rd
# trait_effect = W * gamma per species

gamma_cols <- paste0("gamma[", 1:8, "]")
gamma_draws <- draws[, gamma_cols]  # N_draws x 8
trait_effect_draws <- gamma_draws %*% t(W_heldout)  # N_draws x N_heldout_sp
colnames(trait_effect_draws) <- heldout_traits$SPCD

# Map each held-out observation to its trait_effect column
heldout_dat$sp_chr <- as.character(heldout_dat$SPCD)
te_lookup <- function(sp_chr_vec, draw_idx) {
  trait_effect_draws[draw_idx, sp_chr_vec]
}

# Predict observation by observation
b_draws <- draws[, paste0("b", 0:8)]   # N_draws x 9
sigma_draws <- draws[, "sigma"]

eta_mat <- matrix(0, nrow = nrow(draws), ncol = nrow(heldout_dat))
for (d in seq_len(nrow(draws))) {
  te_vec <- te_lookup(heldout_dat$sp_chr, d)
  eta_mat[d, ] <-
    b_draws[d, 1] + te_vec +
    b_draws[d, 2] * heldout_dat$ln_dbh +
    b_draws[d, 3] * heldout_dat$DBH1 +
    b_draws[d, 4] * heldout_dat$ln_cr_adj +
    b_draws[d, 5] * heldout_dat$ln_bal_sw_adj +
    b_draws[d, 6] * heldout_dat$BAL_HW1 +
    b_draws[d, 7] * heldout_dat$ln_csi +
    b_draws[d, 8] * heldout_dat$ba_x_rd +
    b_draws[d, 9] * heldout_dat$bal_x_rd
}

# Lognormal back-transform: predicted mean = exp(eta + sigma^2 / 2)
pred_mat <- exp(eta_mat + matrix(sigma_draws^2 / 2,
                                  nrow = nrow(eta_mat), ncol = ncol(eta_mat)))

heldout_dat$pred_mean  <- colMeans(pred_mat)
heldout_dat$pred_q5    <- apply(pred_mat, 2, quantile, 0.05)
heldout_dat$pred_q50   <- apply(pred_mat, 2, quantile, 0.50)
heldout_dat$pred_q95   <- apply(pred_mat, 2, quantile, 0.95)
heldout_dat$resid_mean <- heldout_dat$dg_obs_a - heldout_dat$pred_mean

# --- Save per-tree predictions ----------------------------------------------
write_csv(heldout_dat[, .(SPCD, DBH1, dg_obs_a, pred_mean, pred_q5, pred_q50,
                           pred_q95, resid_mean)],
          file.path(OUT_DIR, "heldout_predictions.csv"))

# --- Summary by species -----------------------------------------------------
sp_summary <- heldout_dat[, .(
  n         = .N,
  obs_mean  = mean(dg_obs_a),
  pred_mean = mean(pred_mean),
  bias      = mean(resid_mean),
  rmse      = sqrt(mean(resid_mean^2)),
  mae       = mean(abs(resid_mean)),
  r2        = 1 - sum(resid_mean^2) / sum((dg_obs_a - mean(dg_obs_a))^2),
  coverage_90 = mean(dg_obs_a >= pred_q5 & dg_obs_a <= pred_q95)
), by = SPCD]
sp_summary <- merge(sp_summary, holdout[, c("SPCD", "GENUS", "softwood")],
                    by = "SPCD")
sp_summary <- sp_summary[order(-n)]

write_csv(sp_summary, file.path(OUT_DIR, "heldout_species_summary.csv"))
cat("\n=== Held-out species prediction summary ===\n")
print(sp_summary)

# --- Headline metrics -------------------------------------------------------
overall <- tibble(
  n_species = nrow(sp_summary),
  n_obs     = sum(sp_summary$n),
  rmse_mean = mean(sp_summary$rmse),
  rmse_pooled = sqrt(mean(heldout_dat$resid_mean^2)),
  bias_pooled = mean(heldout_dat$resid_mean),
  coverage_90_pooled = mean(heldout_dat$dg_obs_a >= heldout_dat$pred_q5 &
                            heldout_dat$dg_obs_a <= heldout_dat$pred_q95),
  r2_pooled = 1 - sum(heldout_dat$resid_mean^2) /
                  sum((heldout_dat$dg_obs_a - mean(heldout_dat$dg_obs_a))^2)
)
cat("\n=== Overall held-out prediction performance ===\n")
print(overall)
write_csv(overall, file.path(OUT_DIR, "heldout_overall_summary.csv"))

# --- Figures ----------------------------------------------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold")
  )

# Panel A: observed vs predicted, colored by species
p_obs_pred <- ggplot(heldout_dat, aes(pred_mean, dg_obs_a,
                                       color = factor(SPCD))) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_color_aaas(name = "SPCD") +
  labs(x = "Predicted annual DG (cm/yr)",
       y = "Observed annual DG (cm/yr)",
       title = "A) Observed vs predicted on held-out species") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# Panel B: per-species RMSE bar plot
p_rmse <- ggplot(sp_summary, aes(reorder(factor(SPCD), -rmse), rmse,
                                  fill = factor(softwood))) +
  geom_col() +
  scale_fill_manual(values = c("0" = "#0072B5", "1" = "#E18727"),
                    name = "Softwood?", labels = c("Hardwood", "Softwood")) +
  labs(x = "Species (SPCD)", y = "RMSE (cm/yr)",
       title = "B) Per-species prediction RMSE") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel C: coverage of 90 percent interval per species
p_cov <- ggplot(sp_summary, aes(reorder(factor(SPCD), -coverage_90),
                                 coverage_90)) +
  geom_col(fill = "#3B4992") +
  geom_hline(yintercept = 0.90, color = "red", linetype = "dashed") +
  labs(x = "Species (SPCD)",
       y = "Coverage of 90% interval",
       title = "C) Predictive interval coverage by species",
       subtitle = "Target: 0.90 (red dashed)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylim(0, 1)

panel_fig <- (p_obs_pred) / (p_rmse | p_cov) +
  plot_annotation(
    title = "Held-out species validation: DG_Kuehne B1 species-free",
    subtitle = sprintf("%d species held out, %d observations, predicted via trait vector only.",
                       overall$n_species, overall$n_obs)
  )

ggsave(file.path(OUT_DIR, "heldout_validation_panel.png"), panel_fig,
       width = 22, height = 22, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "heldout_validation_panel.pdf"), panel_fig,
       width = 22, height = 22, units = "cm")

cat("\n=== Verdict ===\n")
if (overall$coverage_90_pooled > 0.85) {
  cat("Held-out coverage achieves at least 85 percent. Trait-based prediction\n")
  cat("transfers cleanly to species the model has never seen.\n")
} else if (overall$coverage_90_pooled > 0.70) {
  cat("Held-out coverage is moderate. Trait prediction partially works;\n")
  cat("expect some species to be poorly predicted.\n")
} else {
  cat("Held-out coverage is poor. Trait-based prediction does not generalize\n")
  cat("well to unseen species. Investigate whether trait scaling or the\n")
  cat("8 selected traits are sufficient.\n")
}

cat("\nDone. Outputs in", OUT_DIR, "\n")
