# =============================================================================
# Title: Residual diagnostics for a fitted FVS-CONUS base model
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: Produces observed vs predicted, residuals vs predicted, QQ,
#              residual histogram, residuals by species, residuals by EPA L1,
#              and residuals by DBH size class. Designed to run on one banked
#              fit at a time; loops over banked fits via a wrapper.
#
# Memory budget: ~64 GB. Loads one fit (5 to 14 GB) and the matching
# prepared data, computes posterior mean predictions, then drops the fit.
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_residuals.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/20_residual_diagnostics.R --model=hg --variant=b2
#
# CLI flags:
#   --model={hg,htdbh,cr,hcb,mort,dg_kue}
#   --variant={b1,b2}
#   --subsample=N      compute residuals on a random N-row subsample (default: 50000)
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggsci)
library(patchwork)
library(posterior)

# --- CLI parsing ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

MODEL     <- get_arg("model",     "hg")
VARIANT   <- get_arg("variant",   "b2")
SUBSAMPLE <- as.integer(get_arg("subsample", "50000"))

cat(sprintf("Residual diagnostics: model=%s variant=%s subsample=%d\n",
            MODEL, VARIANT, SUBSAMPLE))

# --- Paths and lookups ------------------------------------------------------
PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR   <- file.path(PROJ_ROOT, "calibration/output/evaluation/residuals",
                       paste(MODEL, VARIANT, sep = "_"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Map model + variant to fit.rds, data.rds, and response transform.
fit_paths <- tibble(
  model   = c("hg",       "hg",       "dg_kue",   "dg_kue",
              "htdbh",    "cr",       "hcb",      "mort"),
  variant = c("b2",       "b1",       "b2",       "b1",
              "b2",       "b2",       "b2",       "b2"),
  fit_path = c(
    "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds",
    "calibration/output/conus/hg/speciesfree_pilot/hg_organon_fixedK_cspi_traits1_fit.rds",
    "calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds",
    "calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds",
    "calibration/output/conus/ht_dbh/htdbh_wykoff_lognormal_cspi_traits1_fit.rds",
    "calibration/output/conus/crown_recession/cr_recession_cspi_traits1_fit.rds",
    "calibration/output/conus/hcb/hcb_organon_cspi_traits1_fit.rds",
    "calibration/output/conus/mortality/mort_logit_simple_cspi_traits1_fit.rds"
  ),
  response_label = c(
    "Annual height growth (m yr-1, log scale)",
    "Annual height growth (m yr-1, log scale)",
    "Annual diameter growth (cm yr-1, log scale)",
    "Annual diameter growth (cm yr-1, log scale)",
    "Total height (m, log scale)",
    "Crown recession (Beta link)",
    "Height to crown base (Beta link)",
    "Annual mortality probability (logit)"
  )
)

lookup <- fit_paths %>% filter(model == MODEL, variant == VARIANT)
stopifnot(nrow(lookup) == 1)
FIT_PATH <- file.path(PROJ_ROOT, lookup$fit_path)
stopifnot(file.exists(FIT_PATH))

# --- Load fit + meta --------------------------------------------------------
cat("Loading fit:", FIT_PATH, "\n"); flush.console()
fit <- readRDS(FIT_PATH)

META_PATH <- sub("_fit\\.rds$", "_meta.rds", FIT_PATH)
meta <- if (file.exists(META_PATH)) readRDS(META_PATH) else NULL

# --- Compute posterior mean predictions on a subsample ---------------------
#
# Strategy: use the posterior mean of mu_a (or eta_gq) stored in
# generated_quantities, if present. Otherwise fall back to draws of the
# linear predictor. Subsample observations to keep memory under control.

cat("Pulling posterior summaries of mu_a or eta_gq ...\n"); flush.console()

draws_mu <- tryCatch({
  fit$draws(variables = "mu_a", format = "draws_matrix")
}, error = function(e) {
  fit$draws(variables = "eta_gq", format = "draws_matrix")
})

n_obs <- ncol(draws_mu)
cat("N_obs in fit:", n_obs, "\n")

set.seed(42)
idx <- sort(sample.int(n_obs, min(SUBSAMPLE, n_obs)))

pred_mean <- colMeans(draws_mu[, idx, drop = FALSE])
pred_q5   <- apply(draws_mu[, idx, drop = FALSE], 2, quantile, 0.05)
pred_q95  <- apply(draws_mu[, idx, drop = FALSE], 2, quantile, 0.95)

# Read the prepared response and grouping vectors from the data list saved
# alongside the fit. If not saved, reconstruct from the source data.
stan_data_path <- sub("_fit\\.rds$", "_standata.rds", FIT_PATH)
if (file.exists(stan_data_path)) {
  sd_list <- readRDS(stan_data_path)
  obs_vec <- sd_list[["dg_obs_a"]]
  if (is.null(obs_vec)) obs_vec <- sd_list[["y"]]
  sp_idx  <- sd_list[["sp_idx"]]
  L1_idx  <- sd_list[["L1_idx"]]
  dbh_vec <- sd_list[["dbh"]]
  cat("Loaded standata from", stan_data_path, "\n")
} else {
  cat("No standata.rds alongside fit. Predicted vs observed comparison will\n")
  cat("be limited to predicted distribution only. Save standata next time.\n")
  obs_vec <- NULL
  sp_idx <- NULL
  L1_idx <- NULL
  dbh_vec <- NULL
}

# --- Build residual data.table ----------------------------------------------

resid_df <- tibble(
  obs_idx  = idx,
  pred     = pred_mean,
  pred_q5  = pred_q5,
  pred_q95 = pred_q95
)
if (!is.null(obs_vec)) {
  resid_df$obs   <- obs_vec[idx]
  resid_df$resid <- resid_df$obs - resid_df$pred
}
if (!is.null(sp_idx)) resid_df$sp_idx <- sp_idx[idx]
if (!is.null(L1_idx)) resid_df$L1_idx <- L1_idx[idx]
if (!is.null(dbh_vec)) {
  resid_df$dbh <- dbh_vec[idx]
  resid_df$dbh_class <- cut(resid_df$dbh,
    breaks = c(0, 10, 20, 30, 40, 60, 100, Inf),
    labels = c("<10 cm", "10-20 cm", "20-30 cm", "30-40 cm",
               "40-60 cm", "60-100 cm", ">100 cm"))
}

write_csv(resid_df, file.path(OUT_DIR, "residuals_sample.csv"))

# --- Drop fit to free memory ------------------------------------------------
rm(fit, draws_mu)
gc()

# --- Figures ----------------------------------------------------------------

theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold")
  )

if (!"obs" %in% names(resid_df)) {
  cat("No observations available; skipping residual diagnostics.\n")
  quit(save = "no")
}

# Panel A: observed vs predicted
p_obs_pred <- ggplot(resid_df, aes(pred, obs)) +
  geom_point(alpha = 0.25, size = 0.7, color = "#3B4992") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "#EE0000",
              linewidth = 0.5, span = 0.5) +
  labs(x = "Predicted", y = "Observed",
       title = paste0("A) Observed vs predicted (", MODEL, " ", VARIANT, ")"),
       subtitle = lookup$response_label) +
  theme_pub

# Panel B: residuals vs predicted
p_resid_pred <- ggplot(resid_df, aes(pred, resid)) +
  geom_point(alpha = 0.25, size = 0.7, color = "#3B4992") +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "#EE0000",
              linewidth = 0.5, span = 0.5) +
  labs(x = "Predicted", y = "Residual (obs - pred)",
       title = "B) Residuals vs predicted") +
  theme_pub

# Panel C: QQ
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq(alpha = 0.4, size = 0.7, color = "#3B4992") +
  stat_qq_line(color = "red", linewidth = 0.5) +
  labs(x = "Theoretical quantile", y = "Sample quantile",
       title = "C) Normal Q-Q") +
  theme_pub

# Panel D: histogram of residuals
p_hist <- ggplot(resid_df, aes(resid)) +
  geom_histogram(bins = 40, fill = "#3B4992", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "red", linewidth = 0.5) +
  labs(x = "Residual", y = "Count",
       title = "D) Residual distribution") +
  theme_pub

# Panel E: residuals by EPA L1 (boxplot, ordered by median)
if ("L1_idx" %in% names(resid_df)) {
  L1_order <- resid_df %>%
    group_by(L1_idx) %>%
    summarise(med = median(resid), .groups = "drop") %>%
    arrange(med) %>%
    pull(L1_idx)
  p_l1 <- resid_df %>%
    mutate(L1_idx = factor(L1_idx, levels = L1_order)) %>%
    ggplot(aes(L1_idx, resid)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.4) +
    geom_boxplot(fill = "#EE000022", outlier.size = 0.4, outlier.alpha = 0.3,
                 linewidth = 0.3) +
    labs(x = "EPA L1 ecoregion (sorted by median)", y = "Residual",
         title = "E) Residuals by EPA L1 ecoregion") +
    theme_pub
} else {
  p_l1 <- ggplot() + theme_void() + labs(title = "E) (L1 not available)")
}

# Panel F: residuals by DBH size class
if ("dbh_class" %in% names(resid_df)) {
  p_dbh <- ggplot(resid_df, aes(dbh_class, resid)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.4) +
    geom_boxplot(fill = "#0072B5", outlier.size = 0.4, outlier.alpha = 0.3,
                 linewidth = 0.3, alpha = 0.4) +
    labs(x = "DBH size class", y = "Residual",
         title = "F) Residuals by size class") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p_dbh <- ggplot() + theme_void() + labs(title = "F) (DBH not available)")
}

# --- Combine and save -------------------------------------------------------
panel_fig <- (p_obs_pred + p_resid_pred) /
             (p_qq + p_hist) /
             (p_l1 + p_dbh) +
  plot_annotation(
    title = sprintf("Residual diagnostics: %s %s", MODEL, VARIANT),
    subtitle = sprintf("Posterior mean predictions on %d random observations",
                       nrow(resid_df))
  )

fig_png <- file.path(OUT_DIR, "residual_diagnostics.png")
fig_pdf <- file.path(OUT_DIR, "residual_diagnostics.pdf")
ggsave(fig_png, panel_fig, width = 22, height = 28, units = "cm",
       dpi = 300, bg = "white")
ggsave(fig_pdf, panel_fig, width = 22, height = 28, units = "cm")
cat("Wrote", fig_png, "\n")
cat("Wrote", fig_pdf, "\n")

# --- Summary fit statistics -------------------------------------------------
fit_stats <- tibble(
  model   = MODEL,
  variant = VARIANT,
  n_resid = nrow(resid_df),
  rmse    = sqrt(mean(resid_df$resid^2, na.rm = TRUE)),
  bias    = mean(resid_df$resid, na.rm = TRUE),
  mae     = mean(abs(resid_df$resid), na.rm = TRUE),
  r2_classic = 1 - sum(resid_df$resid^2, na.rm = TRUE) /
                   sum((resid_df$obs - mean(resid_df$obs, na.rm = TRUE))^2, na.rm = TRUE)
) %>% mutate(across(where(is.numeric), ~round(., 4)))

stats_path <- file.path(OUT_DIR, "fit_statistics.csv")
write_csv(fit_stats, stats_path)
cat("Wrote", stats_path, "\n")

cat("\n=== Fit statistics ===\n")
print(fit_stats)
cat("\nDone.\n")
