# =============================================================================
# Title: Extrapolation diagnostics: behavior at edge-of-range and beyond
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: Generates a synthetic Cartesian grid of covariates spanning
#              within-range and out-of-range values, predicts using the fit,
#              and flags predictions that violate biological constraints
#              (negative growth where positive expected, mortality outside
#              [0, 1], crown ratios outside [0, 1], growth values an order of
#              magnitude beyond what is physically plausible).
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_extrap.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/70_extrapolation_diagnostics.R --model=dg_kue --variant=b2
#
# Outputs:
#   calibration/output/evaluation/extrapolation/{model}_{variant}/
#     extrapolation_summary.csv     - rate of biological violations per region
#     extrapolation_heatmap.png/pdf - heatmap showing "danger zones"
#     edge_behavior.png/pdf         - predictions at quantile edges
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
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
MODEL    <- get_arg("model", "dg_kue")
VARIANT  <- get_arg("variant", "b2")
N_DRAWS  <- as.integer(get_arg("n_draws", "100"))

cat(sprintf("Extrapolation diagnostics: model=%s variant=%s n_draws=%d\n",
            MODEL, VARIANT, N_DRAWS))

# --- Paths -------------------------------------------------------------------
PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR <- file.path(PROJ_ROOT, "calibration/output/evaluation/extrapolation",
                     paste(MODEL, VARIANT, sep = "_"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Model fit + plausibility bounds ----------------------------------------

model_spec <- list(
  hg = list(
    fit_b2 = "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds",
    fit_b1 = "calibration/output/conus/hg/speciesfree_pilot/hg_organon_fixedK_cspi_traits1_fit.rds",
    response_max = 2.0,   # annual height growth in meters, biologically plausible upper
    response_min = 0.0,
    response_unit = "m yr-1",
    cov_specs = list(
      dbh        = list(in_range = c(2.5, 80),  out_range = c(0.5, 150)),
      ln_cr_adj  = list(in_range = c(-1.5, 0),  out_range = c(-3, 0)),
      bal        = list(in_range = c(0, 60),    out_range = c(0, 100)),
      cspi       = list(in_range = c(-2, 4),    out_range = c(-4, 6))
    )
  ),
  dg_kue = list(
    fit_b2 = "calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds",
    fit_b1 = "calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds",
    response_max = 5.0,   # cm/yr, biologically plausible upper
    response_min = 0.0,
    response_unit = "cm yr-1",
    cov_specs = list(
      dbh           = list(in_range = c(2.5, 80),  out_range = c(0.5, 150)),
      ln_cr_adj     = list(in_range = c(-1.5, 0),  out_range = c(-3, 0)),
      ln_bal_sw_adj = list(in_range = c(-4, 4),    out_range = c(-6, 6)),
      bal_hw        = list(in_range = c(0, 60),    out_range = c(0, 100)),
      ln_csi        = list(in_range = c(-1, 3),    out_range = c(-3, 5)),
      ba_x_rd       = list(in_range = c(0, 60),    out_range = c(0, 100)),
      bal_x_rd      = list(in_range = c(0, 60),    out_range = c(0, 100))
    )
  )
)

if (!MODEL %in% names(model_spec)) {
  cat("This script currently supports:", paste(names(model_spec), collapse = ", "),
      "\nExtend model_spec for", MODEL, ". Exiting.\n")
  quit(save = "no", status = 0)
}
spec <- model_spec[[MODEL]]
FIT_PATH <- file.path(PROJ_ROOT, ifelse(VARIANT == "b1", spec$fit_b1, spec$fit_b2))
stopifnot(file.exists(FIT_PATH))

# --- Load fit, sample posterior draws ---------------------------------------
cat("Loading fit ..."); flush.console()
fit <- readRDS(FIT_PATH)
cat(" done\n")

coef_prefix <- if (MODEL == "hg") "a" else "b"
coef_names <- c(paste0(coef_prefix, 0:8), "sigma")
all_draws <- fit$draws(variables = coef_names, format = "draws_matrix")
sel <- sample.int(nrow(all_draws), min(N_DRAWS, nrow(all_draws)))
draws <- all_draws[sel, , drop = FALSE]
trait_avg <- mean(rowMeans(fit$draws("trait_effect", format = "draws_matrix")[sel, ]))
rm(fit, all_draws); gc()

# --- Predict on a Cartesian extrapolation grid ------------------------------
# Build a sparse grid: 2 levels per covariate (in vs out of range edges) plus
# the in-range midpoint. Adds up to 3^P points. For 7 covariates that is
# 2187 cells which we predict in one matrix operation.

grid_per_cov <- lapply(spec$cov_specs, function(cv) {
  c(cv$out_range[1], mean(cv$in_range), cv$out_range[2])
})

cart <- expand.grid(grid_per_cov, stringsAsFactors = FALSE)
names(cart) <- names(spec$cov_specs)

# Mark each cell as in-range vs out-of-range
out_of_range <- map_dfc(names(spec$cov_specs), function(cn) {
  cv <- spec$cov_specs[[cn]]
  tibble(
    !!paste0("oor_", cn) := cart[[cn]] < cv$in_range[1] | cart[[cn]] > cv$in_range[2]
  )
})
cart <- bind_cols(cart, out_of_range)
cart$n_oor_covariates <- rowSums(cart[, paste0("oor_", names(spec$cov_specs))])

# Predict eta for each grid cell using the same predict_eta as in 60
predict_eta <- function(draws, trait_eff, cart) {
  N_pred  <- nrow(cart)
  N_draws <- nrow(draws)
  as_v <- function(mat, col) as.numeric(mat[, col])
  b0 <- as_v(draws, 1)
  eta_mat <- matrix(b0 + trait_eff, nrow = N_draws, ncol = N_pred)
  if (MODEL == "dg_kue") {
    eta_mat <- eta_mat +
      outer(as_v(draws, 2), log(cart$dbh))            +
      outer(as_v(draws, 3), cart$dbh)                 +
      outer(as_v(draws, 4), cart$ln_cr_adj)           +
      outer(as_v(draws, 5), cart$ln_bal_sw_adj)       +
      outer(as_v(draws, 6), cart$bal_hw)              +
      outer(as_v(draws, 7), cart$ln_csi)              +
      outer(as_v(draws, 8), cart$ba_x_rd)             +
      outer(as_v(draws, 9), cart$bal_x_rd)
  } else if (MODEL == "hg") {
    eta_mat <- eta_mat +
      outer(as_v(draws, 2), log(cart$dbh))  +
      outer(as_v(draws, 3), cart$dbh)       +
      outer(as_v(draws, 4), cart$ln_cr_adj) +
      outer(as_v(draws, 5), cart$bal)       +
      outer(as_v(draws, 6), cart$cspi)
  }
  eta_mat
}

cat("Predicting across", nrow(cart), "grid cells with", nrow(draws), "draws ...\n")
eta_mat <- predict_eta(draws, trait_avg, cart)
sigma_draws <- as.numeric(draws[, "sigma"])
pred_mat <- exp(eta_mat + matrix(sigma_draws^2 / 2,
                                  nrow = nrow(eta_mat), ncol = ncol(eta_mat)))

cart <- cart %>%
  mutate(
    pred_mean = colMeans(pred_mat),
    pred_q5   = apply(pred_mat, 2, quantile, 0.05, na.rm = TRUE),
    pred_q95  = apply(pred_mat, 2, quantile, 0.95, na.rm = TRUE),
    # Biological violation flags
    violates_min      = pred_q5  < spec$response_min,
    violates_max      = pred_q95 > spec$response_max,
    violates_any      = violates_min | violates_max,
    extreme_extrapol  = pred_mean > 2 * spec$response_max | pred_mean < -spec$response_max
  )

# --- Summary stats by extrapolation distance --------------------------------
summary_by_dist <- cart %>%
  group_by(n_oor_covariates) %>%
  summarise(
    n_cells           = n(),
    pct_violates_min  = round(100 * mean(violates_min), 1),
    pct_violates_max  = round(100 * mean(violates_max), 1),
    pct_violates_any  = round(100 * mean(violates_any), 1),
    pct_extreme       = round(100 * mean(extreme_extrapol), 1),
    pred_mean_q5      = round(quantile(pred_mean, 0.05), 4),
    pred_mean_q95     = round(quantile(pred_mean, 0.95), 4),
    pred_mean_max     = round(max(pred_mean), 4),
    pred_mean_min     = round(min(pred_mean), 4),
    .groups = "drop"
  )

write_csv(summary_by_dist, file.path(OUT_DIR, "extrapolation_summary.csv"))
cat("\n=== Extrapolation summary by number of out-of-range covariates ===\n")
print(summary_by_dist)

# --- Cell-level CSV ---------------------------------------------------------
write_csv(cart, file.path(OUT_DIR, "extrapolation_cells.csv"))

# --- Heatmap: predicted mean response by DBH x BAL_HW (or comparable) -------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    strip.text        = element_text(face = "bold"),
    plot.title        = element_text(face = "bold")
  )

heatmap_vars <- if (MODEL == "dg_kue") c("dbh", "bal_hw") else c("dbh", "bal")
h_df <- cart
hm <- h_df %>%
  group_by(.data[[heatmap_vars[1]]], .data[[heatmap_vars[2]]]) %>%
  summarise(pred = mean(pred_mean),
            violates = mean(violates_any),
            .groups = "drop")

p_heat_pred <- ggplot(hm, aes(.data[[heatmap_vars[1]]], .data[[heatmap_vars[2]]],
                              fill = pred)) +
  geom_tile() +
  scale_fill_gradientn(colors = pal_gsea()(12)) +
  labs(x = heatmap_vars[1], y = heatmap_vars[2], fill = "Predicted",
       title = sprintf("A) Predicted response over (%s, %s) plane",
                       heatmap_vars[1], heatmap_vars[2])) +
  theme_pub

p_heat_violate <- ggplot(hm, aes(.data[[heatmap_vars[1]]], .data[[heatmap_vars[2]]],
                                  fill = violates)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#EE0000",
                      limits = c(0, 1), name = "P(violation)") +
  labs(x = heatmap_vars[1], y = heatmap_vars[2],
       title = "B) Rate of biological-violation predictions") +
  theme_pub

heatmap_fig <- p_heat_pred / p_heat_violate +
  plot_annotation(
    title = sprintf("Extrapolation heatmap: %s %s", MODEL, VARIANT),
    subtitle = sprintf("Plausible response range: [%.2f, %.2f] %s",
                       spec$response_min, spec$response_max, spec$response_unit)
  )

ggsave(file.path(OUT_DIR, "extrapolation_heatmap.png"), heatmap_fig,
       width = 18, height = 22, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "extrapolation_heatmap.pdf"), heatmap_fig,
       width = 18, height = 22, units = "cm")

# --- Edge behavior: predictions at the four corners + center ----------------
edge_cells <- cart %>%
  filter(n_oor_covariates == 0 |
         n_oor_covariates == length(spec$cov_specs)) %>%
  arrange(n_oor_covariates) %>%
  slice(c(1, n())) %>%
  mutate(label = c("All covariates in-range",
                   "All covariates out-of-range"))

cat("\n=== Edge behavior comparison ===\n")
print(edge_cells %>% select(label, pred_mean, pred_q5, pred_q95,
                            violates_any, extreme_extrapol))

# --- Final verdict ----------------------------------------------------------
in_range_pct  <- summary_by_dist$pct_violates_any[summary_by_dist$n_oor_covariates == 0]
out_range_pct <- summary_by_dist$pct_violates_any[summary_by_dist$n_oor_covariates ==
                                                    max(summary_by_dist$n_oor_covariates)]

cat("\n=== Verdict ===\n")
cat(sprintf("In-range cells: %.1f%% violate plausibility bounds\n", in_range_pct))
cat(sprintf("Fully out-of-range cells: %.1f%% violate plausibility bounds\n",
            out_range_pct))
if (in_range_pct < 5) {
  cat("In-range behavior is biologically consistent.\n")
} else {
  cat("WARNING: in-range cells produce violations. Investigate the model.\n")
}
if (out_range_pct < 25) {
  cat("Out-of-range behavior degrades gracefully.\n")
} else if (out_range_pct < 60) {
  cat("Out-of-range behavior degrades meaningfully. Predictions at extreme\n")
  cat("covariate values should be flagged in any operational use.\n")
} else {
  cat("WARNING: Severe out-of-range degradation. The model should not be\n")
  cat("extrapolated to these covariate combinations.\n")
}

cat("\nDone. Outputs in", OUT_DIR, "\n")
