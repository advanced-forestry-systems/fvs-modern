# =============================================================================
# Title: Marginal response curves with biological-expectation overlays
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: For a fitted FVS-CONUS base model, generate one curve per
#              covariate showing the predicted response as that covariate
#              varies (all others held at typical values). Overlay reference
#              expectations from forest biometrics literature (e.g., growth
#              should decrease with BAL, mortality should be U-shaped with
#              DBH). Posterior credible intervals come from the full draws.
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_response.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/60_response_curves.R --model=hg --variant=b2
#
# CLI flags:
#   --model={hg,dg_kue,htdbh,cr,hcb,mort}
#   --variant={b1,b2}
#   --species=SPCD     focus species (default: 10 most common in training)
#   --n_curves=N       number of posterior draws to use for CI (default: 200)
#
# Outputs:
#   calibration/output/evaluation/response/{model}_{variant}/
#     response_curves_<covariate>.png
#     response_curves_<covariate>.pdf
#     summary.csv
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
N_CURVES <- as.integer(get_arg("n_curves", "200"))

cat(sprintf("Response curves: model=%s variant=%s n_curves=%d\n",
            MODEL, VARIANT, N_CURVES))

# --- Paths -------------------------------------------------------------------
PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR <- file.path(PROJ_ROOT, "calibration/output/evaluation/response",
                     paste(MODEL, VARIANT, sep = "_"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Model file lookups + covariate definitions per model
model_spec <- list(
  hg = list(
    fit_b2 = "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds",
    fit_b1 = "calibration/output/conus/hg/speciesfree_pilot/hg_organon_fixedK_cspi_traits1_fit.rds",
    response = "ln_hg_annual",
    link = "lognormal",
    cov_specs = list(
      dbh        = list(label = "DBH (cm)",         range = c(2.5, 80),  typical = 25),
      ln_cr_adj  = list(label = "CR (adjusted, log)", range = c(-1.5, 0), typical = -0.4),
      bal        = list(label = "BAL (m2/ha)",      range = c(0, 60),    typical = 15),
      cspi       = list(label = "CSPI (z-score)",   range = c(-2, 4),    typical = 0)
    ),
    biological_expectation = list(
      dbh        = "Hump-shaped (juvenile peak then decline)",
      ln_cr_adj  = "Monotonic increase",
      bal        = "Monotonic decrease",
      cspi       = "Monotonic increase"
    )
  ),
  dg_kue = list(
    fit_b2 = "calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds",
    fit_b1 = "calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds",
    response = "ln_dg_annual",
    link = "lognormal",
    cov_specs = list(
      dbh           = list(label = "DBH (cm)",            range = c(2.5, 80),  typical = 25),
      ln_cr_adj     = list(label = "CR (adjusted, log)",  range = c(-1.5, 0),  typical = -0.4),
      ln_bal_sw_adj = list(label = "BAL softwood (log)",  range = c(-4, 4),    typical = 0),
      bal_hw        = list(label = "BAL hardwood (m2/ha)", range = c(0, 60),   typical = 15),
      ln_csi        = list(label = "CSPI (log)",          range = c(-1, 3),    typical = 1),
      ba_x_rd       = list(label = "BA x RD interaction", range = c(0, 60),    typical = 15),
      bal_x_rd      = list(label = "BAL x RD interaction", range = c(0, 60),   typical = 15)
    ),
    biological_expectation = list(
      dbh           = "Hump-shaped (juvenile peak then decline at large DBH)",
      ln_cr_adj     = "Monotonic increase",
      ln_bal_sw_adj = "Monotonic decrease",
      bal_hw        = "Decrease, weaker than BAL_SW",
      ln_csi        = "Monotonic increase",
      ba_x_rd       = "Monotonic decrease (density effect)",
      bal_x_rd      = "Monotonic decrease"
    )
  )
)

if (!MODEL %in% names(model_spec)) {
  cat("Currently this script supports: ", paste(names(model_spec), collapse = ", "),
      "\nFor", MODEL, ", extend the model_spec list. Exiting.\n")
  quit(save = "no", status = 0)
}

spec <- model_spec[[MODEL]]
FIT_PATH <- file.path(PROJ_ROOT,
                      ifelse(VARIANT == "b1", spec$fit_b1, spec$fit_b2))
stopifnot(file.exists(FIT_PATH))

# --- Load fit + draws of all fixed-effect coefficients ----------------------
cat("Loading fit:", FIT_PATH, "\n"); flush.console()
fit <- readRDS(FIT_PATH)

# For Kuehne form, parameters are b0..b8 + gamma + sigma + sigma_L1/L2/L3.
# For HG/ORGANON form, a0..a8 + gamma + sigma + sigma_L1/L2/L3.
coef_prefix <- if (MODEL == "hg") "a" else "b"
coef_names <- c(paste0(coef_prefix, 0:8), "sigma")
draws <- fit$draws(variables = coef_names, format = "draws_matrix")
n_total_draws <- nrow(draws)
sel <- sample.int(n_total_draws, min(N_CURVES, n_total_draws))
draws_sel <- draws[sel, , drop = FALSE]
cat("Pulled", nrow(draws_sel), "posterior draws of",
    length(coef_names), "coefficients\n")

# trait_effect is a vector indexed by species
trait_effect_draws <- tryCatch({
  fit$draws(variables = "trait_effect", format = "draws_matrix")[sel, , drop = FALSE]
}, error = function(e) NULL)

rm(fit, draws); gc()

# --- Build prediction function ----------------------------------------------
# For Kuehne form:
#   eta = b0 + trait_effect[sp] + b1 ln(dbh) + b2 dbh + b3 ln_cr_adj
#       + b4 ln_bal_sw_adj + b5 bal_hw + b6 ln_csi + b7 ba_x_rd + b8 bal_x_rd
# For HG/ORGANON form:
#   eta = a0 + trait_effect[sp] + a1 ln(dbh) + a2 dbh + a3 ln_cr_adj + a4 bal
#       + a5 cspi + a6 ... (model-specific)

predict_eta <- function(draws_sel, trait_eff, covariate_grid) {
  # covariate_grid: data frame with one row per evaluation point, columns named
  # to match cov_specs keys.
  N_pred  <- nrow(covariate_grid)
  N_draws <- nrow(draws_sel)
  eta_mat <- matrix(0, nrow = N_draws, ncol = N_pred)
  # Coerce each posterior column to a plain numeric vector. draws_matrix
  # subsetting returns a draws_matrix not a vector, which breaks outer().
  as_v <- function(mat, col) as.numeric(mat[, col])
  b0_draws <- as_v(draws_sel, 1)
  eta_mat <- eta_mat + matrix(b0_draws + trait_eff, nrow = N_draws, ncol = N_pred)
  if (MODEL == "dg_kue") {
    eta_mat <- eta_mat +
      outer(as_v(draws_sel, 2), log(covariate_grid$dbh))      +
      outer(as_v(draws_sel, 3), covariate_grid$dbh)           +
      outer(as_v(draws_sel, 4), covariate_grid$ln_cr_adj)     +
      outer(as_v(draws_sel, 5), covariate_grid$ln_bal_sw_adj) +
      outer(as_v(draws_sel, 6), covariate_grid$bal_hw)        +
      outer(as_v(draws_sel, 7), covariate_grid$ln_csi)        +
      outer(as_v(draws_sel, 8), covariate_grid$ba_x_rd)       +
      outer(as_v(draws_sel, 9), covariate_grid$bal_x_rd)
  } else if (MODEL == "hg") {
    eta_mat <- eta_mat +
      outer(as_v(draws_sel, 2), log(covariate_grid$dbh))  +
      outer(as_v(draws_sel, 3), covariate_grid$dbh)       +
      outer(as_v(draws_sel, 4), covariate_grid$ln_cr_adj) +
      outer(as_v(draws_sel, 5), covariate_grid$bal)       +
      outer(as_v(draws_sel, 6), covariate_grid$cspi)
  }
  eta_mat
}

# --- Marginal response curves: vary one covariate, hold others typical ------
# Use mean trait_effect (averaged across species) as a "typical species" anchor.

typical_trait <- if (!is.null(trait_effect_draws)) {
  rowMeans(trait_effect_draws)
} else {
  rep(0, nrow(draws_sel))
}

theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold"),
    plot.title.position = "plot"
  )

summary_rows <- list()

for (cov_name in names(spec$cov_specs)) {
  cv <- spec$cov_specs[[cov_name]]
  cat(sprintf("Building response curve for %s ...\n", cov_name))

  grid <- as_tibble(map(spec$cov_specs, ~.x$typical))
  names(grid) <- names(spec$cov_specs)
  grid_full <- grid[rep(1, 100), ]
  grid_full[[cov_name]] <- seq(cv$range[1], cv$range[2], length.out = 100)

  eta_mat <- predict_eta(draws_sel, typical_trait, grid_full)

  # Lognormal back-transform: mean = exp(mu + sigma^2/2). For visualization,
  # we plot on the natural response scale.
  sigma_draws <- as.numeric(draws_sel[, "sigma"])
  mu_pred <- eta_mat
  pred_mat <- exp(mu_pred + matrix(sigma_draws^2 / 2,
                                    nrow = nrow(mu_pred), ncol = ncol(mu_pred)))

  summary_df <- tibble(
    cov_value = grid_full[[cov_name]],
    pred_mean = colMeans(pred_mat),
    pred_q5   = apply(pred_mat, 2, quantile, 0.05),
    pred_q50  = apply(pred_mat, 2, quantile, 0.50),
    pred_q95  = apply(pred_mat, 2, quantile, 0.95)
  )

  p <- ggplot(summary_df, aes(cov_value, pred_q50)) +
    geom_ribbon(aes(ymin = pred_q5, ymax = pred_q95),
                fill = "#3B4992", alpha = 0.2) +
    geom_line(color = "#3B4992", linewidth = 0.7) +
    labs(x = cv$label,
         y = "Predicted annual response (natural scale)",
         title = sprintf("Marginal response: %s vs %s",
                         spec$response, cov_name),
         subtitle = sprintf("Biological expectation: %s",
                            spec$biological_expectation[[cov_name]])) +
    theme_pub

  fig_png <- file.path(OUT_DIR, paste0("response_", cov_name, ".png"))
  fig_pdf <- file.path(OUT_DIR, paste0("response_", cov_name, ".pdf"))
  ggsave(fig_png, p, width = 14, height = 10, units = "cm",
         dpi = 300, bg = "white")
  ggsave(fig_pdf, p, width = 14, height = 10, units = "cm")

  summary_rows[[cov_name]] <- summary_df %>%
    mutate(covariate = cov_name,
           biological_expectation = spec$biological_expectation[[cov_name]])
}

# --- Save combined summary --------------------------------------------------
combined <- bind_rows(summary_rows) %>%
  mutate(model = MODEL, variant = VARIANT) %>%
  select(model, variant, covariate, biological_expectation,
         cov_value, pred_mean, pred_q5, pred_q50, pred_q95)
write_csv(combined, file.path(OUT_DIR, "response_curves_summary.csv"))

# --- Multi-panel summary figure ---------------------------------------------
p_panels <- lapply(names(spec$cov_specs), function(cn) {
  d <- combined %>% filter(covariate == cn)
  cv <- spec$cov_specs[[cn]]
  ggplot(d, aes(cov_value, pred_q50)) +
    geom_ribbon(aes(ymin = pred_q5, ymax = pred_q95),
                fill = "#3B4992", alpha = 0.2) +
    geom_line(color = "#3B4992", linewidth = 0.7) +
    labs(x = cv$label, y = "Predicted",
         subtitle = spec$biological_expectation[[cn]]) +
    theme_pub +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))
})
panels_fig <- wrap_plots(p_panels, ncol = 2) +
  plot_annotation(
    title = sprintf("Response curves: %s %s", MODEL, VARIANT),
    subtitle = "Posterior median (line), 90% credible interval (ribbon). All other covariates held at typical values."
  )

n_panels <- length(p_panels)
fig_height <- 8 + 6 * ceiling(n_panels / 2)
ggsave(file.path(OUT_DIR, "response_panels.png"), panels_fig,
       width = 20, height = fig_height, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "response_panels.pdf"), panels_fig,
       width = 20, height = fig_height, units = "cm")
cat("Wrote response curves panel figure.\n")

# --- Biological-realism scorecard -------------------------------------------
# For each covariate, classify the actual model behavior:
#   monotone_inc, monotone_dec, hump_shaped, u_shaped, flat
# and compare to the expectation in spec$biological_expectation.

classify_shape <- function(y) {
  if (any(is.na(y)) || length(y) < 5) return("undefined")
  d <- diff(y)
  if (all(d > 0))       "monotone_inc"
  else if (all(d < 0))  "monotone_dec"
  else if (mean(d > 0) > 0.9) "weakly_inc"
  else if (mean(d < 0) > 0.9) "weakly_dec"
  else if (which.max(y) > 5 && which.max(y) < length(y) - 5) "hump_shaped"
  else if (which.min(y) > 5 && which.min(y) < length(y) - 5) "u_shaped"
  else "non_monotone"
}
shape_match <- function(actual, expected) {
  case_when(
    grepl("[Mm]onotonic increase", expected) & actual %in% c("monotone_inc", "weakly_inc") ~ "match",
    grepl("[Mm]onotonic decrease", expected) & actual %in% c("monotone_dec", "weakly_dec") ~ "match",
    grepl("[Hh]ump", expected) & actual == "hump_shaped" ~ "match",
    grepl("[Uu]-shaped", expected) & actual == "u_shaped" ~ "match",
    TRUE ~ "review"
  )
}
realism <- combined %>%
  group_by(model, variant, covariate, biological_expectation) %>%
  summarise(
    actual_shape = classify_shape(pred_q50),
    .groups = "drop"
  ) %>%
  mutate(match_status = shape_match(actual_shape, biological_expectation))
write_csv(realism, file.path(OUT_DIR, "biological_realism_scorecard.csv"))

cat("\n=== Biological realism scorecard ===\n")
print(realism)

cat("\nDone. Outputs in", OUT_DIR, "\n")
