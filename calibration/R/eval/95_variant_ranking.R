# =============================================================================
# Title: Architecture variant ranking for DG_Kuehne B1 species-free
# Author: A. Weiskittel
# Date: 2026-05-13
# Description: Reads the four architecture variants (v0 baseline + v1 quad +
#              v2 L1-varying site + v3 trait-site + v4 full) and ranks them
#              on:
#                 - sigma posterior (lower = tighter fit)
#                 - LOO ELPD difference (computed if log_lik available)
#                 - biological realism: site coefficient sign + magnitude
#                 - convergence (rhat, ESS)
#
# Run on Cardinal after all variant fits land:
#   sbatch calibration/slurm/run_variant_ranking.sh
#     # or directly:
#     module load gcc/12.3.0 R/4.4.0
#     cd ~/fvs-modern
#     Rscript --vanilla calibration/R/eval/95_variant_ranking.R
# =============================================================================

library(tidyverse)
library(posterior)
library(ggsci)
library(patchwork)

PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
VAR_DIR   <- file.path(PROJ_ROOT, "calibration/output/conus/dg_kue/architecture_variants")
B1_BASE   <- file.path(PROJ_ROOT, "calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds")
OUT_DIR   <- file.path(PROJ_ROOT, "calibration/output/evaluation/variant_ranking")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

variants <- tribble(
  ~variant_id, ~variant_label,           ~fit_path,
  "v0",        "baseline (linear)",      B1_BASE,
  "v1_quad",   "linear + quadratic",     file.path(VAR_DIR, "dg_kuehne_b1_v1_quad_fit.rds"),
  "v2_l1site", "L1-varying site slope",  file.path(VAR_DIR, "dg_kuehne_b1_v2_l1site_fit.rds"),
  "v3_traitsite", "trait x site",        file.path(VAR_DIR, "dg_kuehne_b1_v3_traitsite_fit.rds"),
  "v4_full",   "full combined",          file.path(VAR_DIR, "dg_kuehne_b1_v4_full_fit.rds")
)

# --- Read summaries (sigma, b6, b9, b10, rhat) from each variant ----------

read_variant <- function(variant_id, variant_label, fit_path) {
  if (!file.exists(fit_path)) {
    message(sprintf("  SKIP %s | file not present", variant_id))
    return(NULL)
  }
  summary_csv <- sub("_fit\\.rds$", "_summary.csv", fit_path)
  if (file.exists(summary_csv)) {
    df <- read_csv(summary_csv, show_col_types = FALSE)
  } else {
    message(sprintf("  No summary CSV for %s; reading fit ...", variant_id))
    fit <- readRDS(fit_path)
    vars <- c(paste0("b", 0:10), paste0("gamma[", 1:8, "]"),
              paste0("gamma_site[", 1:8, "]"),
              "sigma", "sigma_L1", "sigma_L2", "sigma_L3", "sigma_L1_csi")
    vars <- intersect(vars, variables(fit$draws()))
    df <- fit$summary(variables = vars, "mean", "median", "sd",
                       ~quantile(.x, c(0.05, 0.95), na.rm = TRUE),
                       "rhat", "ess_bulk", "ess_tail")
    names(df)[names(df) %in% c("5%", "95%")] <- c("q5", "q95")
    rm(fit); gc()
  }
  df %>% mutate(variant_id = variant_id, variant_label = variant_label)
}

all_summaries <- pmap_dfr(variants, read_variant)
if (nrow(all_summaries) == 0) {
  stop("No variant fits found. Wait for the fits to land, then rerun.")
}

# --- Sigma comparison -------------------------------------------------------
sigma_cmp <- all_summaries %>%
  filter(variable == "sigma") %>%
  select(variant_id, variant_label, mean, sd, q5, q95, rhat)
write_csv(sigma_cmp, file.path(OUT_DIR, "variant_sigma_comparison.csv"))
cat("\n=== Sigma posterior by variant ===\n")
print(sigma_cmp)

# --- Site coefficients (b6, b9, b10) ---------------------------------------
site_coefs <- all_summaries %>%
  filter(variable %in% c("b6", "b9", "b10")) %>%
  mutate(coef_role = case_when(
    variable == "b6"  ~ "linear site (b6)",
    variable == "b9"  ~ "quadratic site (b9)",
    variable == "b10" ~ "site x BAL (b10)"
  )) %>%
  select(variant_id, variant_label, coef_role, mean, sd, q5, q95, rhat)
write_csv(site_coefs, file.path(OUT_DIR, "variant_site_coefs.csv"))
cat("\n=== Site coefficient posteriors ===\n")
print(site_coefs)

# --- Biological realism flag: linear site coefficient should be positive ---
realism <- sigma_cmp %>%
  select(variant_id, variant_label, sigma_mean = mean) %>%
  left_join(
    site_coefs %>% filter(coef_role == "linear site (b6)") %>%
      select(variant_id, b6_mean = mean, b6_q5 = q5, b6_q95 = q95),
    by = "variant_id"
  ) %>%
  mutate(
    site_sign_correct = b6_q5 > 0,    # 90% CI strictly positive
    site_status = case_when(
      b6_q5 > 0      ~ "positive (biologically plausible)",
      b6_q95 < 0     ~ "negative (implausible for increment)",
      TRUE           ~ "spans zero"
    )
  )
write_csv(realism, file.path(OUT_DIR, "variant_biological_realism.csv"))
cat("\n=== Biological realism: site coefficient sign ===\n")
print(realism)

# --- LOO comparison (heavy: loads all 5 fits) ------------------------------
if (require("loo", quietly = TRUE)) {
  loo_list <- list()
  for (i in seq_len(nrow(variants))) {
    fp <- variants$fit_path[i]
    if (!file.exists(fp)) next
    message(sprintf("Computing LOO for %s ...", variants$variant_id[i]))
    fit <- readRDS(fp)
    log_lik <- tryCatch(
      fit$draws(variables = "log_lik", format = "draws_matrix"),
      error = function(e) NULL
    )
    if (is.null(log_lik)) { rm(fit); gc(); next }
    r_eff <- loo::relative_eff(exp(log_lik),
              chain_id = rep(seq_len(4), each = nrow(log_lik) / 4))
    loo_list[[variants$variant_id[i]]] <- loo::loo(log_lik, r_eff = r_eff, cores = 4)
    rm(fit, log_lik); gc()
  }
  if (length(loo_list) >= 2) {
    cmp <- loo::loo_compare(loo_list)
    cat("\n=== LOO ELPD comparison ===\n"); print(cmp)
    saveRDS(loo_list, file.path(OUT_DIR, "variant_loo_objects.rds"))
    cmp_df <- as_tibble(cmp, rownames = "variant_id")
    write_csv(cmp_df, file.path(OUT_DIR, "variant_loo_comparison.csv"))
  }
}

# --- Composite ranking ------------------------------------------------------
ranking <- realism %>%
  mutate(
    score_sigma = -as.numeric(scale(sigma_mean)),  # lower sigma = better
    score_site_sign = as.integer(site_sign_correct),
    composite_score = score_sigma + 2 * score_site_sign
  ) %>%
  arrange(desc(composite_score))
write_csv(ranking, file.path(OUT_DIR, "variant_ranking.csv"))
cat("\n=== Composite ranking ===\n")
print(ranking)

# --- Figure: sigma + b6 sign per variant -----------------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold"),
    strip.text         = element_text(face = "bold")
  )

p_sigma <- ggplot(sigma_cmp, aes(reorder(variant_label, -mean), mean)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.25, color = "#3B4992") +
  geom_point(size = 3, color = "#3B4992") +
  coord_flip() +
  labs(x = NULL, y = "Posterior sigma (90% CI)",
       title = "A) Residual sigma by architecture variant") +
  theme_pub

p_site <- site_coefs %>%
  filter(coef_role == "linear site (b6)") %>%
  ggplot(aes(reorder(variant_label, mean), mean)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.25, color = "#EE0000") +
  geom_point(size = 3, color = "#EE0000") +
  coord_flip() +
  labs(x = NULL, y = "b6 linear site coef (90% CI)",
       title = "B) Linear site coefficient (positive = biologically plausible)") +
  theme_pub

combined <- p_sigma / p_site +
  plot_annotation(
    title = "DG_Kuehne species-free architecture variant comparison",
    subtitle = "Lower sigma + positive site coefficient is the goal"
  )

ggsave(file.path(OUT_DIR, "variant_comparison.png"), combined,
       width = 22, height = 18, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "variant_comparison.pdf"), combined,
       width = 22, height = 18, units = "cm")

cat("\n=== Done. Outputs in", OUT_DIR, "===\n")
