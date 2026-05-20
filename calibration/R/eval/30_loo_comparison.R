# =============================================================================
# Title: PSIS-LOO model comparison (B1 species-free vs B2 species-aware)
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: Computes Pareto-smoothed importance sampling leave-one-out
#              cross-validation (PSIS-LOO) for B1 and B2 of a given response,
#              then reports the ELPD difference and approximate LOO model
#              weights. Pareto-k diagnostics flag observations where the LOO
#              estimate may be unreliable.
#
# Stan models save log_lik[i] for each observation in generated quantities,
# which is the input needed for loo::loo().
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_loo.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/30_loo_comparison.R --response=hg
#     # or --response=dg_kue (when DG_Kue B1 is banked)
#
# Memory: 96-128 GB recommended. Both fits must be loaded simultaneously.
# Walltime: 30 min to 2 hours depending on N_obs.
# =============================================================================

# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(loo)
library(posterior)
library(ggsci)

# --- CLI parsing -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
RESPONSE  <- get_arg("response", "hg")
N_CORES   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

cat(sprintf("LOO comparison: response=%s (cores=%d)\n", RESPONSE, N_CORES))

# --- Paths -------------------------------------------------------------------
PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR   <- file.path(PROJ_ROOT, "calibration/output/evaluation/loo")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

paths_hg <- list(
  b2 = file.path(PROJ_ROOT, "calibration/output/conus/hg/hg_organon_fixedK_cspi_traits1_fit.rds"),
  b1 = file.path(PROJ_ROOT, "calibration/output/conus/hg/speciesfree_pilot/hg_organon_fixedK_cspi_traits1_fit.rds")
)
paths_dg_kue <- list(
  b2 = file.path(PROJ_ROOT, "calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds"),
  b1 = file.path(PROJ_ROOT, "calibration/output/conus/dg/speciesfree_pilot/dg_kuehne_cspi_traits1_b1_fit.rds")
)

paths <- switch(RESPONSE,
  "hg"     = paths_hg,
  "dg_kue" = paths_dg_kue,
  stop("Unknown response: ", RESPONSE))

stopifnot(file.exists(paths$b1), file.exists(paths$b2))

# --- Compute LOO for one model ----------------------------------------------

compute_loo <- function(fit_path, label) {
  cat(sprintf("\n=== LOO for %s (%s) ===\n", label, basename(fit_path)))
  cat("Loading fit ..."); flush.console()
  fit <- readRDS(fit_path)
  cat(" done\n")

  cat("Extracting log_lik draws ..."); flush.console()
  log_lik <- fit$draws(variables = "log_lik", format = "draws_matrix")
  cat(" done. Dim:", paste(dim(log_lik), collapse = " x "), "\n")

  # Drop fit object before LOO (free memory before chunked computation)
  rm(fit); gc()

  cat("Computing relative ESS ..."); flush.console()
  r_eff <- relative_eff(exp(log_lik), chain_id = rep(seq_len(4),
                       each = nrow(log_lik) / 4))
  cat(" done\n")

  cat("Computing PSIS-LOO (this can take a while) ..."); flush.console()
  loo_result <- loo(log_lik, r_eff = r_eff, cores = N_CORES)
  cat(" done\n")

  rm(log_lik); gc()
  loo_result
}

loo_b1 <- compute_loo(paths$b1, "B1 species-free")
loo_b2 <- compute_loo(paths$b2, "B2 species-aware")

# --- Comparison and Pareto-k diagnostics ------------------------------------

cmp <- loo_compare(list(B1 = loo_b1, B2 = loo_b2))
cat("\n=== LOO comparison (B1 vs B2) ===\n")
print(cmp)

# Approximate model weights via stacking
weights <- loo_model_weights(list(B1 = loo_b1, B2 = loo_b2))
cat("\n=== LOO model weights ===\n")
print(weights)

# Pareto-k diagnostics
pk_b1 <- pareto_k_table(loo_b1)
pk_b2 <- pareto_k_table(loo_b2)
cat("\n=== Pareto-k summary ===\n")
cat("B1:\n"); print(pk_b1)
cat("\nB2:\n"); print(pk_b2)

# --- Save results -----------------------------------------------------------
saveRDS(list(
  loo_b1 = loo_b1, loo_b2 = loo_b2,
  comparison = cmp, weights = weights
), file.path(OUT_DIR, paste0("loo_", RESPONSE, "_b1_b2.rds")))

cmp_df <- as_tibble(cmp, rownames = "model") %>%
  mutate(across(where(is.numeric), ~round(., 3)))
write_csv(cmp_df, file.path(OUT_DIR, paste0("loo_", RESPONSE, "_b1_b2_table.csv")))

# --- Figure: Pareto-k distribution ------------------------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold")
  )

pk_df <- bind_rows(
  tibble(model = "B1 species-free", pareto_k = loo_b1$diagnostics$pareto_k),
  tibble(model = "B2 species-aware", pareto_k = loo_b2$diagnostics$pareto_k)
)
p_pk <- ggplot(pk_df, aes(pareto_k, fill = model)) +
  geom_histogram(bins = 60, alpha = 0.65, position = "identity") +
  geom_vline(xintercept = c(0.5, 0.7, 1.0), linetype = "dashed",
             color = c("grey30", "grey50", "red"), linewidth = 0.4) +
  scale_fill_aaas() +
  labs(x = "Pareto k", y = "Count",
       fill = "Variant",
       title = sprintf("Pareto-k diagnostic for %s LOO", toupper(RESPONSE)),
       subtitle = "k < 0.5 (good), 0.5 to 0.7 (ok), > 0.7 (problematic)") +
  theme_pub

fig_png <- file.path(OUT_DIR, paste0("pareto_k_", RESPONSE, ".png"))
fig_pdf <- file.path(OUT_DIR, paste0("pareto_k_", RESPONSE, ".pdf"))
ggsave(fig_png, p_pk, width = 18, height = 10, units = "cm",
       dpi = 300, bg = "white")
ggsave(fig_pdf, p_pk, width = 18, height = 10, units = "cm")
cat("Wrote", fig_png, "\n")
cat("Wrote", fig_pdf, "\n")

# --- Headline narrative -----------------------------------------------------
elpd_diff <- cmp_df$elpd_diff
elpd_se   <- cmp_df$se_diff
winner    <- cmp_df$model[1]
elpd_ratio <- abs(elpd_diff[2]) / elpd_se[2]

cat("\n=== Verdict ===\n")
cat(sprintf("Winner by ELPD: %s\n", winner))
cat(sprintf("ELPD difference (loser - winner): %.2f +/- %.2f (%.1f SE)\n",
            elpd_diff[2], elpd_se[2], elpd_ratio))
if (elpd_ratio < 2) {
  cat("The two models are indistinguishable in LOO predictive accuracy.\n")
  cat("This is STRONG evidence in favor of the species-free architecture if\n")
  cat("B1 is the winner or tied: traits replace species RE without cost.\n")
} else if (elpd_ratio > 4) {
  cat("There is a clear LOO preference for the winner.\n")
} else {
  cat("There is a mild LOO preference for the winner.\n")
}

cat("\nDone.\n")
