#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: Extract Posteriors from Existing RDS Files
#
# Recovery script for when model fitting succeeded (RDS saved) but posterior
# extraction crashed (e.g., due to posterior::summarise_draws column naming).
#
# Loads mortality_samples.rds and crown_ratio_samples.rds, extracts point
# estimates and summaries, saves CSV outputs.
#
# Usage: Rscript calibration/R/extract_posteriors.R --variant ne
# =============================================================================

library(tidyverse)
library(posterior)
library(logger)

args <- commandArgs(trailingOnly = TRUE)
variant <- "ne"

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variant <- args[i + 1]
    }
  }
}

project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (!nzchar(project_root)) {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
output_dir <- file.path(project_root, "calibration", "output", "variants", variant)

logger::log_info("Extracting posteriors for variant: {variant}")
logger::log_info("Output dir: {output_dir}")

# =============================================================================
# Helper: robust summarise_draws with correct column naming
# =============================================================================
safe_summarise <- function(draws, model_name) {
  summ <- tryCatch({
    posterior::summarise_draws(draws,
      p50 = function(x) unname(median(x)),
      p05 = function(x) unname(quantile(x, 0.05)),
      p95 = function(x) unname(quantile(x, 0.95)),
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk
    )
  }, error = function(e) {
    logger::log_warn("{model_name}: summarise_draws failed ({e$message}), using manual")
    draws_mat <- as.matrix(draws)
    tibble(
      variable = colnames(draws_mat),
      p50 = apply(draws_mat, 2, function(x) unname(median(x))),
      p05 = apply(draws_mat, 2, function(x) unname(quantile(x, 0.05))),
      p95 = apply(draws_mat, 2, function(x) unname(quantile(x, 0.95))),
      rhat = NA_real_,
      ess_bulk = NA_real_
    )
  })

  # Belt and suspenders: fix column names
  names(summ) <- gsub("^5%$", "p05", names(summ))
  names(summ) <- gsub("^95%$", "p95", names(summ))
  if (!"variable" %in% names(summ)) names(summ)[1] <- "variable"

  summ %>%
    mutate(
      ci_width = p95 - p05,
      converged = if_else(!is.na(rhat), rhat < 1.05, NA)
    )
}

# =============================================================================
# Extract Mortality Posteriors
# =============================================================================
mort_rds <- file.path(output_dir, "mortality_samples.rds")
if (file.exists(mort_rds)) {
  logger::log_info("Loading mortality posterior draws...")
  draws_mort <- readRDS(mort_rds)

  summ_mort <- safe_summarise(draws_mort, "mortality")

  write_csv(summ_mort, file.path(output_dir, "mortality_posterior.csv"))
  logger::log_info("Saved mortality_posterior.csv ({nrow(summ_mort)} parameters)")

  map_mort <- summ_mort %>% select(variable, estimate = p50)
  write_csv(map_mort, file.path(output_dir, "mortality_map.csv"))
  logger::log_info("Saved mortality_map.csv")
} else {
  logger::log_warn("No mortality_samples.rds found for {variant}")
}

# =============================================================================
# Extract Crown Ratio Posteriors
# =============================================================================
cr_rds <- file.path(output_dir, "crown_ratio_samples.rds")
if (file.exists(cr_rds)) {
  logger::log_info("Loading crown ratio posterior draws...")
  draws_cr <- readRDS(cr_rds)

  summ_cr <- safe_summarise(draws_cr, "crown_ratio")

  write_csv(summ_cr, file.path(output_dir, "crown_ratio_posterior.csv"))
  logger::log_info("Saved crown_ratio_posterior.csv ({nrow(summ_cr)} parameters)")

  map_cr <- summ_cr %>% select(variable, estimate = p50)
  write_csv(map_cr, file.path(output_dir, "crown_ratio_map.csv"))
  logger::log_info("Saved crown_ratio_map.csv")
} else {
  logger::log_warn("No crown_ratio_samples.rds found for {variant}")
}

# =============================================================================
# Extract Height-Diameter Posteriors (if brms RDS exists)
# =============================================================================
hd_rds <- file.path(output_dir, "height_diameter_fit.rds")
if (file.exists(hd_rds)) {
  logger::log_info("Loading height-diameter model fit...")
  fit_hd <- readRDS(hd_rds)

  draws_hd <- tryCatch(as_draws_df(fit_hd), error = function(e) NULL)
  if (!is.null(draws_hd)) {
    summ_hd <- safe_summarise(draws_hd, "height_diameter")
    write_csv(summ_hd, file.path(output_dir, "height_diameter_posterior.csv"))
    map_hd <- summ_hd %>% select(variable, estimate = p50)
    write_csv(map_hd, file.path(output_dir, "height_diameter_map.csv"))
    logger::log_info("Saved height_diameter posterior and MAP")
  }
} else {
  logger::log_info("No height_diameter_fit.rds (Step 03 may have failed)")
}

logger::log_info("Extraction complete for {variant}")
