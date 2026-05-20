#!/usr/bin/env Rscript
# =============================================================================
# 61_extract_conus_summaries.R
#
# Extract posterior summaries from a CONUS Phase 4 component fit into the
# CSVs used by 62_conus_to_variant_json.R for per-variant integration.
#
# Reads:
#   calibration/output/conus/{component}/{model}_cspi_traits1_fit.rds
#   calibration/output/conus/{component}/modifier/{model}_modifier_lambda{lam}_fit.rds
#
# Writes (under calibration/output/conus/):
#   {model}_fixed_draws.csv         full posterior of fixed effects
#   {model}_fixed_summary.csv       mean / sd / q025 / q500 / q975 per fixed param
#   {model}_species_intercepts.csv  SPCD-keyed posterior mean random intercept
#   {model}_ecodiv_intercepts.csv   ecodivision-keyed posterior mean random intercept
#   {model}_modifier_summary.csv    per-lambda modifier coefficient summary
#
# Usage:
#   Rscript 61_extract_conus_summaries.R --component dg --model dg_kuehne
#   Rscript 61_extract_conus_summaries.R --component hg --model hg_organon_fixedK
#   Rscript 61_extract_conus_summaries.R --component all   # iterate all 6
#
# Tested model patterns per component:
#   dg              dg_kuehne, dg_organon
#   hg              hg_organon_fixedK
#   ht_dbh          htdbh_wykoff_lognormal
#   hcb             hcb_organon
#   mortality       mort_logit_simple, mort_gompit
#   crown_recession cr_recession
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(posterior)
  library(jsonlite)
  library(logger)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

project_root    <- Sys.getenv("FVS_PROJECT_ROOT", normalizePath(".."))
calibration_dir <- file.path(project_root, "calibration")
conus_dir       <- file.path(calibration_dir, "output", "conus")

# Component -> default model mapping (override via --model)
DEFAULT_MODELS <- list(
  dg              = "dg_kuehne_cspi_traits1",
  hg              = "hg_organon_fixedK_cspi_traits1",
  ht_dbh          = "htdbh_wykoff_lognormal_cspi_traits1",
  hcb             = "hcb_organon_cspi_traits1",
  mortality       = "mort_logit_simple_cspi_traits1",
  crown_recession = "cr_recession_cspi_traits1"
)

LAMBDAS_TO_SUMMARIZE <- c(5, 10, 20)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  out <- list(component = NULL, model = NULL, all = FALSE)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--component" && i < length(args)) {
      out$component <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--model" && i < length(args)) {
      out$model <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--all") {
      out$all <- TRUE; i <- i + 1
    } else {
      i <- i + 1
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Core extractors
# ---------------------------------------------------------------------------

extract_fixed_effects <- function(fit, model_id, out_dir) {
  ## Identify fixed effect parameters by name pattern. Adjust per model family.
  draws_df <- as_draws_df(fit)
  param_names <- variables(draws_df)

  ## Rule 1: name must not contain a bracket — that excludes all indexed
  ## random effects (z_sp[i], z_eco[i], b_sp[i,j], etc.) which would
  ## otherwise blow up the draws CSV to multi-GB.
  no_index <- !grepl("\\[", param_names)

  ## Rule 2: name must match a top-level fixed-effect pattern. The CONUS
  ## Phase 4 models use varying conventions, so the pattern is broad but
  ## still excludes a few well-known non-parameter scalars (lp__, log_lik,
  ## any "raw" parameterizations from non-centered draws).
  ##
  ## Names that should match (across the 6 components, verified from real fits
  ## via SLURM inspection job 9347190 on 2026-05-11):
  ## - DG Kuehne / Organon: b1..b10, K1, K2, mu_b0 etc., sigma*, gamma
  ## - HG Organon fixedK:   a0..a8, mu_a, sigma_sp, sigma_L1..L3, sigma
  ## - HT-DBH Wykoff:       b0, b1, a_bal, a_ba, a_cspi, a_bard, a_blrd, gamma,
  ##                        trait_effect, sigma_sp, sigma_L1, sigma
  ## - HCB Organon:         h0..h5, gamma, sigma_sp, sigma_L1..L3, phi
  ## - Mortality:           m0..m6, gamma, sigma_sp, sigma_L1
  ## - Crown recession:     r0..r6, gamma, sigma_sp, sigma_L1..L3, phi
  ##
  ## The per-component prefixes h/m/r reflect the unified parameterization
  ## introduced in the cspi_traits1 fitting round.
  pattern_ok <- grepl(
    paste0(
      "^(",
      "mu_[a-z_][a-z0-9_]*$|",          # mu_b0, mu_a, mu_alpha, mu_sp
      "[abhmr][0-9]+$|",                # a0..a8, b0..b10, h0..h5, m0..m6, r0..r6
      "K[0-9]+[a-z]?$|",                # K1, K2, K1h, K4h
      "sigma$|sigma_[a-zA-Z0-9_]+$|",   # sigma, sigma_sp, sigma_eco, sigma_L1
      "alpha$|alpha_[a-z]+$|",
      "beta_[a-zA-Z0-9_]+$|",           # beta_dbh, beta_ba, beta_csi
      "gamma$|gamma_[a-zA-Z0-9_]+$|",   # bare gamma (cspi_traits1 fits) or suffixed
      "delta_[a-zA-Z0-9_]+$|",
      "tau$|tau_[a-z]+$|",
      "phi$|nu$|kappa$|",
      "a_[a-z][a-z0-9_]*$|",            # a_bal, a_ba, a_cspi, a_bard, a_blrd (HT-DBH)
      "trait_effect$",                  # cspi_traits1 main effect coefficient
      ")"
    ),
    param_names
  )

  ## Explicit blocklist of non-parameter scalars and raw/auxiliary terms.
  ## Note: trait_effect is a real fixed effect in cspi_traits1 fits and is
  ## explicitly allowed by pattern_ok above; only the "_raw" non-centered
  ## reparam variables and likelihood arrays are blocked here.
  is_aux <- param_names %in% c("lp__", "lp_approx__") |
            grepl("_raw$|^log_lik|^mu_pred|^p_die", param_names)

  is_fixed <- no_index & pattern_ok & !is_aux
  fixed_names <- param_names[is_fixed]
  log_info("Fixed-effect filter: kept {length(fixed_names)} of {length(param_names)} parameters")

  if (length(fixed_names) == 0) {
    log_warn("No fixed effects matched for {model_id}; check parameter names")
    return(invisible(NULL))
  }

  fixed_draws <- draws_df %>%
    select(all_of(fixed_names), .chain, .iteration, .draw)

  fixed_draws_path <- file.path(out_dir, paste0(model_id, "_fixed_draws.csv"))
  write_csv(fixed_draws, fixed_draws_path)
  log_info("Wrote {fixed_draws_path} ({nrow(fixed_draws)} draws x {ncol(fixed_draws)-3} params)")

  fixed_summary <- summarise_draws(
    select(fixed_draws, all_of(fixed_names)),
    mean, sd, ~quantile(.x, c(0.025, 0.5, 0.975)),
    rhat, ess_bulk, ess_tail
  )
  fixed_summary_path <- file.path(out_dir, paste0(model_id, "_fixed_summary.csv"))
  write_csv(fixed_summary, fixed_summary_path)
  log_info("Wrote {fixed_summary_path}")

  invisible(fixed_summary)
}

extract_species_intercepts <- function(fit, model_id, out_dir, spcd_lookup) {
  ## Look for parameters of form b0_sp[i] or z_sp[i] or alpha_sp[i]
  draws_df <- as_draws_df(fit)
  pn <- variables(draws_df)
  is_sp <- grepl("^(b0_sp|z_sp|alpha_sp|sp_int)\\[[0-9]+\\]$", pn)

  if (!any(is_sp)) {
    log_warn("No species random intercepts found for {model_id}")
    return(invisible(NULL))
  }

  sp_summary <- draws_df %>%
    select(all_of(pn[is_sp])) %>%
    summarise_draws(mean, sd, rhat) %>%
    mutate(
      idx = as.integer(gsub(".*\\[(\\d+)\\].*", "\\1", variable))
    )

  ## Map index to SPCD via spcd_lookup if provided
  if (!is.null(spcd_lookup)) {
    sp_summary <- sp_summary %>%
      left_join(spcd_lookup, by = "idx")
  }

  sp_path <- file.path(out_dir, paste0(model_id, "_species_intercepts.csv"))
  write_csv(sp_summary, sp_path)
  log_info("Wrote {sp_path} ({nrow(sp_summary)} species)")

  invisible(sp_summary)
}

extract_ecodiv_intercepts <- function(fit, model_id, out_dir, ecodiv_lookup) {
  draws_df <- as_draws_df(fit)
  pn <- variables(draws_df)
  is_eco <- grepl("^(b0_eco|z_eco|alpha_eco|eco_int)\\[[0-9]+\\]$", pn)

  if (!any(is_eco)) {
    log_warn("No ecodivision random intercepts found for {model_id}")
    return(invisible(NULL))
  }

  eco_summary <- draws_df %>%
    select(all_of(pn[is_eco])) %>%
    summarise_draws(mean, sd, rhat) %>%
    mutate(
      idx = as.integer(gsub(".*\\[(\\d+)\\].*", "\\1", variable))
    )

  if (!is.null(ecodiv_lookup)) {
    eco_summary <- eco_summary %>%
      left_join(ecodiv_lookup, by = "idx")
  }

  eco_path <- file.path(out_dir, paste0(model_id, "_ecodiv_intercepts.csv"))
  write_csv(eco_summary, eco_path)
  log_info("Wrote {eco_path} ({nrow(eco_summary)} ecodivisions)")

  invisible(eco_summary)
}

extract_modifier_summary <- function(component_dir, base_model, model_id, out_dir) {
  modifier_dir <- file.path(component_dir, "modifier")
  if (!dir.exists(modifier_dir)) {
    log_info("No modifier dir for {model_id}; skipping")
    return(invisible(NULL))
  }

  modifier_summary <- map_dfr(LAMBDAS_TO_SUMMARIZE, function(lam) {
    pattern <- sprintf("%s_modifier_lambda%d_(meta|fit)\\.rds$", base_model, lam)
    candidates <- list.files(modifier_dir, pattern = pattern, full.names = TRUE)
    fit_path <- candidates[grepl("_fit\\.rds$", candidates)][1]
    if (is.na(fit_path) || !file.exists(fit_path)) {
      log_warn("Modifier lambda={lam} fit missing for {model_id}")
      return(tibble())
    }
    fit <- readRDS(fit_path)
    draws <- as_draws_df(fit)
    pn <- variables(draws)
    coef_names <- pn[grepl("^(beta|gamma|mod)_", pn)]
    if (length(coef_names) == 0) return(tibble(lambda = lam))
    summary <- draws %>%
      select(all_of(coef_names)) %>%
      summarise_draws(mean, sd, rhat) %>%
      mutate(lambda = lam)
    summary
  })

  if (nrow(modifier_summary) > 0) {
    out_path <- file.path(out_dir, paste0(model_id, "_modifier_summary.csv"))
    write_csv(modifier_summary, out_path)
    log_info("Wrote {out_path} ({nrow(modifier_summary)} rows across lambdas)")
  }

  invisible(modifier_summary)
}

# ---------------------------------------------------------------------------
# Per-component driver
# ---------------------------------------------------------------------------

extract_component <- function(component, model = NULL,
                              spcd_lookup_fallback = NULL,
                              ecodiv_lookup_fallback = NULL) {
  if (is.null(model)) model <- DEFAULT_MODELS[[component]]
  if (is.null(model)) stop("Unknown component: ", component)

  component_dir <- file.path(conus_dir, component)
  fit_path <- file.path(component_dir, paste0(model, "_fit.rds"))

  if (!file.exists(fit_path)) {
    log_error("Missing fit file: {fit_path}")
    return(invisible(FALSE))
  }

  ## Prefer meta.rds lookups (per-component, current); fall back to global CSVs
  meta_lookups <- load_meta_lookups(component_dir, model)
  spcd_lookup   <- meta_lookups$spcd %||% spcd_lookup_fallback
  ecodiv_lookup <- meta_lookups$ecodiv %||% ecodiv_lookup_fallback

  if (is.null(spcd_lookup))   log_warn("No SPCD lookup; species output will use raw Stan index for {model}")
  if (is.null(ecodiv_lookup)) log_warn("No ecodiv lookup; ecodiv output will use raw Stan index for {model}")

  log_info("=== Extracting {component} / {model} ===")
  log_info("Loading {fit_path} ...")
  fit <- readRDS(fit_path)

  out_dir <- conus_dir  # write CSVs at top of conus/ for easy discovery

  extract_fixed_effects(fit, model, out_dir)
  extract_species_intercepts(fit, model, out_dir, spcd_lookup)
  extract_ecodiv_intercepts(fit, model, out_dir, ecodiv_lookup)
  extract_modifier_summary(component_dir,
                           sub("_cspi_traits1$", "", model),
                           model, out_dir)

  rm(fit); gc()

  log_info("=== Done {component} / {model} ===\n")
  invisible(TRUE)
}

## Null-coalesce helper (since base R has no %||%)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Lookup tables
#
# Map the integer Stan index in z_sp[i] / z_eco[i] back to the real codes
# (SPCD, ecodivision). Two sources, in priority order:
#
# 1. Per-component meta.rds — the fitting scripts (31, 32, 33, 34, ...) save
#    a list that includes the sp_levels and ecodiv_levels vectors. This is
#    the source of truth when present.
#
# 2. Fallback CSVs at calibration/data/conus_{species,ecodiv}_index.csv —
#    only used if the meta.rds doesn't include the lookup vectors.
#
# The meta.rds layout per the existing fitting scripts is approximately:
#   meta$prep$sp_levels   # integer vector of sorted unique SPCDs
#   meta$prep$ecodiv_levels  # character vector of sorted unique ecodiv codes
# Some older scripts use meta$species / meta$ecodivs instead; we try both.
# ---------------------------------------------------------------------------

load_meta_lookups <- function(component_dir, model_id) {
  meta_path <- file.path(component_dir, paste0(model_id, "_meta.rds"))
  if (!file.exists(meta_path)) {
    log_info("No meta.rds at {meta_path}")
    return(list(spcd = NULL, ecodiv = NULL))
  }
  meta <- readRDS(meta_path)

  ## Walk all known nesting conventions used by the fitting scripts. The
  ## cspi_traits1 round stores lookups under meta$prep_meta$sp / $L1 / $L2 / $L3
  ## (verified via SLURM job 9357983 on 2026-05-11). Earlier scripts may use
  ## meta$prep$sp_levels or top-level meta$species. Try each in order.
  pull_field <- function(obj, names) {
    nests <- list(obj, obj$prep_meta, obj$prep, obj$data)
    for (nest in nests) {
      if (is.null(nest)) next
      for (n in names) if (!is.null(nest[[n]])) return(nest[[n]])
    }
    NULL
  }

  spcd_vec <- pull_field(meta, c("sp", "sp_levels", "species", "spcd_levels", "SPCD"))
  ## Ecodivision RE is absent in cspi_traits1 fits (they use L1/L2/L3 nested
  ## plot/landscape grouping instead). Keep ecodiv lookup wired for older fits.
  eco_vec  <- pull_field(meta, c("ecodiv_levels", "ecodivs", "ecodivision_levels", "ecodiv"))

  spcd_lookup <- if (!is.null(spcd_vec)) {
    tibble(idx = seq_along(spcd_vec), SPCD = as.integer(spcd_vec))
  } else NULL

  ecodiv_lookup <- if (!is.null(eco_vec)) {
    tibble(idx = seq_along(eco_vec), ecodiv = as.character(eco_vec))
  } else NULL

  list(spcd = spcd_lookup, ecodiv = ecodiv_lookup)
}

load_spcd_lookup_fallback <- function() {
  path <- file.path(calibration_dir, "data", "conus_species_index.csv")
  if (file.exists(path)) {
    read_csv(path, col_types = cols())
  } else {
    NULL
  }
}

load_ecodiv_lookup_fallback <- function() {
  path <- file.path(calibration_dir, "data", "conus_ecodiv_index.csv")
  if (file.exists(path)) {
    read_csv(path, col_types = cols())
  } else {
    NULL
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  ## Fallbacks for cases where meta.rds doesn't include the lookup vectors
  spcd_fallback <- load_spcd_lookup_fallback()
  ecodiv_fallback <- load_ecodiv_lookup_fallback()

  if (isTRUE(args$all)) {
    components <- names(DEFAULT_MODELS)
  } else if (!is.null(args$component)) {
    components <- args$component
  } else {
    stop("Specify --component <name> or --all")
  }

  results <- map_lgl(components, function(c) {
    extract_component(c, model = args$model,
                      spcd_lookup_fallback = spcd_fallback,
                      ecodiv_lookup_fallback = ecodiv_fallback)
  })

  ok <- sum(results)
  log_info("Extraction complete: {ok}/{length(components)} components ok")
}

if (!interactive()) main()
