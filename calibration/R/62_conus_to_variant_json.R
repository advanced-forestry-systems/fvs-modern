#!/usr/bin/env Rscript
# =============================================================================
# 62_conus_to_variant_json.R
#
# Land CONUS Phase 4 calibration parameters into per-variant config JSON files
# under a new `categories_conus.{component}` block. Reads the summary CSVs
# produced by 61_extract_conus_summaries.R, filters to the variant's species
# and ecodivision coverage, and writes the merged JSON back.
#
# Reads:
#   calibration/output/conus/{model}_fixed_summary.csv
#   calibration/output/conus/{model}_species_intercepts.csv
#   calibration/output/conus/{model}_ecodiv_intercepts.csv
#   calibration/output/conus/{model}_modifier_summary.csv
#   config/calibrated/{variant}.json (existing)
#
# Writes:
#   config/calibrated/{variant}.json with categories_conus.{component} block
#
# Usage:
#   Rscript 62_conus_to_variant_json.R --variant ne --component dg
#   Rscript 62_conus_to_variant_json.R --variant all --component all
#   Rscript 62_conus_to_variant_json.R --variant ne --component all --modifier-lambda 10
#
# Notes
# -----
# - Existing categories.* blocks are preserved unchanged. categories_conus is
#   the new alternate parameter set; the runtime layer in config_loader.py
#   chooses between them.
# - Posterior draws are NOT embedded in the JSON; the JSON references the
#   draws CSV by relative path. Runtime samplers read the CSV on demand.
# - Species intercepts are filtered to the variant's FIAJSP set. SPCDs that
#   do not appear in the CONUS posterior are flagged in `species_missing`.
# - Ecodivision intercepts are filtered by the variant's ecodivision coverage
#   list (variant_ecodiv_coverage.csv). Coverage weights are stored alongside
#   so the runtime sampler can compute a coverage-weighted ecodiv RE.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(logger)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

project_root    <- Sys.getenv("FVS_PROJECT_ROOT", normalizePath(".."))
calibration_dir <- file.path(project_root, "calibration")
conus_dir       <- file.path(calibration_dir, "output", "conus")
config_dir      <- file.path(project_root, "config")
calibrated_dir  <- file.path(config_dir, "calibrated")

# Component -> model mapping (matches DEFAULT_MODELS in 61_extract_conus_summaries.R)
DEFAULT_MODELS <- list(
  diameter_growth   = "dg_kuehne_cspi_traits1",
  height_growth     = "hg_organon_fixedK_cspi_traits1",
  height_diameter   = "htdbh_wykoff_lognormal_cspi_traits1",
  height_crown_base = "hcb_organon_cspi_traits1",
  mortality         = "mort_logit_simple_cspi_traits1",
  crown_recession   = "cr_recession_cspi_traits1"
)

# Aliases so users can pass short component names
COMPONENT_ALIASES <- list(
  dg              = "diameter_growth",
  hg              = "height_growth",
  ht_dbh          = "height_diameter",
  hcb             = "height_crown_base",
  mort            = "mortality",
  mortality       = "mortality",
  cr              = "crown_recession",
  crown_recession = "crown_recession"
)

DEFAULT_LAMBDA <- 10  # selected from 5/10/20 sweep

ALL_VARIANTS <- c("acd","ak","bm","ca","ci","cr","cs","ec","em","ie","kt",
                  "ls","nc","ne","oc","on","op","pn","sn","so","tt","ut",
                  "wc","ws","bc")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  out <- list(
    variant = NULL, component = NULL, model = NULL,
    modifier_lambda = DEFAULT_LAMBDA, dry_run = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--variant" && i < length(args)) {
      out$variant <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--component" && i < length(args)) {
      out$component <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--model" && i < length(args)) {
      out$model <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--modifier-lambda" && i < length(args)) {
      out$modifier_lambda <- as.integer(args[i + 1]); i <- i + 2
    } else if (args[i] == "--dry-run") {
      out$dry_run <- TRUE; i <- i + 1
    } else {
      i <- i + 1
    }
  }
  out
}

resolve_component <- function(name) {
  if (name == "all") return(names(DEFAULT_MODELS))
  if (!is.null(COMPONENT_ALIASES[[name]])) return(COMPONENT_ALIASES[[name]])
  if (name %in% names(DEFAULT_MODELS)) return(name)
  stop("Unknown component: ", name)
}

resolve_variant <- function(name) {
  if (name == "all") return(ALL_VARIANTS)
  if (name %in% ALL_VARIANTS) return(name)
  stop("Unknown variant: ", name)
}

# ---------------------------------------------------------------------------
# CSV loaders (from 61 outputs)
# ---------------------------------------------------------------------------

load_fixed_summary <- function(model_id) {
  path <- file.path(conus_dir, paste0(model_id, "_fixed_summary.csv"))
  if (!file.exists(path)) {
    log_error("Missing fixed summary: {path}. Run 61_extract_conus_summaries.R first.")
    return(NULL)
  }
  read_csv(path, col_types = cols())
}

load_species_intercepts <- function(model_id) {
  path <- file.path(conus_dir, paste0(model_id, "_species_intercepts.csv"))
  if (!file.exists(path)) {
    log_error("Missing species intercepts: {path}")
    return(NULL)
  }
  read_csv(path, col_types = cols())
}

load_ecodiv_intercepts <- function(model_id) {
  ## Ecodivision random intercepts are model-design optional. Several CONUS
  ## Phase 4 fits (e.g., dg_kuehne_cspi_traits1, hg_organon_fixedK_cspi_traits1)
  ## omit the ecodivision RE entirely because species + climate covariates
  ## absorbed the regional variance. Missing CSV is therefore not fatal; we
  ## return an empty tibble that the block builder converts to an empty
  ## ecodiv_intercepts block.
  path <- file.path(conus_dir, paste0(model_id, "_ecodiv_intercepts.csv"))
  if (!file.exists(path)) {
    log_info("No ecodiv intercepts at {path}; model has no ecodivision RE")
    return(tibble::tibble(variable = character(), mean = numeric(),
                          sd = numeric(), rhat = numeric(),
                          idx = integer(), ecodiv = character()))
  }
  read_csv(path, col_types = cols())
}

load_modifier_summary <- function(model_id, lambda) {
  path <- file.path(conus_dir, paste0(model_id, "_modifier_summary.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, col_types = cols()) %>%
    filter(lambda == !!lambda)
}

load_variant_ecodiv_coverage <- function(variant) {
  path <- file.path(calibration_dir, "data", "variant_ecodiv_coverage.csv")
  if (!file.exists(path)) {
    log_warn("Variant ecodiv coverage table not found at {path}; using all ecodivisions equally weighted")
    return(NULL)
  }
  read_csv(path, col_types = cols()) %>%
    filter(variant == !!variant)
}

# ---------------------------------------------------------------------------
# Block builder
# ---------------------------------------------------------------------------

build_conus_block <- function(component, model_id, lambda,
                              fixed_summary, species_int, ecodiv_int,
                              modifier_summary,
                              variant_spcds, variant_ecodiv_cov) {

  # Fixed effects: store mean / sd / quantile summary by parameter name
  fixed_block <- fixed_summary %>%
    transmute(
      param  = variable,
      mean   = mean,
      sd     = sd,
      q025   = `2.5%`,
      median = `50%`,
      q975   = `97.5%`,
      rhat   = if ("rhat" %in% names(.)) rhat else NA_real_
    ) %>%
    as.list()

  # Species intercepts filtered to variant's species set
  if ("SPCD" %in% names(species_int) && !is.null(variant_spcds)) {
    species_filtered <- species_int %>% filter(SPCD %in% variant_spcds)
    missing_spcds <- setdiff(variant_spcds, species_int$SPCD)
  } else {
    species_filtered <- species_int
    missing_spcds <- integer(0)
  }
  species_block <- species_filtered %>%
    transmute(
      SPCD = if ("SPCD" %in% names(.)) SPCD else idx,
      mean = mean,
      sd   = sd
    ) %>%
    as.list()

  # Ecodivision intercepts: keep all, store coverage weights alongside if known
  ecodiv_block <- ecodiv_int %>%
    transmute(
      ecodiv = if ("ecodiv" %in% names(.)) ecodiv else as.character(idx),
      mean   = mean,
      sd     = sd
    ) %>%
    as.list()

  ecodiv_weights <- if (!is.null(variant_ecodiv_cov)) {
    setNames(as.list(variant_ecodiv_cov$weight), variant_ecodiv_cov$ecodiv)
  } else {
    list()
  }

  # Modifier summary at chosen lambda
  modifier_block <- if (!is.null(modifier_summary) && nrow(modifier_summary) > 0) {
    modifier_summary %>%
      transmute(coef = variable, mean = mean, sd = sd) %>%
      as.list()
  } else {
    list()
  }

  list(
    model           = model_id,
    modifier_lambda = lambda,
    fixed_effects   = fixed_block,
    species_intercepts = species_block,
    species_missing = as.integer(missing_spcds),
    ecodiv_intercepts = ecodiv_block,
    ecodiv_weights   = ecodiv_weights,
    modifier         = modifier_block,
    draws_csv        = paste0("calibration/output/conus/", model_id, "_fixed_draws.csv"),
    notes            = "Posterior summary; full draws referenced via draws_csv."
  )
}

# ---------------------------------------------------------------------------
# Per-variant integration
# ---------------------------------------------------------------------------

integrate_variant <- function(variant, components, model_overrides,
                              modifier_lambda, dry_run) {
  json_path <- file.path(calibrated_dir, paste0(variant, ".json"))
  if (!file.exists(json_path)) {
    log_error("Variant config not found: {json_path}")
    return(invisible(FALSE))
  }
  cfg <- fromJSON(json_path, simplifyVector = FALSE)

  variant_spcds <- cfg$categories$species_definitions$FIAJSP
  if (is.null(variant_spcds)) {
    log_warn("Variant {variant} has no FIAJSP species list; species filtering disabled")
    variant_spcds <- NULL
  } else {
    variant_spcds <- as.integer(variant_spcds[!is.na(variant_spcds) & variant_spcds != "NA"])
  }

  variant_ecodiv_cov <- load_variant_ecodiv_coverage(variant)

  if (is.null(cfg$categories_conus)) cfg$categories_conus <- list()
  cfg$categories_conus$metadata <- list(
    integration_date = format(Sys.Date()),
    pipeline_version = "fvs-conus phase 4",
    selected_lambda  = modifier_lambda,
    components_present = c(names(cfg$categories_conus), components) %>% unique()
  )

  for (component in components) {
    model_id <- if (!is.null(model_overrides[[component]])) {
      model_overrides[[component]]
    } else {
      DEFAULT_MODELS[[component]]
    }
    log_info("--- {variant} / {component} ({model_id}) ---")

    fs <- load_fixed_summary(model_id)
    si <- load_species_intercepts(model_id)
    ei <- load_ecodiv_intercepts(model_id)
    ## ei may be an empty tibble for models with no ecodivision RE — that's
    ## allowed. Only fs and si are mandatory.
    if (is.null(fs) || is.null(si) || is.null(ei)) {
      log_error("Skipping {component} for {variant} due to missing summaries")
      next
    }
    ms <- load_modifier_summary(model_id, modifier_lambda)

    block <- build_conus_block(
      component, model_id, modifier_lambda,
      fs, si, ei, ms,
      variant_spcds, variant_ecodiv_cov
    )
    cfg$categories_conus[[component]] <- block

    log_info("  fixed_effects: {length(block$fixed_effects$param)} params; species kept: {length(block$species_intercepts$SPCD)}; species missing: {length(block$species_missing)}; ecodivs: {length(block$ecodiv_intercepts$ecodiv)}")
  }

  if (dry_run) {
    log_info("[dry-run] would write {json_path}")
    return(invisible(TRUE))
  }

  ## Backup before writing
  backup_path <- paste0(json_path, ".pre_conus_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(json_path, backup_path)
  write_json(cfg, json_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  log_info("Wrote {json_path} (backup at {backup_path})")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  if (is.null(args$variant) || is.null(args$component)) {
    stop("Usage: --variant <name|all> --component <name|all> [--modifier-lambda N] [--model M] [--dry-run]")
  }

  variants   <- resolve_variant(args$variant)
  components <- resolve_component(args$component)
  model_ov   <- list()
  if (!is.null(args$model) && length(components) == 1) {
    model_ov[[components[1]]] <- args$model
  }

  log_info("Integrating {length(components)} component(s) into {length(variants)} variant(s) at lambda={args$modifier_lambda}")

  ok <- 0
  for (v in variants) {
    if (integrate_variant(v, components, model_ov, args$modifier_lambda, args$dry_run)) {
      ok <- ok + 1
    }
  }
  log_info("{ok}/{length(variants)} variants updated")
}

if (!interactive()) main()
