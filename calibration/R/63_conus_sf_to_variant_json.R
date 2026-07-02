#!/usr/bin/env Rscript
# =============================================================================
# 63_conus_sf_to_variant_json.R
#
# Land CONUS species-free (Leg B, trait-driven) component fits into per-variant
# config JSON files under a new `categories_conus_sf.{component}` block. This is
# the Leg B counterpart to 62_conus_to_variant_json.R (Leg A, per-species).
#
# Where Leg A gives every species its own random intercept, Leg B replaces the
# per-species intercept with a trait fixed-effect block (standardized traits
# times gamma) plus nested ecoregion random effects (L1 > L2 > L3) and, where
# fit, a forest-type random effect. The result generalizes to species that have
# no per-species fit, which is the point of the "species-independent equations"
# option in the CONUS single variant.
#
# Unlike Leg A, this script does the extract and the serialize in one pass: the
# species-free blocks are compact (one trait_effect per species, small RE
# tables, no multi-GB draws CSV), so there is no separate 61-style extractor.
#
# Reads (per component, on Cardinal where the fits live):
#   {OUT_DIR}/{model}_fit.rds        cmdstanr fit object
#   {OUT_DIR}/{model}_meta.rds       prep_meta (sp, L1, L2, L3), trait_cols, ...
#   calibration/traits/species_traits.rds   raw species trait table
#   config/calibrated/{variant}.json (existing, must already carry categories)
#
# Writes:
#   config/calibrated/{variant}.json with categories_conus_sf.{component} block,
#   in the exact shape config_loader.py get_conus_sf_runtime_block() decodes:
#     model
#     fixed_effects     {param, mean, sd}
#     trait_gamma       {trait_col, gamma_mean, scale_mean, scale_sd}
#     species           {SPCD, trait_effect_mean, raw_<col>..., std_<col>...}
#     re_L1 / re_L2 / re_L3 / re_FT   {level, mean}
#     hybrid_source_map {SPCD, source}   leg_a where a reliable per-species fit
#                                        exists in categories_conus, else leg_b
#
# Usage:
#   Rscript 63_conus_sf_to_variant_json.R --variant all --component dg
#   Rscript 63_conus_sf_to_variant_json.R --variant ne  --component all --dry-run
#   Rscript 63_conus_sf_to_variant_json.R --variant all --component dg \
#       --outdir calibration/output/conus/dg/speciesfree_pilot \
#       --model  dg_kuehne_cspi_traits1_b1
#
# Notes
# -----
# - The trait standardization is recomputed here from species_traits.rds using
#   the SAME logic as the fitting script (median-impute NA, then z-score over
#   the fit's species set). This must stay in lockstep with
#   32_fit_dg_kuehne_speciesfree.R; if the fit changes its standardization,
#   change it here too. scale_mean / scale_sd are stored so the loader can
#   standardize traits for species outside the fit at runtime.
# - Existing categories.* and categories_conus.* blocks are left untouched.
# - Only components whose species-free fit is actually present are written; a
#   missing fit is logged and skipped (Leg B is being rolled out component by
#   component, DG first).
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
traits_file     <- file.path(calibration_dir, "traits", "species_traits.rds")

# Component -> (default speciesfree output dir, default model name). Override
# either via --outdir / --model. DG is the landed pilot; the others follow the
# same OUT_DIR / OUT_NAME convention as their fitting scripts once fit.
SF_DEFAULTS <- list(
  diameter_growth = list(
    outdir = file.path(conus_dir, "dg", "speciesfree_pilot"),
    model  = "dg_kuehne_cspi_traits1_b1"),
  height_growth = list(
    outdir = file.path(conus_dir, "hg", "speciesfree"),
    model  = "hg_organon_fixedK_cspi_traits1_sf"),
  height_diameter = list(
    outdir = file.path(conus_dir, "ht_dbh", "speciesfree"),
    model  = "htdbh_wykoff_lognormal_cspi_traits1_sf"),
  height_crown_base = list(
    outdir = file.path(conus_dir, "hcb", "speciesfree"),
    model  = "hcb_organon_cspi_traits1_sf"),
  mortality = list(
    outdir = file.path(conus_dir, "mortality", "speciesfree"),
    model  = "mort_logit_simple_cspi_traits1_sf"),
  crown_recession = list(
    outdir = file.path(conus_dir, "crown_recession", "speciesfree"),
    model  = "cr_recession_cspi_traits1_sf")
)

COMPONENT_ALIASES <- list(
  dg = "diameter_growth", hg = "height_growth", ht_dbh = "height_diameter",
  hcb = "height_crown_base", mort = "mortality", mortality = "mortality",
  cr = "crown_recession", crown_recession = "crown_recession"
)

ALL_VARIANTS <- c("acd","ak","bm","ca","ci","cr","cs","ec","em","ie","kt",
                  "ls","nc","ne","oc","on","op","pn","sn","so","tt","ut",
                  "wc","ws","bc")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  out <- list(variant = NULL, component = NULL, model = NULL,
              outdir = NULL, dry_run = FALSE)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--variant"   && i < length(args)) { out$variant   <- args[i+1]; i <- i+2 }
    else if (a == "--component" && i < length(args)) { out$component <- args[i+1]; i <- i+2 }
    else if (a == "--model" && i < length(args)) { out$model     <- args[i+1]; i <- i+2 }
    else if (a == "--outdir" && i < length(args)) { out$outdir   <- args[i+1]; i <- i+2 }
    else if (a == "--dry-run") { out$dry_run <- TRUE; i <- i+1 }
    else i <- i+1
  }
  out
}

resolve_component <- function(name) {
  if (name == "all") return(names(SF_DEFAULTS))
  if (!is.null(COMPONENT_ALIASES[[name]])) return(COMPONENT_ALIASES[[name]])
  if (name %in% names(SF_DEFAULTS)) return(name)
  stop("Unknown component: ", name)
}

resolve_variant <- function(name) {
  if (name == "all") return(ALL_VARIANTS)
  if (name %in% ALL_VARIANTS) return(name)
  stop("Unknown variant: ", name)
}

# ---------------------------------------------------------------------------
# Trait standardization (must match 32_fit_dg_kuehne_speciesfree.R lines 148-161)
# ---------------------------------------------------------------------------

build_trait_scaling <- function(sp_levels, trait_cols) {
  if (!file.exists(traits_file)) {
    stop("species_traits.rds not found at ", traits_file,
         " — required to recompute trait standardization for Leg B.")
  }
  traits <- as.data.frame(readRDS(traits_file))
  sub <- traits[match(sp_levels, traits$SPCD), c("SPCD", trait_cols), drop = FALSE]
  W_raw <- as.matrix(sub[, trait_cols, drop = FALSE])

  scale_mean <- numeric(length(trait_cols))
  scale_sd   <- numeric(length(trait_cols))
  W_std <- W_raw
  for (j in seq_along(trait_cols)) {
    col <- W_raw[, j]
    na  <- is.na(col)
    if (any(na)) col[na] <- median(col[!na], na.rm = TRUE)  # median-impute
    m <- mean(col); s <- sd(col)
    scale_mean[j] <- m
    scale_sd[j]   <- s
    W_std[, j] <- (col - m) / s
    W_raw[, j] <- col  # keep the imputed raw value so raw_/std_ stay consistent
  }
  list(SPCD = sp_levels, trait_cols = trait_cols,
       scale_mean = scale_mean, scale_sd = scale_sd,
       W_raw = W_raw, W_std = W_std)
}

# ---------------------------------------------------------------------------
# Fit readers
# ---------------------------------------------------------------------------

# Summarise a set of posterior variables to mean/sd. Uses the fit's own
# $summary() so this works on any cmdstanr fit without loading all draws.
summ <- function(fit, vars) {
  s <- fit$summary(variables = vars, "mean", "sd")
  tibble(variable = s$variable, mean = s$mean, sd = s$sd)
}

read_sf_component <- function(component, outdir, model) {
  fit_path  <- file.path(outdir, paste0(model, "_fit.rds"))
  meta_path <- file.path(outdir, paste0(model, "_meta.rds"))
  if (!file.exists(fit_path) || !file.exists(meta_path)) {
    log_warn("Species-free fit/meta missing for {component}: {fit_path}")
    return(NULL)
  }
  log_info("Reading {component} SF fit: {model}")
  fit  <- readRDS(fit_path)
  meta <- readRDS(meta_path)

  sp_levels  <- meta$prep_meta$sp
  L1_levels  <- meta$prep_meta$L1 %||% character(0)
  L2_levels  <- meta$prep_meta$L2 %||% character(0)
  L3_levels  <- meta$prep_meta$L3 %||% character(0)
  FT_levels  <- meta$prep_meta$FT %||% character(0)
  trait_cols <- meta$trait_cols

  # Fixed effects: intercept + covariate coefficients + sigmas. Grab everything
  # that is not indexed and not the per-species trait_effect vector.
  all_vars <- fit$metadata()$model_params %||% fit$summary()$variable
  fe_vars <- all_vars[!grepl("\\[", all_vars) &
                      !all_vars %in% c("lp__") &
                      !grepl("^(z_|trait_effect|log_lik|mu_|y_rep)", all_vars)]
  fixed <- summ(fit, fe_vars)

  # Trait gamma vector, aligned to trait_cols order (gamma[1]..gamma[P]).
  P <- length(trait_cols)
  gamma <- summ(fit, paste0("gamma[", seq_len(P), "]"))$mean

  # Per-species trait effect (posterior mean of W_std %*% gamma).
  te <- summ(fit, paste0("trait_effect[", seq_along(sp_levels), "]"))$mean

  # Nested ecoregion random-effect means. Level order matches *_levels.
  re_means <- function(prefix, levels) {
    if (length(levels) == 0) return(tibble(level = character(0), mean = numeric(0)))
    v <- paste0(prefix, "[", seq_along(levels), "]")
    tibble(level = as.character(levels), mean = summ(fit, v)$mean)
  }
  re_L1 <- re_means("z_L1", L1_levels)
  re_L2 <- re_means("z_L2", L2_levels)
  re_L3 <- re_means("z_L3", L3_levels)
  re_FT <- re_means("z_FT", FT_levels)

  scaling <- build_trait_scaling(sp_levels, trait_cols)

  list(model = model, sp = sp_levels, trait_cols = trait_cols,
       fixed = fixed, gamma = gamma, trait_effect = te,
       re_L1 = re_L1, re_L2 = re_L2, re_L3 = re_L3, re_FT = re_FT,
       scaling = scaling)
}

# ---------------------------------------------------------------------------
# Block builder (per variant, per component)
# ---------------------------------------------------------------------------

build_sf_block <- function(sf, variant_spcds, legA_spcds) {
  # Keep the national species table but flag which SPCDs the variant actually
  # uses via the hybrid source map; runtime looks up by SPCD, not by position.
  sp <- sf$sp
  sc <- sf$scaling

  species_block <- list(SPCD = as.integer(sp),
                        trait_effect_mean = as.numeric(sf$trait_effect))
  for (j in seq_along(sf$trait_cols)) {
    col <- sf$trait_cols[j]
    species_block[[paste0("raw_", col)]] <- as.numeric(sc$W_raw[, j])
    species_block[[paste0("std_", col)]] <- as.numeric(sc$W_std[, j])
  }

  # leg_a where the variant has a reliable per-species fit for this component
  # (SPCD present in categories_conus.{component}.species_intercepts), else leg_b.
  src <- ifelse(sp %in% legA_spcds, "leg_a", "leg_b")
  # Restrict the source map to species the variant actually carries, when known.
  if (!is.null(variant_spcds)) {
    keep <- sp %in% variant_spcds
    hybrid_map <- list(SPCD = as.integer(sp[keep]), source = src[keep])
  } else {
    hybrid_map <- list(SPCD = as.integer(sp), source = src)
  }

  re_to_list <- function(re) list(level = as.character(re$level),
                                  mean = as.numeric(re$mean))

  list(
    model = sf$model,
    fixed_effects = list(param = sf$fixed$variable,
                         mean  = as.numeric(sf$fixed$mean),
                         sd    = as.numeric(sf$fixed$sd)),
    trait_gamma = list(trait_col  = sf$trait_cols,
                       gamma_mean = as.numeric(sf$gamma),
                       scale_mean = as.numeric(sc$scale_mean),
                       scale_sd   = as.numeric(sc$scale_sd)),
    species = species_block,
    re_L1 = re_to_list(sf$re_L1),
    re_L2 = re_to_list(sf$re_L2),
    re_L3 = re_to_list(sf$re_L3),
    re_FT = re_to_list(sf$re_FT),
    hybrid_source_map = hybrid_map,
    notes = "Species-free (Leg B): species effect = standardized traits x gamma."
  )
}

# ---------------------------------------------------------------------------
# Per-variant integration
# ---------------------------------------------------------------------------

integrate_variant <- function(variant, sf_by_component, dry_run) {
  json_path <- file.path(calibrated_dir, paste0(variant, ".json"))
  if (!file.exists(json_path)) {
    log_error("Variant config not found: {json_path}"); return(invisible(FALSE))
  }
  cfg <- fromJSON(json_path, simplifyVector = FALSE)

  variant_spcds <- cfg$categories$species_definitions$FIAJSP
  if (!is.null(variant_spcds)) {
    variant_spcds <- as.integer(unlist(variant_spcds))
    variant_spcds <- variant_spcds[!is.na(variant_spcds)]
  }

  if (is.null(cfg$categories_conus_sf)) cfg$categories_conus_sf <- list()

  for (component in names(sf_by_component)) {
    sf <- sf_by_component[[component]]
    if (is.null(sf)) next

    # Leg A species set for this variant/component (for the hybrid source map).
    legA <- cfg$categories_conus[[component]]$species_intercepts$SPCD
    legA_spcds <- if (!is.null(legA)) as.integer(unlist(legA)) else integer(0)

    block <- build_sf_block(sf, variant_spcds, legA_spcds)
    cfg$categories_conus_sf[[component]] <- block
    log_info("  {variant}/{component}: {length(sf$sp)} species, {length(sf$trait_cols)} traits, L1/L2/L3 = {nrow(sf$re_L1)}/{nrow(sf$re_L2)}/{nrow(sf$re_L3)}")
  }

  cfg$categories_conus_sf$metadata <- list(
    integration_date = format(Sys.Date()),
    pipeline_version = "fvs-conus phase 4 species-free (Leg B)",
    components_present = setdiff(names(cfg$categories_conus_sf), "metadata")
  )

  if (dry_run) { log_info("[dry-run] would write {json_path}"); return(invisible(TRUE)) }

  backup <- paste0(json_path, ".pre_sf_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(json_path, backup)
  write_json(cfg, json_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  log_info("Wrote {json_path} (backup {backup})")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$variant) || is.null(args$component)) {
    stop("Usage: --variant <name|all> --component <name|all> [--outdir D] [--model M] [--dry-run]")
  }

  variants   <- resolve_variant(args$variant)
  components <- resolve_component(args$component)

  # Read each component's SF fit once (fits are national; variants only differ
  # in which species they carry and their Leg A coverage for the source map).
  sf_by_component <- list()
  for (component in components) {
    d <- SF_DEFAULTS[[component]]
    outdir <- args$outdir %||% d$outdir
    model  <- if (!is.null(args$model) && length(components) == 1) args$model else d$model
    sf_by_component[[component]] <- read_sf_component(component, outdir, model)
  }
  if (all(map_lgl(sf_by_component, is.null))) {
    stop("No species-free fits found for the requested component(s). ",
         "Check --outdir / --model or run the Leg B fits first.")
  }

  log_info("Landing Leg B for {length(components)} component(s) into {length(variants)} variant(s)")
  ok <- 0
  for (v in variants) {
    if (integrate_variant(v, sf_by_component, args$dry_run)) ok <- ok + 1
  }
  log_info("{ok}/{length(variants)} variants updated")
}

if (!interactive()) main()
