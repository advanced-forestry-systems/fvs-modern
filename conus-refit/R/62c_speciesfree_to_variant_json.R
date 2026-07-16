#!/usr/bin/env Rscript
# =============================================================================
# 62b_speciesfree_to_variant_json.R
#
# Species-free (Leg B) companion to 62_conus_to_variant_json.R. Reads a 61b
# bundle for one component and lands a `categories_conus_sf.{component}` block
# into per-variant config JSON, parallel to the existing per-species
# `categories_conus.{component}` (Leg A) block.
#
# Leg A stores per-species intercepts. Leg B stores instead:
#   - fixed_effects        global fixed + covariate coefficients
#   - trait_gamma          trait coefficients + standardization constants
#                          (engine computes a species effect as W_species . gamma)
#   - re_L1/re_L2/re_L3/re_FT   random-effect tables keyed by real codes
#   - species              raw + standardized traits + trait_effect per species
#   - hybrid_source_map    per species: "leg_a" if the species has a reliable
#                          per-species fit in categories_conus (use it), else
#                          "leg_b" (trait fallback). This is Aaron's chosen
#                          hybrid-per-component-per-species default.
#
# SAFETY: defaults to --dry_run=TRUE, which writes {variant}.sf_preview.json
# (a full copy with the new block) instead of modifying the real config. Run
# with --dry_run=FALSE only when all component bundles are ready and you want
# to land into the production configs (a .pre_sf_<timestamp> backup is made).
#
# Usage:
#   Rscript calibration/R/62b_speciesfree_to_variant_json.R \
#     --component hg \
#     --bundle_dir calibration/output/conus/sf_integration \
#     --bundle_prefix hg_v5_prod \
#     --config_dir config/calibrated \
#     --variant ne            # omit for all variants
#     --dry_run TRUE
#
# Author: A. Weiskittel + Claude
# Date: 2026-05-21
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) > 0) return(sub(paste0("^--", name, "="), "", m[1]))
  i <- which(args == paste0("--", name))
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}

COMPONENT    <- get_arg("component")
BUNDLE_DIR   <- get_arg("bundle_dir", "calibration/output/conus/sf_integration")
BUNDLE_PREF  <- get_arg("bundle_prefix")
CONFIG_DIR   <- get_arg("config_dir", "config/calibrated")
VARIANT      <- get_arg("variant", NULL)        # NULL = all
DRY_RUN      <- toupper(get_arg("dry_run", "TRUE")) %in% c("TRUE", "T", "1", "YES")

stopifnot(!is.null(COMPONENT), !is.null(BUNDLE_PREF))

# component short name -> variant-json category key (matches Leg A keys)
KEY_MAP <- c(dg = "diameter_growth", hg = "height_growth",
             htdbh = "height_diameter", hcb = "height_crown_base",
             mort = "mortality", cr = "crown_recession",
             ingrowth = "ingrowth")
CKEY <- KEY_MAP[[COMPONENT]]
if (is.null(CKEY)) stop("Unknown component: ", COMPONENT)

cat("== 62b_speciesfree_to_variant_json.R ==\n")
cat("component:", COMPONENT, "-> category key:", CKEY, "\n")
cat("dry_run:  ", DRY_RUN, "\n\n")

# ---- Read the 61b bundle ---------------------------------------------------
bp <- function(suffix) file.path(BUNDLE_DIR, paste0(BUNDLE_PREF, suffix))
manifest <- fromJSON(paste(readLines(bp("_sf_manifest.json")), collapse="\n"))
fixed    <- fread(bp("_sf_fixed.csv"))
gamma    <- fread(bp("_sf_gamma.csv"))
species  <- fread(bp("_sf_species.csv"))
fixed    <- fixed[variable != "lp__"]   # drop log-posterior if present

read_re <- function(tag) {
  f <- bp(paste0("_sf_re_", tag, ".csv"))
  if (file.exists(f)) fread(f) else NULL
}
re_L1 <- read_re("L1"); re_L2 <- read_re("L2")
re_L3 <- read_re("L3"); re_FT <- read_re("FT")

# Build the reusable (variant-independent) part of the SF block once.
sf_common <- list(
  model = manifest$form,
  stan_file = manifest$stan_file,
  n_obs = manifest$n_obs,
  fixed_effects = list(
    param = fixed$variable, mean = fixed$mean, sd = fixed$sd,
    q5 = fixed$q5, q95 = fixed$q95),
  trait_gamma = list(
    trait_col = gamma$trait_col, gamma_mean = gamma$gamma_mean,
    gamma_sd = gamma$gamma_sd,
    scale_mean = gamma$scale_mean, scale_sd = gamma$scale_sd),
  re_L1 = if (!is.null(re_L1)) list(level = re_L1$level, mean = re_L1$mean) else NULL,
  re_L2 = if (!is.null(re_L2)) list(level = re_L2$level, mean = re_L2$mean) else NULL,
  re_L3 = if (!is.null(re_L3)) list(level = re_L3$level, mean = re_L3$mean) else NULL,
  re_FT = if (!is.null(re_FT)) list(level = re_FT$level, mean = re_FT$mean) else NULL,
  notes = paste0("Species-free (B1) coefficients. Species effect computed at ",
                 "runtime as standardized traits times gamma plus ecoregion / ",
                 "forest-type random effects. Source bundle: ", BUNDLE_PREF, ".")
)

bundle_spcd <- species$SPCD

# ---- Per-variant landing ---------------------------------------------------
variant_files <- if (!is.null(VARIANT)) {
  file.path(CONFIG_DIR, paste0(VARIANT, ".json"))
} else {
  fs <- list.files(CONFIG_DIR, pattern = "\\.json$", full.names = TRUE)
  fs[!grepl("(_draws|pre_conus|pre_sf|sf_preview)", fs)]
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
n_done <- 0L
for (vf in variant_files) {
  if (!file.exists(vf)) { cat("skip (missing):", vf, "\n"); next }
  d <- fromJSON(paste(readLines(vf), collapse="\n"), simplifyVector = TRUE)

  # Leg A per-species reliability for the hybrid source map: a species is
  # "leg_a" if it has a per-species intercept in the existing categories_conus.
  legA_spcd <- integer(0)
  if (!is.null(d$categories_conus) && !is.null(d$categories_conus[[CKEY]]) &&
      !is.null(d$categories_conus[[CKEY]]$species_intercepts)) {
    legA_spcd <- as.integer(d$categories_conus[[CKEY]]$species_intercepts$SPCD)
  }
  # Hybrid map over the union of species either leg knows.
  all_spcd <- sort(unique(c(bundle_spcd, legA_spcd)))
  src <- ifelse(all_spcd %in% legA_spcd, "leg_a", "leg_b")

  sf_block <- sf_common
  # restrict the species trait table to species this variant / either leg uses
  sp_keep <- species[SPCD %in% all_spcd]
  sf_block$species <- as.list(sp_keep)
  sf_block$hybrid_source_map <- list(SPCD = all_spcd, source = src)

  if (is.null(d$categories_conus_sf)) d$categories_conus_sf <- list()
  d$categories_conus_sf[[CKEY]] <- sf_block
  d$categories_conus_sf$metadata <- list(
    pipeline_version = "fvs-conus species-free (Leg B)",
    integration_date = as.character(Sys.Date()),
    default_policy = "hybrid_per_species: leg_a where per-species fit reliable, else leg_b",
    components_present = names(d$categories_conus_sf)[names(d$categories_conus_sf) != "metadata"]
  )

  out_path <- if (DRY_RUN) sub("\\.json$", ".sf_preview.json", vf) else vf
  if (!DRY_RUN) file.copy(vf, paste0(vf, ".pre_sf_", stamp), overwrite = FALSE)
  write_json(d, out_path, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null")
  n_done <- n_done + 1L
  cat(sprintf("%s  %s  (leg_a=%d, leg_b=%d species)\n",
              if (DRY_RUN) "preview" else "LANDED",
              basename(out_path), sum(src == "leg_a"), sum(src == "leg_b")))
}

cat(sprintf("\nDone. %d variant file(s) %s.\n", n_done,
            if (DRY_RUN) "previewed" else "updated"))
quit(save = "no", status = 0)
