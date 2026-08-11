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
# RELIABILITY GATE (added 2026-07-02)
#   Previously a species routed to leg_a whenever it merely had a per-species
#   intercept in categories_conus, with no quality threshold. The held-out
#   benchmark showed Leg A alone is usually the worst-informed arm, so routing
#   thinly-sampled species to leg_a imports noise. This version gates leg_a on
#   reliability:
#     --reliability presence   legacy behaviour (present => leg_a)
#     --reliability shrinkage  (default) leg_a only if the per-species posterior
#                              carries enough of its own information, measured by
#                              the partial-pooling shrinkage weight
#                                 w_j = tau^2 / (tau^2 + s_j^2)
#                              where s_j is the posterior sd of the species
#                              intercept (in categories_conus.species_intercepts)
#                              and tau is the between-species spread of the
#                              posterior means for that component within the
#                              variant (comp-level fallback when a variant has
#                              < 2 leg_a species). Route to leg_a iff w_j >= W_MIN.
#   w is unit-free and per-component, so it handles the very different intercept
#   scales across components (e.g. HT-DBH sd ~0.29 vs HCB sd ~0.05) automatically.
#   Default W_MIN = 0.5 (the species' own data supplies at least half the
#   information). The gate only changes hybrid_source_map; SF coefficients are
#   untouched, so a re-land against the same bundles produces a source-map-only
#   diff.
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
#     --reliability shrinkage --w_min 0.5 \
#     --dry_run TRUE
#
# Author: A. Weiskittel + Claude
# Date: 2026-05-21 (reliability gate 2026-07-02)
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
RELIABILITY  <- tolower(get_arg("reliability", "shrinkage"))  # presence | shrinkage
W_MIN        <- as.numeric(get_arg("w_min", "0.5"))

stopifnot(!is.null(COMPONENT), !is.null(BUNDLE_PREF))
if (!RELIABILITY %in% c("presence", "shrinkage"))
  stop("--reliability must be 'presence' or 'shrinkage'")

# component short name -> variant-json category key (matches Leg A keys)
KEY_MAP <- c(dg = "diameter_growth", hg = "height_growth",
             htdbh = "height_diameter", hcb = "height_crown_base",
             mort = "mortality", cr = "crown_recession",
             ingrowth = "ingrowth")
CKEY <- KEY_MAP[[COMPONENT]]
if (is.null(CKEY)) stop("Unknown component: ", COMPONENT)

cat("== 62b_speciesfree_to_variant_json.R ==\n")
cat("component:  ", COMPONENT, "-> category key:", CKEY, "\n")
cat("reliability:", RELIABILITY, if (RELIABILITY == "shrinkage") paste0("(w_min=", W_MIN, ")") else "", "\n")
cat("dry_run:    ", DRY_RUN, "\n\n")

# ---- Read the 61b bundle ---------------------------------------------------
bp <- function(suffix) file.path(BUNDLE_DIR, paste0(BUNDLE_PREF, suffix))
manifest <- fromJSON(bp("_sf_manifest.json"))
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

# ---- Reliability gate helper ----------------------------------------------
# Given the Leg A species_intercepts (SPCD, mean, sd), return the SPCD kept in
# leg_a plus a per-species reliability table. tau is estimated per component
# within the variant as the between-species sd of the posterior means; a
# comp-level tau (passed in) is used as a fallback when a variant has < 2
# species (posterior means alone cannot estimate a spread).
reliability_gate <- function(si_spcd, si_mean, si_sd, tau_fallback) {
  n <- length(si_spcd)
  if (n == 0) return(list(keep = integer(0), tab = NULL))
  if (RELIABILITY == "presence" || is.null(si_sd) || all(!is.finite(si_sd))) {
    return(list(keep = as.integer(si_spcd),
                tab = data.table(SPCD = as.integer(si_spcd),
                                 shrink_w = NA_real_, reliable = TRUE)))
  }
  tau <- if (n >= 2 && is.finite(sd(si_mean))) sd(si_mean) else tau_fallback
  if (!is.finite(tau) || tau <= 0) tau <- tau_fallback
  w <- tau^2 / (tau^2 + si_sd^2)              # partial-pooling shrinkage weight
  reliable <- is.finite(w) & w >= W_MIN
  list(keep = as.integer(si_spcd[reliable]),
       tab = data.table(SPCD = as.integer(si_spcd),
                        shrink_w = round(w, 4), reliable = reliable))
}

# ---- Per-variant landing ---------------------------------------------------
variant_files <- if (!is.null(VARIANT)) {
  file.path(CONFIG_DIR, paste0(VARIANT, ".json"))
} else {
  fs <- list.files(CONFIG_DIR, pattern = "\\.json$", full.names = TRUE)
  fs[!grepl("(_draws|pre_conus|pre_sf|sf_preview)", fs)]
}

# comp-level tau fallback: pooled between-species spread of Leg A means for this
# component across all variant files (computed once, before the per-variant loop).
comp_means <- c()
for (vf0 in variant_files) {
  if (!file.exists(vf0)) next
  d0 <- tryCatch(fromJSON(vf0), error = function(e) NULL); if (is.null(d0)) next
  si0 <- d0$categories_conus[[CKEY]]$species_intercepts
  if (!is.null(si0) && !is.null(si0$mean)) comp_means <- c(comp_means, as.numeric(si0$mean))
}
TAU_FALLBACK <- if (length(comp_means) >= 2 && is.finite(sd(comp_means))) sd(comp_means) else 1.0

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
n_done <- 0L
tot_present <- 0L; tot_legA <- 0L; tot_legB <- 0L; tot_gated <- 0L
for (vf in variant_files) {
  if (!file.exists(vf)) { cat("skip (missing):", vf, "\n"); next }
  d <- fromJSON(vf, simplifyVector = TRUE)

  # Leg A per-species entries for this component (SPCD, mean, sd).
  si_spcd <- integer(0); si_mean <- numeric(0); si_sd <- numeric(0)
  if (!is.null(d$categories_conus) && !is.null(d$categories_conus[[CKEY]]) &&
      !is.null(d$categories_conus[[CKEY]]$species_intercepts)) {
    si <- d$categories_conus[[CKEY]]$species_intercepts
    si_spcd <- as.integer(si$SPCD)
    si_mean <- if (!is.null(si$mean)) as.numeric(si$mean) else rep(NA_real_, length(si_spcd))
    si_sd   <- if (!is.null(si$sd))   as.numeric(si$sd)   else rep(NA_real_, length(si_spcd))
  }
  present_spcd <- si_spcd

  gate <- reliability_gate(si_spcd, si_mean, si_sd, TAU_FALLBACK)
  legA_spcd <- gate$keep

  # Hybrid map over the union of species either leg knows.
  all_spcd <- sort(unique(c(bundle_spcd, present_spcd)))
  src <- ifelse(all_spcd %in% legA_spcd, "leg_a", "leg_b")

  n_present <- length(present_spcd)
  n_legA <- sum(src == "leg_a"); n_legB <- sum(src == "leg_b")
  n_gated <- length(setdiff(present_spcd, legA_spcd))  # present but demoted to leg_b
  tot_present <- tot_present + n_present; tot_legA <- tot_legA + n_legA
  tot_legB <- tot_legB + n_legB; tot_gated <- tot_gated + n_gated

  sf_block <- sf_common
  # restrict the species trait table to species this variant / either leg uses
  sp_keep <- species[SPCD %in% all_spcd]
  sf_block$species <- as.list(sp_keep)
  sf_block$hybrid_source_map <- list(SPCD = all_spcd, source = src)
  # audit trail: reliability of each Leg A species and the gate applied
  sf_block$hybrid_reliability <- list(
    gate = RELIABILITY, w_min = if (RELIABILITY == "shrinkage") W_MIN else NA,
    n_present = n_present, n_leg_a = n_legA, n_gated_to_leg_b = n_gated,
    SPCD = if (!is.null(gate$tab)) gate$tab$SPCD else integer(0),
    shrink_w = if (!is.null(gate$tab)) gate$tab$shrink_w else numeric(0))

  if (is.null(d$categories_conus_sf)) d$categories_conus_sf <- list()
  d$categories_conus_sf[[CKEY]] <- sf_block
  d$categories_conus_sf$metadata <- list(
    pipeline_version = "fvs-conus species-free (Leg B)",
    integration_date = as.character(Sys.Date()),
    default_policy = paste0("hybrid_per_species: leg_a where per-species fit reliable (",
                            RELIABILITY,
                            if (RELIABILITY == "shrinkage") paste0(", w>=", W_MIN) else "",
                            "), else leg_b"),
    components_present = names(d$categories_conus_sf)[names(d$categories_conus_sf) != "metadata"]
  )

  out_path <- if (DRY_RUN) sub("\\.json$", ".sf_preview.json", vf) else vf
  if (!DRY_RUN) file.copy(vf, paste0(vf, ".pre_sf_", stamp), overwrite = FALSE)
  write_json(d, out_path, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null")
  n_done <- n_done + 1L
  cat(sprintf("%s  %-18s present=%3d  leg_a=%3d  leg_b=%3d  (gated present->leg_b: %d)\n",
              if (DRY_RUN) "preview" else "LANDED",
              basename(out_path), n_present, n_legA, n_legB, n_gated))
}

cat(sprintf("\nDone. %d variant file(s) %s. component=%s reliability=%s%s\n",
            n_done, if (DRY_RUN) "previewed" else "updated", COMPONENT, RELIABILITY,
            if (RELIABILITY == "shrinkage") paste0(" w_min=", W_MIN) else ""))
cat(sprintf("TOTALS  present(leg_a-eligible)=%d  ->  leg_a=%d  leg_b=%d  (demoted by gate: %d)\n",
            tot_present, tot_legA, tot_legB, tot_gated))
quit(save = "no", status = 0)
