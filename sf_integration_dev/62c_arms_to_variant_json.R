#!/usr/bin/env Rscript
# =============================================================================
# 62c_arms_to_variant_json.R
#
# Companion to 62b. Lands the two NEW blocks that complete the three-arm CONUS
# design into the per-variant config JSONs:
#
#   --block=mod       categories_conus_mod.{component}      (Bayesian modifier
#                     layer: management / disturbance / driver multipliers, from
#                     70_fit_modifiers.R). Shared across ALL three arms.
#   --block=organon   categories_conus_organon.{component}  (arm 1, ORGANON-form
#                     coefficients). Selected by version=conus_organon.
#
# Arms 2 (conus / Leg A) and 3 (conus_sf / Leg B) already exist in the configs;
# this script adds arm 1 and the shared modifier layer, so all three arms are
# selectable via config_loader version and carry the same modifier layer.
#
# SAFETY: --dry_run=TRUE (default) writes {variant}.armpreview.json instead of
# the real config. --dry_run=FALSE makes a .pre_arm_<ts> backup then edits.
#
# Usage:
#   Rscript 62c_arms_to_variant_json.R --block=mod --component=dg \
#     --bundle=<dir>/dg_modifier_bundle.json --config_dir=config/calibrated \
#     --dry_run=TRUE
#   Rscript 62c_arms_to_variant_json.R --block=organon --component=dg \
#     --bundle=<dir>/dg_organon_bundle.json --config_dir=config/calibrated
#
# Author: A. Weiskittel + Claude (OODA autopilot)  Date: 2026-07-02
# =============================================================================
suppressPackageStartupMessages({ library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d = NULL) { m <- grep(paste0("^--", n, "="), args, value = TRUE); if (length(m)) sub(paste0("^--", n, "="), "", m[1]) else d }
BLOCK <- ga("block")                       # mod | organon | stand
COMPONENT <- ga("component", "dg")
STAND_KEY <- ga("stand_key", "survival")   # for --block=stand: survival | bagrowth
BUNDLE <- ga("bundle")
CONFIG_DIR <- ga("config_dir", "config/calibrated")
VARIANT <- ga("variant", NULL)
DRY_RUN <- toupper(ga("dry_run", "TRUE")) %in% c("TRUE","T","1","YES")
stopifnot(!is.null(BLOCK), BLOCK %in% c("mod","organon","stand"), !is.null(BUNDLE), file.exists(BUNDLE))

# --block=stand lands categories_conus_stand.{survival|bagrowth}; the sub-key is
# STAND_KEY (not a tree component). mod/organon land categories_conus_*.{component}.
if (BLOCK == "stand") {
  stopifnot(STAND_KEY %in% c("survival","bagrowth","topht","stems"))
  TOPKEY <- "categories_conus_stand"; CKEY <- STAND_KEY
  PREVIEW_SUFFIX <- ".standpreview.json"; BACKUP_TAG <- "pre_stand"
} else {
  KEY_MAP <- c(dg="diameter_growth", hg="height_growth", htdbh="height_diameter",
               hcb="height_crown_base", mort="mortality", cr="crown_recession")
  CKEY <- KEY_MAP[[COMPONENT]]; if (is.null(CKEY)) stop("unknown component ", COMPONENT)
  TOPKEY <- if (BLOCK == "mod") "categories_conus_mod" else "categories_conus_organon"
  PREVIEW_SUFFIX <- ".armpreview.json"; BACKUP_TAG <- "pre_arm"
}
bundle <- fromJSON(BUNDLE, simplifyVector = TRUE)

cat("== 62c_arms_to_variant_json.R ==\n block:", BLOCK, "-> ", TOPKEY, "\n sub-key:",
    if (BLOCK == "stand") STAND_KEY else COMPONENT, "->", CKEY, "\n dry_run:", DRY_RUN, "\n\n")

variant_files <- if (!is.null(VARIANT)) file.path(CONFIG_DIR, paste0(VARIANT, ".json")) else {
  fs <- list.files(CONFIG_DIR, pattern = "\\.json$", full.names = TRUE)
  fs[!grepl("(_draws|pre_conus|pre_sf|pre_arm|pre_stand|sf_preview|armpreview|standpreview)", fs)]
}
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S"); n_done <- 0L
for (vf in variant_files) {
  if (!file.exists(vf)) { cat("skip (missing):", vf, "\n"); next }
  d <- fromJSON(vf, simplifyVector = TRUE)
  if (is.null(d[[TOPKEY]])) d[[TOPKEY]] <- list()
  d[[TOPKEY]][[CKEY]] <- bundle
  d[[TOPKEY]]$metadata <- list(
    pipeline_version = switch(BLOCK,
      mod = "conus Bayesian modifier layer",
      organon = "conus ORGANON-form arm (arm 1)",
      stand = "conus stand-level constraint layer (survival / BA-growth)"),
    integration_date = as.character(Sys.Date()),
    components_present = names(d[[TOPKEY]])[names(d[[TOPKEY]]) != "metadata"],
    note = switch(BLOCK,
      mod = "Shared multiplicative modifier applied by any arm: growth = base * exp(mod_eta).",
      organon = "Arm 1 ORGANON-form coefficients; select with version=conus_organon.",
      stand = paste("Stand-level target that reconciles the summed tree predictions",
                    "(disaggregation): survival -> M_stand for the kappa hazard solve;",
                    "bagrowth -> stand BA increment for the tree-DG scale;",
                    "topht -> García/GADA H2|H1 target for the tree height-growth scale;",
                    "stems -> García N(t) surviving-stems target (reuses the kappa hazard",
                    "solve). Applied by any arm.")))
  out_path <- if (DRY_RUN) sub("\\.json$", PREVIEW_SUFFIX, vf) else vf
  if (!DRY_RUN) file.copy(vf, paste0(vf, ".", BACKUP_TAG, "_", stamp), overwrite = FALSE)
  write_json(d, out_path, auto_unbox = TRUE, pretty = TRUE, digits = 10, null = "null")
  n_done <- n_done + 1L
  cat(sprintf("%s  %s  [%s.%s]\n", if (DRY_RUN) "preview" else "LANDED", basename(out_path), TOPKEY, CKEY))
}
cat(sprintf("\nDone. %d variant file(s) %s.\n", n_done, if (DRY_RUN) "previewed" else "updated"))
quit(save = "no", status = 0)
