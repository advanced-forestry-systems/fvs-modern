#!/usr/bin/env Rscript
# =============================================================================
# 62g_greg_to_variant_json.R
#
# Greg Johnson (fvs_remodeling) companion to 62b_speciesfree_to_variant_json.R.
# Builds the categories_conus_greg block that gives the native FVS engine the
# THIRD user-selectable CONUS option (version=conus_greg), alongside
# categories_conus (species-dependent) and categories_conus_sf (species-free).
#
# Greg's parameters are PER-SPECIES direct fits (not trait + ecoregion RE), so
# the block is simpler than the Leg A / Leg B blocks: three component tables
# keyed by SPCD, plus metadata. Greg's equations are CONUS-global (not
# per-variant), so the same block is emitted; the runtime filters to a variant's
# species set and applies the documented ~9 to 11 percent species fallback
# (fvs-conus trait / softwood-hardwood median), exactly as
# conus_eq_projector_greg.R does.
#
# Greg's components:
#   diameter_growth  dg_parms.RDS   (84 spp; B0..B7; predictors dbh, cr, ht, bal, elev, EMT)
#   height_growth    hg_parms.RDS   (96 spp; B0..B8; max_ht, ccfl, cch, elev, TD, EMT)
#   survival         mort_parm_base_rate_cr_cch.RDS (per spp; b0..b4; cch-gompit)
# crown: NONE in Greg's repo -> uses the fvs-conus CR recession kernel (disclosed).
#
# Climate dependency: Greg's kernels need EMT and TD = MWMT - MCMT from ClimateNA
# 1991-2020 normals at the stand; the engine path must supply these per stand
# (the projector uses greg_emt_td_lookup.rds). Recorded in the block metadata.
#
# SAFETY: --dry_run=TRUE (default) writes greg_block_preview.json (the standalone
# block) and does NOT touch production configs. With --variant and
# --dry_run=FALSE it merges categories_conus_greg into config/calibrated/{variant}.json
# after a .pre_greg_<timestamp> backup.
#
# Usage:
#   Rscript 62g_greg_to_variant_json.R --greg_rds=/users/PUOM0008/crsfaaron/fvs_remodeling/rds \
#       --out=greg_block_preview.json [--variant=ne --config_dir=config/calibrated --dry_run=FALSE]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d = NULL) { m <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(m)) return(d); sub(paste0("^--", n, "="), "", m[1]) }
GREG_RDS  <- ga("greg_rds", "/users/PUOM0008/crsfaaron/fvs_remodeling/rds")
OUT       <- ga("out", "greg_block_preview.json")
VARIANT   <- ga("variant", NA_character_)
CONFIG_DIR<- ga("config_dir", "config/calibrated")
DRY_RUN   <- toupper(ga("dry_run", "TRUE")) != "FALSE"
REPO_TAG  <- ga("repo_tag", "TBD-confirm-with-author")  # github.com/gregjohnsonbiometrics/fvs_remodeling

read_rds_dt <- function(f) { p <- file.path(GREG_RDS, f); if (!file.exists(p)) stop("missing: ", p)
  as.data.table(readRDS(p)) }

dg <- read_rds_dt("dg_parms.RDS")
hg <- read_rds_dt("hg_parms.RDS")
mo <- read_rds_dt("mort_parm_base_rate_cr_cch.RDS")

# normalize the species-code column name (dg/hg use 'spcd', mort uses 'SPCD')
norm_spcd <- function(d) { sc <- intersect(c("spcd","SPCD","Spcd"), names(d))[1]
  if (is.na(sc)) stop("no species code col in: ", paste(names(d), collapse=",")); setnames(d, sc, "SPCD"); d }
dg <- norm_spcd(dg); hg <- norm_spcd(hg); mo <- norm_spcd(mo)

# keep only converged fits where the flag exists
if ("isConv" %in% names(dg)) dg <- dg[isConv == TRUE | is.na(isConv)]
if ("isConv" %in% names(hg)) hg <- hg[isConv == TRUE | is.na(isConv)]

coef_cols <- function(d) grep("^B[0-9]+$|^b[0-9]+$", names(d), value = TRUE)
sp_table  <- function(d) {
  cc <- coef_cols(d)
  rows <- lapply(seq_len(nrow(d)), function(i) {
    r <- as.list(d[i, c("SPCD", cc), with = FALSE]); r$SPCD <- as.integer(r$SPCD); r })
  rows
}

block <- list(
  parameter_source = "greg_fvs_remodeling",
  repo = list(url = "https://github.com/gregjohnsonbiometrics/fvs_remodeling",
              tag = REPO_TAG,
              note = "Confirm release tag/commit with Greg Johnson before production"),
  climate_dependency = list(
    required = c("EMT", "TD"),
    definition = "EMT and TD = MWMT - MCMT from ClimateNA 1991-2020 normals at the stand; engine must supply per stand",
    lookup = "greg_emt_td_lookup.rds"),
  crown = list(source = "fvs_conus_cr_recession",
               greg_alternative = list(
                 coefficients = "greg_crown_coefficients.csv",
                 form = "dHTLC = (CHT-CHTLC)*(1 - exp(B0 + B1*dHT_ft + B2*dCCH))  (dHT in FEET)",
                 units = "imperial (confirmed by A/B: feet bias +0.13 m vs meters -2.76 m)"),
               note = "A/B on 102k held-out FIA pairs (2026-07-03): fvs-conus CR-recession kernel wins on RMSE (5.61 vs 5.98 m) and is unbiased; Greg's crown (fvs_remodeling 8310b2a) is nearly unbiased with slightly higher r but higher RMSE, and its dCCH term is inert under our CCH scale. Kept kernel as default; Greg retained as a documented alternative pending CCH-scale reconciliation with Greg."),
  species_fallback = list(
    rule = "species outside Greg's fitted sets use the fvs-conus trait / softwood-hardwood median fallback",
    expected_fraction = "~0.09 to 0.11 of trees, logged per stand"),
  components = list(
    diameter_growth = list(form = "greg_est_dg", predictors = c("dbh","cr","ht","bal","elev","EMT"),
                           n_species = nrow(dg), coef_names = coef_cols(dg), species_params = sp_table(dg),
                           # keyword-selectable site driver (DGDRIVER). A/B 2026-07-03: driver is a
                           # minor DG lever (~0.6% RMSE); CSPI marginally best. EMT variant is Greg's
                           # native (imperial, needs ClimateNA); {none,elev,bgi,cspi,esi} refit on CONUS pairs.
                           site_driver = "cspi",
                           site_driver_options = c("none","elev","bgi","cspi","esi","emt"),
                           coefficients_by_driver = list(
                             none = "greg_dg_coefficients_none.csv", elev = "greg_dg_coefficients_elev.csv",
                             bgi  = "greg_dg_coefficients_bgi.csv",  cspi = "greg_dg_coefficients_cspi.csv",
                             esi  = "greg_dg_coefficients_esi.csv",  emt  = "greg_dg_coefficients.csv")),
    height_growth   = list(form = "greg_est_hg", predictors = c("max_ht","ccfl","cch","elev","TD","EMT"),
                           n_species = nrow(hg), coef_names = coef_cols(hg), species_params = sp_table(hg),
                           site_driver = "emt",
                           note = "HG site-driver family deferred: est_hg needs TD/EMT (ClimateNA) + a per-species max-height asymptote; refit in-engine or with Greg, not offline."),
    survival        = list(form = "greg_gompit",
                           equation = "P_surv = 1 - exp(-exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4))  (annual)",
                           predictors = c("cr","cch"),
                           n_species = nrow(mo), coef_names = coef_cols(mo), species_params = sp_table(mo),
                           # keyword-selectable site driver (MORTDRIVER). A/B 2026-07-03: driver MATTERS
                           # for survival (~5.7% log-loss); BGI best. Default BGI.
                           site_driver = "bgi",
                           site_driver_options = c("none","elev","bgi","cspi","esi"),
                           coefficients_by_driver = list(
                             none = "greg_mort_coefficients_none.csv", elev = "greg_mort_coefficients_elev.csv",
                             bgi  = "greg_mort_coefficients_bgi.csv",  cspi = "greg_mort_coefficients_cspi.csv",
                             esi  = "greg_mort_coefficients_esi.csv"))
  ),
  stand_level = list(
    garcia_constraint = FALSE,
    note = "density is emergent from the cch-gompit mortality; no fvs-conus Garcia self-thinning / BA carrying-capacity overlay. Optional softened SDIMAX cap only."),
  components_present = c("diameter_growth", "height_growth", "survival", "crown_change")
)

cat(sprintf("Greg block built: DG %d spp, HG %d spp, survival %d spp.\n",
            nrow(dg), nrow(hg), nrow(mo)))

if (DRY_RUN || is.na(VARIANT)) {
  writeLines(toJSON(list(categories_conus_greg = block), auto_unbox = TRUE, pretty = TRUE, digits = 10), OUT)
  cat("DRY RUN: wrote standalone block ->", OUT, "\n")
  cat("To land into a variant config: --variant=ne --dry_run=FALSE\n")
} else {
  cfgf <- file.path(CONFIG_DIR, paste0(tolower(VARIANT), ".json"))
  if (!file.exists(cfgf)) stop("config not found: ", cfgf)
  bak <- paste0(cfgf, ".pre_greg_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(cfgf, bak); cat("backup ->", bak, "\n")
  cfg <- fromJSON(cfgf, simplifyVector = FALSE)
  cfg$categories_conus_greg <- block
  writeLines(toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = 10), cfgf)
  cat("LANDED categories_conus_greg into", cfgf, "\n")
}
