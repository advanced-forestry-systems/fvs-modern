# =============================================================================
# Title: Held-out species fit (DG_Kuehne B1 species-free)
# Author: A. Weiskittel
# Date: 2026-05-11
# Description: The species-free architecture claims that the eight-trait fixed
#              effects substitute for species random intercepts. The strongest
#              test of that claim is: fit B1 on N-k species, then predict for
#              the k held-out species using ONLY their trait vector W (no
#              species RE exists by construction, so this is a clean
#              extrapolation to new species).
#
# This script is a wrapper around 32_fit_dg_kuehne_speciesfree.R that
# (a) selects k held-out species, (b) writes their identity to a meta file
# alongside the fit, (c) calls the fit driver on the remaining species,
# (d) is read back by 80b_holdout_species_predict.R to score predictions
# on the held-out species.
#
# Selection criteria for held-out species:
#   - At least 3,000 observations in the training data (so we have enough
#     observed values to evaluate prediction error)
#   - Trait values present and non-NA across the 8 W columns
#   - Covers a range of softwood/hardwood, shade tolerance, max height
#
# CLI flags:
#   --n_holdout=K        number of species to hold out (default: 10)
#   --random_seed=N      seed for reproducible selection (default: 42)
#   --outdir=PATH        output dir (default: holdout_pilot_<timestamp>)
#   --submit_only        write the configuration but don't actually submit
#
# Run on Cardinal:
#   sbatch calibration/slurm/eval_holdout_species.sh
#     # which calls:
#     # Rscript --vanilla calibration/R/eval/80_holdout_species_fit.R --n_holdout=10
# =============================================================================

library(tidyverse)
library(data.table)

# --- CLI parsing -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}
has_flag <- function(name) any(grepl(paste0("^--", name, "$"), args))

N_HOLDOUT <- as.integer(get_arg("n_holdout", "10"))
SEED      <- as.integer(get_arg("random_seed", "42"))
OUT_DIR   <- get_arg("outdir",
                     paste0("calibration/output/conus/dg/holdout_pilot_",
                            format(Sys.time(), "%Y%m%d_%H%M%S")))
SUBMIT_ONLY <- has_flag("submit_only")

PROJ_ROOT <- "/users/PUOM0008/crsfaaron/fvs-modern"
OUT_DIR_FULL <- file.path(PROJ_ROOT, OUT_DIR)
dir.create(OUT_DIR_FULL, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Holdout species selection: N_HOLDOUT=%d SEED=%d\n", N_HOLDOUT, SEED))
cat(sprintf("Output dir: %s\n", OUT_DIR_FULL))

# --- Load data + identify candidate holdout species -------------------------

DATA_FILE   <- file.path(PROJ_ROOT,
  "calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds")
TRAITS_FILE <- file.path(PROJ_ROOT, "calibration/traits/species_traits.rds")

cat("Loading data ..."); flush.console()
dat <- as.data.table(readRDS(DATA_FILE))
traits <- as.data.table(readRDS(TRAITS_FILE))
cat(" done\n")

# Apply the same Kuehne filters as the production driver (kept inline so
# this script can be inspected without reading the driver).
MIN_OBS_SPECIES <- 5000
dat <- dat[
  is.finite(DBH1) & DBH1 >= 2.54 &
  is.finite(DBH2) &
  is.finite(CR1) & CR1 > 0 & CR1 <= 1.0 &
  is.finite(YEARS) & YEARS >= 1 & YEARS <= 20 &
  !is.na(EPA_L1_CODE) & !is.na(EPA_L2_CODE) & !is.na(EPA_L3_CODE) &
  EPA_L1_CODE != "" & EPA_L2_CODE != "" & EPA_L3_CODE != "" &
  TREESTATUS1 == 1 & TREESTATUS2 == 1 &
  is.finite(BAL_SW1) & BAL_SW1 >= 0 &
  is.finite(BAL_HW1) & BAL_HW1 >= 0 &
  is.finite(cspi) &
  is.finite(ba_x_rd) & is.finite(bal_x_rd)
]
sp_counts <- dat[, .N, by = SPCD][N >= MIN_OBS_SPECIES]
candidate_species <- sp_counts[N >= 3000][order(-N)]

cat(sprintf("Candidate species (>=3000 obs): %d\n", nrow(candidate_species)))

# Require traits non-NA across the 8 W columns
trait_cols <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                "vulnerability_score", "sensitivity")
traits_ok <- traits[, c("SPCD", trait_cols), with = FALSE][
  complete.cases(.SD), SPCD]
candidate_species <- candidate_species[SPCD %in% traits_ok]
cat(sprintf("Candidates with complete traits: %d\n", nrow(candidate_species)))

# Stratified random selection: include softwoods and hardwoods, range of
# shade tolerance, range of max height. Use trait-based stratification.
candidate_species <- merge(candidate_species,
                           traits[, c("SPCD", trait_cols, "GENUS"), with = FALSE],
                           by = "SPCD")

set.seed(SEED)
# Stratify by softwood (binary). Pick proportionally.
n_softwood <- round(N_HOLDOUT * mean(candidate_species$softwood))
n_hardwood <- N_HOLDOUT - n_softwood

hold_sw <- candidate_species[softwood == 1][sample(.N, min(n_softwood, .N))]
hold_hw <- candidate_species[softwood == 0][sample(.N, min(n_hardwood, .N))]
holdout <- rbind(hold_sw, hold_hw)

if (nrow(holdout) < N_HOLDOUT) {
  cat("WARNING: only", nrow(holdout), "species selected (requested",
      N_HOLDOUT, "). Add more candidates or reduce N_HOLDOUT.\n")
}

cat("\n=== Selected held-out species ===\n")
print(holdout[, .(SPCD, GENUS, softwood, shade_tolerance_num, max_ht_m, N)])

# --- Save selection ----------------------------------------------------------
holdout_meta_path <- file.path(OUT_DIR_FULL, "holdout_species.csv")
write_csv(holdout, holdout_meta_path)
cat("Saved holdout species selection to:", holdout_meta_path, "\n")

# --- Write a stripped data file (no held-out species) -----------------------
# To minimize compute, we save the SPCD vector to pass to the fit driver via
# a custom CLI arg. The fit driver needs to add a filter step.

holdout_spcd_path <- file.path(OUT_DIR_FULL, "holdout_spcd.txt")
writeLines(as.character(holdout$SPCD), holdout_spcd_path)
cat("Saved holdout SPCD list to:", holdout_spcd_path, "\n")

# --- Generate the submit script ---------------------------------------------
submit_path <- file.path(OUT_DIR_FULL, "submit_holdout_fit.sh")

submit_text <- sprintf('#!/bin/bash
#SBATCH --job-name=fvs_dg_kue_sf_holdout
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=%s/logs/fvs_dg_kue_sf_holdout_%%%%j.out
#SBATCH --error=%s/logs/fvs_dg_kue_sf_holdout_%%%%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

## Held-out species pilot for DG_Kuehne B1 species-free.
## Holds out: %s
## See holdout_species.csv for trait values of the held-out species.

set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd %s

mkdir -p %s

Rscript --vanilla calibration/R/32_fit_dg_kuehne_speciesfree.R \\
  --stan_file=calibration/stan/dg_kuehne2022_speciesfree.stan \\
  --outdir=%s \\
  --outname=dg_kuehne_b1_holdout \\
  --subsample=200000 \\
  --holdout_spcd_file=%s
',
  PROJ_ROOT, PROJ_ROOT,
  paste(holdout$SPCD, collapse = ","),
  PROJ_ROOT, OUT_DIR, OUT_DIR, holdout_spcd_path
)

writeLines(submit_text, submit_path)
Sys.chmod(submit_path, "0755")
cat("Generated submit script:", submit_path, "\n")

# --- Submit (or dry run) ----------------------------------------------------
if (SUBMIT_ONLY) {
  cat("\n--submit_only flag set; submit manually with:\n  sbatch", submit_path, "\n")
} else {
  cat("\nNot auto-submitting from this script. Review the submit script,\n")
  cat("ensure 32_fit_dg_kuehne_speciesfree.R supports --holdout_spcd_file,\n")
  cat("then submit manually:\n  sbatch", submit_path, "\n")
}

cat("\nNext step: after the fit lands, run\n")
cat("  Rscript calibration/R/eval/80b_holdout_species_predict.R --outdir=", OUT_DIR, "\n",
    sep = "")
cat("to score predictions on the held-out species using ONLY trait values.\n")

cat("\nDone.\n")
