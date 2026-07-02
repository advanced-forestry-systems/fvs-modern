#!/usr/bin/env Rscript
# =============================================================================
# 61b_extract_speciesfree_summaries.R
#
# Species-free (Leg B) companion to 61_extract_conus_summaries.R. The original
# 61 understands the per-species "cspi_traits1" shape. This script extracts the
# trait-and-hierarchy (B1) shape into a standardized bundle that
# 62b_speciesfree_to_variant_json.R lands into fvs-modern as a
# categories_conus_sf block.
#
# A species-free component's linear predictor is, generically:
#   eta = <fixed/covariate terms> + (W_species . gamma)
#       + z_L1[L1] + z_L2[L2] + z_L3[L3] + z_FT[FT]
# where W_species are the standardized species traits. To evaluate this at
# runtime in the engine for ANY species (the whole point of species-free), the
# engine needs: the global fixed effects, the trait coefficient vector gamma,
# the per-trait standardization constants (so raw traits can be standardized),
# and the random-effect tables keyed by the real EPA / forest-type codes.
#
# This reads a fit's meta (level keys, trait_cols) and the full fit object
# (for the RE realizations, which are not in the summary CSV), recomputes the
# trait standardization exactly as the driver did, and writes:
#   {out}_sf_fixed.csv     global fixed + covariate coefficients + sigmas
#   {out}_sf_gamma.csv      trait_col, gamma mean/sd, standardization mean/sd
#   {out}_sf_re_L1.csv      EPA L1 code, RE mean/sd   (and L2, L3, FT)
#   {out}_sf_species.csv     SPCD, raw + standardized traits, trait_effect mean
#   {out}_sf_manifest.json   component, form, counts, file list, date
#
# Usage (from fvs-modern so calibration/ resolves):
#   Rscript calibration/R/61b_extract_speciesfree_summaries.R \
#     --component hg \
#     --fit calibration/output/conus/hg/v5/hg_sf_v5_100k_prod_fit.rds \
#     --meta calibration/output/conus/hg/v5/hg_sf_v5_100k_prod_meta.rds \
#     --traits calibration/traits/species_traits_v2.rds \
#     --outdir calibration/output/conus/sf_integration
#
# Author: A. Weiskittel + Claude
# Date: 2026-05-21
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cmdstanr)
  library(posterior)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) > 0) return(sub(paste0("^--", name, "="), "", m[1]))
  # also allow space-separated form "--name value"
  i <- which(args == paste0("--", name))
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}

COMPONENT <- get_arg("component", "unknown")
FIT_PATH  <- get_arg("fit")
META_PATH <- get_arg("meta")
TRAITS    <- get_arg("traits", "calibration/traits/species_traits_v2.rds")
OUTDIR    <- get_arg("outdir", "calibration/output/conus/sf_integration")
OUTNAME   <- get_arg("outname", COMPONENT)

stopifnot(!is.null(FIT_PATH), !is.null(META_PATH))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("== 61b_extract_speciesfree_summaries.R ==\n")
cat("component:", COMPONENT, "\n")
cat("fit:      ", FIT_PATH, "\n")
cat("meta:     ", META_PATH, "\n\n")

meta <- readRDS(META_PATH)
trait_cols <- meta$trait_cols
sp_levels  <- meta$sp_levels
L1_levels  <- meta$L1_levels
L2_levels  <- meta$L2_levels
L3_levels  <- meta$L3_levels
FT_levels  <- meta$FT_levels

cat("Reading fit object (may be several GB) ...\n"); flush.console()
fit <- readRDS(FIT_PATH)

# ---- Identify which parameter blocks exist (component-agnostic) -------------
sv <- tryCatch(fit$metadata()$stan_variables, error = function(e) character(0))
# Drop generated quantities and per-observation vectors (they make $summary slow
# and are not part of the coefficient bundle): log_lik, mu_*, eta*, p_*, y_rep,
# anything with a _gq suffix, lp__, and the *_raw non-centered helpers.
# Note: richer components (HG v5 BGI, DG v8) carry extra structure such as
# gamma_site, species_site_slope, z_L1_bgi. The generic bundle keeps the
# standard B1 pieces (scalars, gamma, z_L1/L2/L3/FT); those extras are ignored
# here and need component-specific handling in 62b.
exclude <- c("log_lik", "mu_pred", "mu_raw", "mu_a", "p_mort_annual",
             "p_surv_annual", "y_rep", "eta", "eta_safe", "eta_gq",
             "trait_effect", "lp__")
keep <- sv[!sv %in% exclude & !grepl("_raw$|_gq$", sv)]
re_names    <- intersect(c("z_L1", "z_L2", "z_L3", "z_FT"), keep)
gamma_here  <- "gamma" %in% keep
fixed_names <- setdiff(keep, c(re_names, "gamma"))

cat("stan variables kept:", paste(keep, collapse = ", "), "\n")
cat("RE tables present:   ", paste(re_names, collapse = ", "), "\n\n")

summ <- fit$summary(variables = keep, "mean", "sd",
                    ~quantile(.x, probs = c(0.05, 0.95)))
setDT(summ)
setnames(summ, c("5%", "95%"), c("q5", "q95"), skip_absent = TRUE)

idx_of <- function(v) as.integer(sub(".*\\[(\\d+)\\]$", "\\1", v))

# ---- 1. Fixed / covariate coefficients (scalars) ---------------------------
fixed <- summ[variable %in% fixed_names]
fwrite(fixed, file.path(OUTDIR, paste0(OUTNAME, "_sf_fixed.csv")))
cat("wrote fixed (", nrow(fixed), "params)\n")

# ---- 2. Trait standardization (replicate the driver exactly) ---------------
# Driver: traits[match(sp_levels, SPCD)] -> for each col, NA imputed to median,
# then (x - mean) / sd computed on the imputed column.
traits <- as.data.table(readRDS(TRAITS))
tsub <- traits[match(sp_levels, SPCD), c("SPCD", trait_cols), with = FALSE]
Wraw <- as.matrix(tsub[, trait_cols, with = FALSE])
scale_mean <- numeric(length(trait_cols))
scale_sd   <- numeric(length(trait_cols))
Wstd <- Wraw
for (j in seq_along(trait_cols)) {
  col <- Wraw[, j]
  na  <- is.na(col)
  if (any(na)) col[na] <- median(col[!na], na.rm = TRUE)
  m <- mean(col); s <- sd(col)
  scale_mean[j] <- m; scale_sd[j] <- s
  Wstd[, j] <- (col - m) / s
}

# ---- 3. gamma aligned to trait_cols + standardization constants ------------
if (gamma_here) {
  g <- summ[grepl("^gamma\\[", variable)]
  g[, j := idx_of(variable)]
  setorder(g, j)
  gamma_tbl <- data.table(
    trait_col  = trait_cols[g$j],
    gamma_mean = g$mean, gamma_sd = g$sd, gamma_q5 = g$q5, gamma_q95 = g$q95,
    scale_mean = scale_mean[g$j], scale_sd = scale_sd[g$j]
  )
  fwrite(gamma_tbl, file.path(OUTDIR, paste0(OUTNAME, "_sf_gamma.csv")))
  cat("wrote gamma (", nrow(gamma_tbl), "traits)\n")
  gamma_mean_vec <- g$mean[order(g$j)]
} else {
  gamma_mean_vec <- rep(0, length(trait_cols))
}

# ---- 4. Random-effect tables keyed by real codes ---------------------------
write_re <- function(zname, levels_vec, suffix) {
  if (!(zname %in% re_names)) return(invisible(NULL))
  z <- summ[grepl(paste0("^", zname, "\\["), variable)]
  z[, i := idx_of(variable)]
  setorder(z, i)
  out <- data.table(level = levels_vec[z$i], mean = z$mean, sd = z$sd,
                    q5 = z$q5, q95 = z$q95)
  fwrite(out, file.path(OUTDIR, paste0(OUTNAME, "_sf_re_", suffix, ".csv")))
  cat("wrote RE", suffix, "(", nrow(out), "levels)\n")
}
write_re("z_L1", L1_levels, "L1")
write_re("z_L2", L2_levels, "L2")
write_re("z_L3", L3_levels, "L3")
write_re("z_FT", FT_levels, "FT")

# ---- 5. Per-species trait bundle + trait_effect ----------------------------
trait_effect <- as.numeric(Wstd %*% gamma_mean_vec)
species_tbl <- data.table(SPCD = sp_levels)
for (j in seq_along(trait_cols)) {
  species_tbl[[paste0("raw_", trait_cols[j])]] <- Wraw[, j]
  species_tbl[[paste0("std_", trait_cols[j])]] <- Wstd[, j]
}
species_tbl[, trait_effect_mean := trait_effect]
fwrite(species_tbl, file.path(OUTDIR, paste0(OUTNAME, "_sf_species.csv")))
cat("wrote species (", nrow(species_tbl), "species; trait_effect range",
    round(min(trait_effect), 3), "to", round(max(trait_effect), 3), ")\n")

# ---- 6. Manifest -----------------------------------------------------------
manifest <- list(
  component = COMPONENT,
  form = meta$form,
  stan_file = meta$stan_file,
  n_obs = meta$n_obs,
  extracted = as.character(Sys.Date()),
  trait_cols = trait_cols,
  counts = list(species = length(sp_levels), L1 = length(L1_levels),
                L2 = length(L2_levels), L3 = length(L3_levels),
                FT = length(FT_levels)),
  fixed_params = fixed$variable,
  re_tables = re_names,
  files = list(
    fixed   = paste0(OUTNAME, "_sf_fixed.csv"),
    gamma   = paste0(OUTNAME, "_sf_gamma.csv"),
    species = paste0(OUTNAME, "_sf_species.csv"),
    re      = paste0(OUTNAME, "_sf_re_", sub("z_", "", re_names), ".csv")
  )
)
write_json(manifest, file.path(OUTDIR, paste0(OUTNAME, "_sf_manifest.json")),
           auto_unbox = TRUE, pretty = TRUE)
cat("\nwrote manifest. Done. Output dir:", OUTDIR, "\n")

# Explicit clean exit so SLURM afterok dependencies are not tripped by any
# trailing parser noise once all writes have completed.
quit(save = "no", status = 0)
