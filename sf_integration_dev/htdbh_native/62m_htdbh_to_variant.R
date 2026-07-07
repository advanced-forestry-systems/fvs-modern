#!/usr/bin/env Rscript
# =============================================================================
# 62m_htdbh_to_variant.R
#
# 62-series serializer for the CONUS HT-DBH (height-diameter) native hook.
# Maps the per-species CONUS HT-DBH coefficient table (from 62k best-AIC
# selection) onto each FVS variant using the SAME species crosswalk the other
# 62-series serializers use: the per-variant config JSON block
#   categories.species_definitions.FIAJSP  (FVS species index -> FIA SPCD),
# exactly the source 62_conus_to_variant_json.R reads
# (cfg$categories$species_definitions$FIAJSP) and the same FIAJSP mapping the
# engine loader GREGLOAD (src-converted/base/gregdghg.f90) uses at runtime.
#
# The engine HT-DBH hook (FVS_GREGHTDBH_COEF, GREGLOAD pattern) keys the
# coefficient table by FIA SPCD and resolves FVS species -> FIA at load time via
# FIAJSP, mirroring how greg_dg_coefficients.csv is emitted by 62i and consumed
# by gregdghg.f90. So the engine-facing deliverable is a single flat CSV keyed
# by FIA SPCD (header + SPCD model_id B1 B2 B3), same style as 62i. This script
# also emits per-variant coefficient tables and a per-variant coverage report so
# each variant`s HD species coverage is auditable and native-fallback species
# (variant species with no acceptable CONUS HD fit) are flagged.
#
# Usage:
#   Rscript 62m_htdbh_to_variant.R
#     [--coef <path to conus_htdbh_coefficients.csv>]
#     [--config_dir <path to config/calibrated>]
#     [--out_dir <path to config output dir>]
#
# Author: A. Weiskittel + Claude   Date: 2026-07-07
# =============================================================================
suppressWarnings(suppressMessages({ library(data.table); library(jsonlite) }))

args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d) { m <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (length(m)) sub(paste0("^--", n, "="), "", m[1]) else d }

WT        <- "/fs/scratch/PUOM0008/crsfaaron/wt-htdbh-ser"
COEF      <- ga("coef",       file.path(WT, "sf_integration_dev/htdbh_native/conus_htdbh_coefficients.csv"))
CONFIGDIR <- ga("config_dir", file.path(WT, "config/calibrated"))
OUTDIR    <- ga("out_dir",    file.path(WT, "config"))
PVDIR     <- file.path(OUTDIR, "htdbh_per_variant")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PVDIR,  showWarnings = FALSE, recursive = TRUE)

# The 25 CONUS variants, from ALL_VARIANTS in 62_conus_to_variant_json.R.
ALL_VARIANTS <- c("acd","ak","bm","ca","ci","cr","cs","ec","em","ie","kt",
                  "ls","nc","ne","oc","on","op","pn","sn","so","tt","ut",
                  "wc","ws","bc")

# ---- 1. Load the per-species CONUS HT-DBH coefficient table (62k output) ----
coef <- fread(COEF)
# Keep only rows with a usable model + B1,B2 (B3 optional; model 1 is 2-param).
coef <- coef[is.finite(SPCD) & is.finite(model_id) & is.finite(B1) & is.finite(B2)]
coef[, SPCD := as.integer(SPCD)]
setkey(coef, SPCD)
cov_spcd <- coef$SPCD

# ---- 2. Emit the single engine-facing flat CSV keyed by FIA SPCD -----------
# Same style/consumption as greg_dg_coefficients.csv (62i): header + rows the
# GREGLOAD-pattern loader reads, mapping FIA SPCD -> FVS species via FIAJSP.
engine_tbl <- coef[order(SPCD),
  .(SPCD, model_id = as.integer(model_id), B1, B2, B3 = fifelse(is.na(B3), 0, B3))]
ENGINE_CSV <- file.path(OUTDIR, "conus_htdbh_coefficients.csv")
fwrite(engine_tbl, ENGINE_CSV, quote = FALSE)

# ---- 3. Per-variant mapping via the reused FIAJSP crosswalk ----------------
read_fiajsp <- function(v) {
  f <- file.path(CONFIGDIR, paste0(v, ".json"))
  if (!file.exists(f)) return(NULL)
  cfg <- fromJSON(f, simplifyVector = TRUE)
  sd  <- cfg$categories$species_definitions
  if (is.null(sd) || is.null(sd$FIAJSP)) return(NULL)
  # FVS species index order = row order in species_definitions; keep JSP alpha.
  data.table(
    fvs_index = seq_along(sd$FIAJSP),
    JSP       = if (!is.null(sd$JSP)) trimws(sd$JSP) else NA_character_,
    FIAJSP    = trimws(sd$FIAJSP)
  )
}

cov_rows <- list()
for (v in ALL_VARIANTS) {
  xw <- read_fiajsp(v)
  if (is.null(xw)) {
    cov_rows[[v]] <- data.table(variant = v, n_species = NA_integer_,
      n_covered = NA_integer_, n_native_fallback = NA_integer_, pct_covered = NA_real_)
    next
  }
  # Real variant species = rows with a parseable FIA SPCD.
  xw[, spcd := suppressWarnings(as.integer(FIAJSP))]
  real <- xw[!is.na(spcd)]
  real[, covered := spcd %in% cov_spcd]

  # Per-variant coefficient table: one row per FVS species index, coefficients
  # attached where the CONUS HD fit exists; native-fallback flagged otherwise.
  pv <- merge(real[, .(fvs_index, JSP, SPCD = spcd)], coef,
              by = "SPCD", all.x = TRUE, sort = FALSE)
  setorder(pv, fvs_index)
  pv[, source := fifelse(is.na(model_id), "native_fallback", "conus_htdbh")]
  pv_out <- pv[, .(fvs_index, JSP, SPCD, model_id, B1, B2, B3, source)]
  fwrite(pv_out, file.path(PVDIR, paste0(v, "_htdbh_coefficients.csv")), quote = FALSE)

  n_species <- nrow(real)
  n_cov     <- sum(real$covered)
  cov_rows[[v]] <- data.table(
    variant = v, n_species = n_species, n_covered = n_cov,
    n_native_fallback = n_species - n_cov,
    pct_covered = round(100 * n_cov / n_species, 1))
}
cov <- rbindlist(cov_rows)
COV_CSV <- file.path(OUTDIR, "htdbh_variant_coverage.csv")
fwrite(cov, COV_CSV, quote = FALSE)

# ---- 4. Report -------------------------------------------------------------
cat("== 62m_htdbh_to_variant.R ==\n")
cat(sprintf("CONUS HD species with usable fit (engine table): %d -> %s\n",
            nrow(engine_tbl), ENGINE_CSV))
cat(sprintf("Per-variant tables: %s/{variant}_htdbh_coefficients.csv (%d variants)\n",
            PVDIR, sum(!is.na(cov$n_species))))
cat(sprintf("Coverage report: %s\n\n", COV_CSV))
cat("=== Per-variant HT-DBH coverage ===\n")
print(cov, row.names = FALSE)
cat(sprintf("\nTotals: %d variants; mean coverage %.1f%%\n",
            nrow(cov), mean(cov$pct_covered, na.rm = TRUE)))
