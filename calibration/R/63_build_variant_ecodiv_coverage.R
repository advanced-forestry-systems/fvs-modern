#!/usr/bin/env Rscript
# =============================================================================
# 63_build_variant_ecodiv_coverage.R
#
# Build the variant x ecodivision coverage table consumed by
# 62_conus_to_variant_json.R. The output is a one-time CSV listing, for each
# FVS variant, the Bailey ecodivisions that overlap the variant and the
# weight to apply when assembling the variant's coverage-weighted ecodivision
# random effect from the CONUS posterior.
#
# Input
# -----
#   calibration/data/fvs_variant_ecodivision_map.csv
#     Pre-computed area-weighted overlay of FVS variant polygons and Bailey
#     ecodivision polygons. Columns: FVSVariant, DIVISION, area, area_pct.
#     Covers 19 of 25 variants. Built upstream by the variant boundary
#     spatial pipeline.
#
# Output
# ------
#   calibration/data/variant_ecodiv_coverage.csv
#     Columns:
#       variant     lowercase FVS variant code (acd, ne, pn, ...)
#       ecodiv      Bailey ecodivision code (210, M210, ...)
#       area_m2     overlap area in square meters (NA for proxied variants)
#       area_pct    percent of variant area in that ecodivision (proxies copy from source)
#       weight      normalized weight, 0..1, sums to 1.0 per variant
#       source      "area" for direct overlay; "proxy_<src>" for proxied variants
#
# Coverage handling
# -----------------
# 19 variants have direct area-weighted entries from the spatial overlay.
# 6 variants (acd, bc, kt, oc, on, op) lack variant boundary geometry. For
# these we proxy from a related variant whose state coverage best matches:
#
#   acd -> ne   (ME, NH, VT are subset of NE)
#   bc  -> pn   (Pacific Northwest US proxy for British Columbia)
#   kt  -> em   (Montana proxy for Klamath-Tahoe)
#   oc  -> ca   (California Coast proxy for Oregon-California Other)
#   on  -> ls   (Great Lakes US proxy for Ontario)
#   op  -> pn   (Oregon Pacific proxy)
#
# Water bodies (DIVISION starting with "Wate") are filtered out.
#
# Usage
# -----
#   Rscript calibration/R/63_build_variant_ecodiv_coverage.R
#   Rscript calibration/R/63_build_variant_ecodiv_coverage.R --dry-run
#
# Author: A. Weiskittel
# Date:   2026-05-10
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(logger)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

project_root    <- Sys.getenv("FVS_PROJECT_ROOT", normalizePath(file.path(getwd(), "..")))
calibration_dir <- file.path(project_root, "calibration")
data_dir        <- file.path(calibration_dir, "data")

INPUT_MAP   <- file.path(data_dir, "fvs_variant_ecodivision_map.csv")
OUTPUT_CSV  <- file.path(data_dir, "variant_ecodiv_coverage.csv")

# All 25 variants expected by 62_conus_to_variant_json.R
ALL_VARIANTS <- c("acd","ak","bm","bc","ca","ci","cr","cs","ec","em","ie","kt",
                  "ls","nc","ne","oc","on","op","pn","sn","so","tt","ut","wc","ws")

# Proxy map for variants without direct geometry overlay. Values are
# lowercase source variant codes already present in the area map.
PROXY_SOURCES <- list(
  acd = "ne",
  bc  = "pn",
  kt  = "em",
  oc  = "ca",
  on  = "ls",
  op  = "pn"
)

WATER_PATTERN <- "^Wate"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  out <- list(dry_run = FALSE)
  for (a in args) {
    if (a == "--dry-run") out$dry_run <- TRUE
  }
  out
}

# ---------------------------------------------------------------------------
# Core builders
# ---------------------------------------------------------------------------

build_direct_coverage <- function(map_df) {
  ## Standardize column names and normalize weights per variant.
  map_df %>%
    rename(variant_upper = FVSVariant, ecodiv = DIVISION, area_m2 = area) %>%
    filter(!grepl(WATER_PATTERN, ecodiv)) %>%
    mutate(variant = tolower(variant_upper)) %>%
    group_by(variant) %>%
    mutate(weight = area_pct / sum(area_pct)) %>%
    ungroup() %>%
    transmute(
      variant,
      ecodiv,
      area_m2,
      area_pct,
      weight,
      source = "area"
    ) %>%
    arrange(variant, desc(weight))
}

build_proxied_coverage <- function(direct_cov, proxy_map) {
  ## For each proxied variant, copy the source variant's rows and tag.
  out <- list()
  for (v in names(proxy_map)) {
    src <- proxy_map[[v]]
    src_rows <- direct_cov %>% filter(variant == src)
    if (nrow(src_rows) == 0) {
      log_warn("Proxy source {src} not in direct coverage; skipping {v}")
      next
    }
    out[[v]] <- src_rows %>%
      mutate(
        variant = v,
        area_m2 = NA_real_,
        source  = paste0("proxy_", src)
      )
  }
  bind_rows(out)
}

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

validate_coverage <- function(cov_df) {
  ## Each variant's weight column must sum to 1.0 within rounding.
  sums <- cov_df %>% group_by(variant) %>% summarise(w = sum(weight), .groups = "drop")
  bad  <- sums %>% filter(abs(w - 1.0) > 1e-6)
  if (nrow(bad) > 0) {
    log_error("Weights do not sum to 1.0 for variants: {paste(bad$variant, collapse=', ')}")
    print(bad)
    stop("Weight normalization failed")
  }

  ## Coverage of all 25 variants.
  covered <- sort(unique(cov_df$variant))
  missing <- setdiff(ALL_VARIANTS, covered)
  if (length(missing) > 0) {
    log_warn("Missing variants from coverage: {paste(missing, collapse=', ')}")
  }
  unexpected <- setdiff(covered, ALL_VARIANTS)
  if (length(unexpected) > 0) {
    log_warn("Unexpected variants (not in ALL_VARIANTS): {paste(unexpected, collapse=', ')}")
  }

  ## Ecodiv code sanity: should be 3-digit Bailey code, optionally M-prefixed.
  bad_codes <- cov_df %>%
    filter(!grepl("^M?[0-9]{3}$", ecodiv)) %>%
    distinct(ecodiv)
  if (nrow(bad_codes) > 0) {
    log_warn("Non-standard ecodiv codes: {paste(bad_codes$ecodiv, collapse=', ')}")
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!file.exists(INPUT_MAP)) {
    stop("Input map not found: ", INPUT_MAP)
  }
  log_info("Reading input: {INPUT_MAP}")
  map_df <- read_csv(INPUT_MAP, col_types = cols())
  log_info("Input rows: {nrow(map_df)}; unique variants: {length(unique(map_df$FVSVariant))}")

  direct_cov  <- build_direct_coverage(map_df)
  proxied_cov <- build_proxied_coverage(direct_cov, PROXY_SOURCES)

  full_cov <- bind_rows(direct_cov, proxied_cov) %>%
    arrange(variant, desc(weight))

  log_info("Direct entries:  {nrow(direct_cov)}  ({length(unique(direct_cov$variant))} variants)")
  log_info("Proxied entries: {nrow(proxied_cov)} ({length(unique(proxied_cov$variant))} variants)")
  log_info("Total entries:   {nrow(full_cov)}  ({length(unique(full_cov$variant))} variants)")

  validate_coverage(full_cov)

  ## Per-variant summary table for the console log
  per_variant_summary <- full_cov %>%
    group_by(variant) %>%
    summarise(
      n_ecodiv = n(),
      dominant = ecodiv[1],
      dominant_pct = round(weight[1] * 100, 1),
      source = unique(source)[1],
      .groups = "drop"
    )
  message("\nPer-variant coverage summary:")
  print(as.data.frame(per_variant_summary), row.names = FALSE)

  if (args$dry_run) {
    log_info("[dry-run] would write {OUTPUT_CSV}")
    invisible(full_cov)
  } else {
    write_csv(full_cov, OUTPUT_CSV)
    log_info("Wrote {OUTPUT_CSV} ({nrow(full_cov)} rows)")
    invisible(full_cov)
  }
}

if (!interactive()) main()
