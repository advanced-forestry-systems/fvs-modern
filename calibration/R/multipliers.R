# =============================================================================
# calibration/R/multipliers.R
#
# Precompute schema-free, per-species calibration multipliers that FVS keyword
# records (BAIMULT, MORTMULT, HTGMULT) consume directly.
#
# Why this exists
# ---------------
# The Python keyword emitter (config/config_loader.py :: generate_keywords)
# previously reverse-engineered raw coefficients (B0, MORT_B0, ...) out of the
# calibrated JSON to derive growth/mortality/height multipliers. The writer
# (06_posterior_to_json.R) and the emitter drifted apart: coefficients were
# written into FVS-native slots the emitter never read, so only the SDIMAX block
# survived and at most one keyword block per variant reached FVS at runtime.
#
# This module moves the multiplier computation to where the model knowledge
# already lives (R), and serializes finished per-species arrays. The emitter
# then just formats them, with no per-variant schema coupling.
#
# Definition (per species, relative to the calibrated population/pooled baseline
# of each component fit):
#   diameter growth   dds_multiplier[s]  = exp(b0_dg[s] - mu_b0_dg)   (DDS scale)
#   height growth      htg_multiplier[s] = exp(b0_hg[s] - mu_b0_hg)
#   mortality          mort_multiplier[s] = (1 - p_s) / (1 - p_base),
#                        p = plogis(intercept [+ r_SPCD[s]])  (survival logit)
#
# All arrays have length maxsp, padded with 1.0 (no-op) for species a given
# component did not fit, and clipped to [lo, hi] for runtime safety.
#
# NOTE ON BASELINE: these express each species RELATIVE TO THE CALIBRATED
# POPULATION MEAN of its own fit, i.e. the species-differentiation layer. A
# separate global calibrated-vs-FVS-default shift can be layered on later once a
# default-model baseline is wired in. That modeling choice is isolated to this
# one function on purpose (see issue #54, decision points 2 and 3).
# =============================================================================

.mult_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

# Pull the posterior median column, robust to summary-file naming conventions.
.mult_median_col <- function(df) {
  for (nm in c("p50", "median", "Estimate", "estimate", "mean")) {
    if (nm %in% names(df)) return(df[[nm]])
  }
  df[[2]]
}

# Extract b0[i] (1-indexed) medians into a length-n vector (NA where absent).
.mult_b0_vec <- function(df, n) {
  out <- rep(NA_real_, n)
  if (is.null(df) || nrow(df) == 0) return(out)
  var <- as.character(df[[1]])
  med <- .mult_median_col(df)
  hit <- grepl("^b0\\[[0-9]+\\]$", var)
  if (!any(hit)) return(out)
  idx <- as.integer(sub("^b0\\[([0-9]+)\\]$", "\\1", var[hit]))
  vals <- med[hit]
  for (k in seq_along(idx)) {
    i <- idx[k]
    if (!is.na(i) && i >= 1 && i <= n) out[i] <- vals[k]
  }
  out
}

.mult_clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' Compute per-species calibration multipliers for the keyword emitter.
#'
#' @param output_dir calibration/output/variants/<variant> directory holding the
#'   *_summary.csv and *_map.csv posterior files.
#' @param config the variant config list (as parsed from config/<variant>.json or
#'   the calibrated JSON); must carry maxsp and categories$species_definitions.
#' @param lo,hi multiplier clip bounds.
#' @return list(dds_multiplier, htg_multiplier, mort_multiplier, provenance).
compute_calibration_multipliers <- function(output_dir, config, lo = 0.1, hi = 10) {
  cats <- config$categories
  sd <- if (!is.null(cats)) cats$species_definitions else NULL
  fia <- if (!is.null(sd$FIAJSP)) suppressWarnings(as.integer(sd$FIAJSP)) else integer(0)
  maxsp <- if (!is.null(config$maxsp)) as.integer(config$maxsp) else length(fia)
  if (is.na(maxsp) || maxsp < 1) maxsp <- length(fia)
  if (maxsp < 1) maxsp <- 1

  prov <- list(
    baseline = "calibrated population/pooled mean of each component fit",
    scale = list(
      dds_multiplier  = "DDS scale exp(delta b0); emitter applies sqrt for diameter-growth scale",
      htg_multiplier  = "height-increment ratio exp(delta b0)",
      mort_multiplier = "mortality rate ratio from survival logit incl species random effect"
    ),
    clip = c(lo, hi),
    maxsp = maxsp
  )

  ## ---- diameter growth ----
  # Baseline is the EMPIRICAL mean of the fitted species intercepts, not the
  # hierarchical hyperparameter mu_b0 from the *_map.csv. The map mu_b0 can sit
  # on a non-centered / transformed scale (b0[i] = mu + sigma * z), which does
  # not line up with the reconstructed b0[i] in the summary for some components
  # (notably height increment). The empirical species mean guarantees the
  # multipliers center on 1.0 by construction.
  dds <- rep(1.0, maxsp)
  dg <- .mult_read_csv(file.path(output_dir, "diameter_growth_summary.csv"))
  if (!is.null(dg)) {
    b0 <- .mult_b0_vec(dg, maxsp)
    if (any(!is.na(b0))) {
      mu <- mean(b0, na.rm = TRUE)
      m <- exp(b0 - mu); m[is.na(m)] <- 1.0
      dds <- .mult_clip(m, lo, hi)
    }
    prov$dds_n_species <- sum(!is.na(b0))
  }

  ## ---- height growth (height increment) ----
  htg <- rep(1.0, maxsp)
  hg <- .mult_read_csv(file.path(output_dir, "height_increment_summary.csv"))
  if (!is.null(hg)) {
    b0 <- .mult_b0_vec(hg, maxsp)
    if (any(!is.na(b0))) {
      mu <- mean(b0, na.rm = TRUE)
      m <- exp(b0 - mu); m[is.na(m)] <- 1.0
      htg <- .mult_clip(m, lo, hi)
    }
    prov$htg_n_species <- sum(!is.na(b0))
  }

  ## ---- mortality (survival logit + species random effect), SPCD -> FVS index ----
  mort <- rep(1.0, maxsp)
  msum <- .mult_read_csv(file.path(output_dir, "mortality_summary.csv"))
  mpost <- .mult_read_csv(file.path(output_dir, "mortality_posterior.csv"))
  if (!is.null(msum) && length(fia) > 0) {
    icpt <- suppressWarnings(as.numeric(.mult_median_col(msum)[match("Intercept", msum[[1]])]))
    if (!is.na(icpt)) {
      p_base <- stats::plogis(icpt)
      re_spcd <- integer(0); re_val <- numeric(0)
      if (!is.null(mpost)) {
        rv <- as.character(mpost[[1]])
        rm <- .mult_median_col(mpost)
        hit <- grepl("r_SPCD\\[[0-9]+", rv)
        if (any(hit)) {
          re_spcd <- as.integer(sub('^[^0-9]*r_SPCD\\[([0-9]+).*$', "\\1", rv[hit]))
          re_val <- rm[hit]
        }
      }
      for (i in seq_len(maxsp)) {
        r <- 0.0
        if (i <= length(fia) && length(re_spcd) > 0) {
          j <- match(fia[i], re_spcd)
          if (!is.na(j)) r <- re_val[j]
        }
        p_s <- stats::plogis(icpt + r)
        mort[i] <- (1 - p_s) / (1 - p_base)
      }
      mort <- .mult_clip(mort, lo, hi)
      prov$mort_n_re <- length(re_spcd)
      prov$mort_intercept <- icpt
    }
  }

  list(
    dds_multiplier  = as.numeric(dds),
    htg_multiplier  = as.numeric(htg),
    mort_multiplier = as.numeric(mort),
    provenance      = prov
  )
}
