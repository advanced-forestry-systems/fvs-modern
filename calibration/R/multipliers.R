# =============================================================================
# calibration/R/multipliers.R
#
# Comprehensive per-species calibration factors for all five tree component
# equations, written into the `calibration_multipliers` block that the FVS
# keyword emitter (config/config_loader.py) and downstream tooling consume.
#
# Components and how each species-specific factor is derived
# ----------------------------------------------------------
#   height_diameter (HD)  asymptote ratio: (a_pooled + r_SPCD__a[s]) / a_pooled
#                         from height_diameter_summary.csv (a is on the natural
#                         ft scale; H = 4.5 + a*(1-exp(-b*DBH))^c).
#   mortality (MORT)      mortality-rate ratio from the survival logit:
#                         (1 - plogis(b0 + r_SPCD[s])) / (1 - plogis(b0)).
#   crown_ratio (CR)      relative factor exp(r_SPCD_grouped[s]) from the CR
#                         change fit; an "OTHER" group is used for species the
#                         fit pooled into the catch-all.
#   diameter_growth (DG)  intercept-shift factor exp(b0[s] - mean b0). DG uses a
#                         dense Stan species index that is NOT SPCD-keyed and was
#                         not saved with a crosswalk, so DG factors are emitted
#                         ONLY where DG is adopted and are flagged approximate.
#   height_increment (HI) same dense-index situation as DG; gated identically.
#
# Species mapping
# ---------------
# HD, MORT, CR carry SPCD-keyed random effects (brms r_SPCD[...] terms), mapped
# to the FVS species slot through config categories.species_definitions.FIAJSP
# (FIA SPCD per FVS index). This is exact. DG/HI use a dense index and are
# best-effort only.
#
# Availability gating
# -------------------
# Every component is gated by the authoritative
# calibration/data/equation_availability_full.csv (columns variant, HD, MORT,
# CR, DG, SDI, HI). Where a component is FALSE for a variant, its factor array
# is all 1.0 (no-op) and the provenance records available = FALSE. This keeps
# the runtime config faithful to the published, adopted calibration: HD/MORT/CR
# are adopted for all 25 variants; DG for 7; HI for 6.
#
# All arrays are length maxsp, clipped to [lo, hi] for runtime safety.
# =============================================================================

.mult_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}

.mult_median_col <- function(df) {
  for (nm in c("p50", "median", "Estimate", "estimate", "mean")) {
    if (nm %in% names(df)) return(df[[nm]])
  }
  df[[2]]
}

.mult_clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# Parse SPCD-keyed brms random effects of the form  <prefix>[<SPCD>,<term>]
# Returns a named list: values = c(SPCD -> median), other = median for the
# "OTHER" catch-all group if present (else NA).
.mult_parse_re <- function(df, prefix) {
  out <- list(values = numeric(0), other = NA_real_)
  if (is.null(df) || nrow(df) == 0) return(out)
  v <- as.character(df[[1]])
  med <- .mult_median_col(df)
  # prefixes are plain identifiers (r_SPCD, r_SPCD__a, r_SPCD_grouped): no regex
  # specials, so match the literal prefix followed by a bracketed key.
  pat <- paste0("^\"?", prefix, "\\[")
  hit <- grepl(pat, v)
  if (!any(hit)) return(out)
  inside <- sub(paste0(".*", prefix, "\\[([^,\\]]+).*"), "\\1", v[hit], perl = TRUE)
  vals <- med[hit]
  oth <- which(toupper(trimws(inside)) == "OTHER")
  if (length(oth)) out$other <- vals[oth[1]]
  num <- suppressWarnings(as.integer(inside))
  keep <- !is.na(num)
  if (any(keep)) {
    nv <- vals[keep]; names(nv) <- as.character(num[keep])
    out$values <- nv
  }
  out
}

# Convert a dense-index intercept vector b0[i] (DG / HI) to a length-maxsp factor
# array. If a {component}_species_index.csv crosswalk (written by the fit script)
# is present, b0[i] is mapped exactly to the FVS species slot via SPCD; otherwise
# it falls back to placing b0[i] in slot i (approximate, flagged).
.mult_dense_factor <- function(b0, lo, hi, fia, maxsp, xwalk) {
  if (length(b0) == 0 || !any(!is.na(b0)))
    return(list(arr = rep(1.0, maxsp), mapping = "none", n = 0))
  mu <- mean(b0, na.rm = TRUE)
  m <- exp(b0 - mu); m[is.na(m)] <- 1.0
  arr <- rep(1.0, maxsp)
  ix <- if (file.exists(xwalk)) .mult_read_csv(xwalk) else NULL
  if (!is.null(ix) && all(c("species_idx", "SPCD") %in% names(ix))) {
    for (k in seq_len(nrow(ix))) {
      di <- suppressWarnings(as.integer(ix$species_idx[k]))
      sp <- suppressWarnings(as.integer(ix$SPCD[k]))
      if (!is.na(di) && di >= 1 && di <= length(m)) {
        slot <- which(fia == sp)
        if (length(slot)) arr[slot[1]] <- m[di]
      }
    }
    mapping <- "exact (species_index crosswalk)"
  } else {
    n <- min(length(m), maxsp); arr[seq_len(n)] <- m[seq_len(n)]
    mapping <- "approximate (no crosswalk; refit to save species index)"
  }
  list(arr = .mult_clip(arr, lo, hi), mapping = mapping, n = sum(!is.na(b0)))
}

# b0[i] dense-index intercept vector (DG / HI), NA where absent.
.mult_b0_vec <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(numeric(0))
  v <- as.character(df[[1]]); med <- .mult_median_col(df)
  hit <- grepl("^b0\\[[0-9]+\\]$", v)
  if (!any(hit)) return(numeric(0))
  idx <- as.integer(sub("^b0\\[([0-9]+)\\]$", "\\1", v[hit]))
  out <- rep(NA_real_, max(idx)); out[idx] <- med[hit]; out
}

# Load the availability table once; return named logical vector for a variant.
.mult_availability <- function(calibration_dir, variant) {
  f <- file.path(calibration_dir, "data", "equation_availability_full.csv")
  default <- c(HD = TRUE, MORT = TRUE, CR = TRUE, DG = TRUE, SDI = TRUE, HI = TRUE)
  tab <- .mult_read_csv(f)
  if (is.null(tab)) return(default)            # fail-open if table missing
  row <- tab[toupper(tab$variant) == toupper(variant), , drop = FALSE]
  if (nrow(row) == 0) return(default)
  as_lgl <- function(x) isTRUE(x) || (is.character(x) && toupper(x) == "TRUE")
  c(HD = as_lgl(row$HD), MORT = as_lgl(row$MORT), CR = as_lgl(row$CR),
    DG = as_lgl(row$DG), SDI = as_lgl(row$SDI), HI = as_lgl(row$HI))
}

#' Compute comprehensive per-species calibration factors for a variant.
#'
#' @param output_dir calibration/output/variants/<variant>
#' @param config parsed variant config (needs maxsp + species_definitions$FIAJSP)
#' @param calibration_dir calibration/ root (for the availability table)
#' @param lo,hi clip bounds
compute_calibration_multipliers <- function(output_dir, config,
                                             calibration_dir = NULL,
                                             lo = 0.1, hi = 10) {
  cats <- config$categories
  sd <- if (!is.null(cats)) cats$species_definitions else NULL
  fia <- if (!is.null(sd$FIAJSP)) suppressWarnings(as.integer(sd$FIAJSP)) else integer(0)
  maxsp <- if (!is.null(config$maxsp)) as.integer(config$maxsp) else length(fia)
  if (is.na(maxsp) || maxsp < 1) maxsp <- max(length(fia), 1L)
  variant <- if (!is.null(config$variant)) config$variant else ""
  if (is.null(calibration_dir)) calibration_dir <- dirname(dirname(output_dir))
  avail <- .mult_availability(calibration_dir, variant)

  # SPCD-keyed factor -> length-maxsp array on the FVS species slots.
  spcd_to_array <- function(re, default_other = NA_real_) {
    arr <- rep(1.0, maxsp)
    if (length(re$values) == 0 && is.na(re$other)) return(arr)
    for (i in seq_len(maxsp)) {
      if (i > length(fia)) next
      key <- as.character(fia[i])
      if (!is.na(re$values[key])) arr[i] <- re$values[[key]]
      else if (!is.na(re$other)) arr[i] <- re$other
    }
    arr
  }

  prov <- list(
    variant = variant, maxsp = maxsp, clip = c(lo, hi),
    available = as.list(avail),
    definitions = list(
      htdbh_multiplier = "HD asymptote ratio (a_pooled + r_SPCD__a)/a_pooled",
      mort_multiplier  = "mortality-rate ratio from survival logit incl species RE",
      cr_multiplier    = "exp(r_SPCD_grouped) relative crown-ratio factor",
      dds_multiplier   = "DG intercept-shift exp(delta b0); dense index, approximate, adopted variants only",
      htg_multiplier   = "HI intercept-shift exp(delta b0); dense index, approximate, adopted variants only"
    )
  )

  ## ---- height_diameter (HD): asymptote ratio, SPCD-keyed ----
  htdbh <- rep(1.0, maxsp)
  if (isTRUE(avail["HD"])) {
    hd <- .mult_read_csv(file.path(output_dir, "height_diameter_summary.csv"))
    if (!is.null(hd)) {
      a_pooled <- suppressWarnings(as.numeric(
        .mult_median_col(hd)[match("b_a_Intercept", hd[[1]])]))
      re <- .mult_parse_re(hd, "r_SPCD__a")
      if (!is.na(a_pooled) && a_pooled != 0) {
        # convert additive species deviation on a to a height-asymptote ratio
        re$values <- (a_pooled + re$values) / a_pooled
        if (!is.na(re$other)) re$other <- (a_pooled + re$other) / a_pooled
        htdbh <- .mult_clip(spcd_to_array(re), lo, hi)
        prov$htdbh_n_species <- length(re$values)
      }
    }
  }

  ## ---- mortality (MORT): survival-logit rate ratio, SPCD-keyed ----
  mort <- rep(1.0, maxsp)
  if (isTRUE(avail["MORT"])) {
    msum <- .mult_read_csv(file.path(output_dir, "mortality_summary.csv"))
    mpost <- .mult_read_csv(file.path(output_dir, "mortality_posterior.csv"))
    if (!is.null(msum)) {
      icpt <- suppressWarnings(as.numeric(
        .mult_median_col(msum)[match("Intercept", msum[[1]])]))
      if (!is.na(icpt)) {
        p_base <- stats::plogis(icpt)
        re <- .mult_parse_re(mpost, "r_SPCD")
        re$values <- (1 - stats::plogis(icpt + re$values)) / (1 - p_base)
        if (!is.na(re$other)) re$other <- (1 - stats::plogis(icpt + re$other)) / (1 - p_base)
        mort <- .mult_clip(spcd_to_array(re), lo, hi)
        prov$mort_n_species <- length(re$values)
      }
    }
  }

  ## ---- crown_ratio (CR): relative factor, SPCD-keyed ----
  cr <- rep(1.0, maxsp)
  if (isTRUE(avail["CR"])) {
    crp <- .mult_read_csv(file.path(output_dir, "crown_ratio_v2_posterior.csv"))
    if (is.null(crp)) crp <- .mult_read_csv(file.path(output_dir, "crown_ratio_posterior.csv"))
    if (!is.null(crp)) {
      re <- .mult_parse_re(crp, "r_SPCD_grouped")
      if (length(re$values) == 0) re <- .mult_parse_re(crp, "r_SPCD")
      re$values <- exp(re$values)
      if (!is.na(re$other)) re$other <- exp(re$other)
      cr <- .mult_clip(spcd_to_array(re), lo, hi)
      prov$cr_n_species <- length(re$values)
    }
  }

  ## ---- diameter_growth (DG): dense index, adopted variants only ----
  dds <- rep(1.0, maxsp)
  if (isTRUE(avail["DG"])) {
    dg <- .mult_read_csv(file.path(output_dir, "diameter_growth_summary.csv"))
    res <- .mult_dense_factor(.mult_b0_vec(dg), lo, hi, fia, maxsp,
                              file.path(output_dir, "diameter_growth_species_index.csv"))
    dds <- res$arr; prov$dds_n_species <- res$n; prov$dds_mapping <- res$mapping
  }

  ## ---- height_increment (HI): dense index, adopted variants only ----
  htg <- rep(1.0, maxsp)
  if (isTRUE(avail["HI"])) {
    hg <- .mult_read_csv(file.path(output_dir, "height_increment_summary.csv"))
    res <- .mult_dense_factor(.mult_b0_vec(hg), lo, hi, fia, maxsp,
                              file.path(output_dir, "height_increment_species_index.csv"))
    htg <- res$arr; prov$htg_n_species <- res$n; prov$htg_mapping <- res$mapping
  }

  list(
    htdbh_multiplier = as.numeric(htdbh),
    mort_multiplier  = as.numeric(mort),
    cr_multiplier    = as.numeric(cr),
    dds_multiplier   = as.numeric(dds),
    htg_multiplier   = as.numeric(htg),
    provenance       = prov
  )
}
