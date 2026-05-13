## =============================================================================
## 40b_enrich_cond_modifiers.R
##
## Phase 2 COND-level enrichment, revised. The original
## 40_enrich_cond_modifiers.R assumed PLT_CN1 and CONDID1 were carried on
## the pairs file, but the current metric RDS does not have them. This
## script joins on the composite FIA plot-measurement key that IS present:
## (STATECD, UNITCD, COUNTYCD, PLOT, INVYR1).
##
## Outputs:
##   conus_remeasurement_pairs_metric_cond.rds
##
## Also derives compact binary / time-since modifier variables that the
## Phase 2 Stan models can consume directly:
##
##   is_plantation                 STDORGCD == 1
##   had_fire_t1, had_insect_t1,
##   had_disease_t1, had_wind_t1,
##   had_harvest_t1                DSTRBCD1 family indicators at t1
##   had_cutting_t1,
##   had_site_prep_t1              TRTCD1 indicators at t1
##   years_since_dstrb             INVYR1 - DSTRBYR1, clipped to [0, 50]
##   years_since_trt               INVYR1 - TRTYR1,   clipped to [0, 50]
##   dstrb_decay_5yr, dstrb_decay_10yr, dstrb_decay_20yr
##                                 exp(-years_since_dstrb / lambda) time-
##                                 since decay envelopes at lambda = 5, 10,
##                                 20 for use as modifier covariates.
##   trt_decay_5yr, trt_decay_10yr, trt_decay_20yr
##                                 same, for treatment.
##
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opts <- OptionParser(option_list = list(
  make_option("--pairs", type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric.rds",
              dest    = "pairs_path"),
  make_option("--cond",  type = "character",
              default = "calibration/data/fia_cond.rds",
              dest    = "cond_path"),
  make_option("--out",   type = "character",
              default = "calibration/data/conus_remeasurement_pairs_metric_cond.rds",
              dest    = "out_path")
)) |> parse_args()

cat("==============================================================\n")
cat("FVS-CONUS Phase 2 COND enrichment (composite-key join)\n")
cat("  pairs :", opts$pairs_path, "\n")
cat("  cond  :", opts$cond_path, "\n")
cat("  out   :", opts$out_path, "\n")
cat("==============================================================\n")

stopifnot(file.exists(opts$pairs_path))
stopifnot(file.exists(opts$cond_path))

pairs <- readRDS(opts$pairs_path)
setDT(pairs)
cat(sprintf("Pairs: %s rows, %d cols\n",
            format(nrow(pairs), big.mark = ","), ncol(pairs)))

cond <- readRDS(opts$cond_path)
setDT(cond)
cat(sprintf("COND : %s rows, %d cols\n",
            format(nrow(cond), big.mark = ","), ncol(cond)))

## ---- 1. Align keys --------------------------------------------------------

key_pairs <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR1")
key_cond  <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR")

missing_p <- setdiff(key_pairs, names(pairs))
missing_c <- setdiff(key_cond,  names(cond))
if (length(missing_p) > 0)
  stop("pairs missing keys: ", paste(missing_p, collapse = ", "))
if (length(missing_c) > 0)
  stop("cond missing keys: ",  paste(missing_c, collapse = ", "))

## Rename COND -> INVYR1 for the join and suffix non-key columns
setnames(cond, "INVYR", "INVYR1")
suffix_cols <- setdiff(names(cond), key_pairs)
new_names   <- paste0(suffix_cols, "_cond1")
setnames(cond, suffix_cols, new_names)

cat(sprintf("Joining %d COND columns on %s\n",
            length(new_names), paste(key_pairs, collapse = ", ")))

pairs_out <- cond[pairs, on = key_pairs]  # right join: all pairs rows retained
setkeyv(pairs_out, NULL)

cat(sprintf("  Joined: %s rows (pairs had %s)\n",
            format(nrow(pairs_out), big.mark = ","),
            format(nrow(pairs),     big.mark = ",")))

match_rate <- 100 * mean(!is.na(pairs_out$STDORGCD_cond1))
cat(sprintf("  Match rate (STDORGCD non-NA): %.1f%%\n", match_rate))

## ---- 2. Derive compact modifier variables --------------------------------

if ("STDORGCD_cond1" %in% names(pairs_out)) {
  pairs_out[, is_plantation := fifelse(STDORGCD_cond1 == 1, 1L, 0L,
                                       na = NA_integer_)]
}
if ("DSTRBCD1_cond1" %in% names(pairs_out)) {
  pairs_out[, had_fire_t1    := as.integer(DSTRBCD1_cond1 == 30)]
  pairs_out[, had_insect_t1  := as.integer(DSTRBCD1_cond1 %in% 10:19)]
  pairs_out[, had_disease_t1 := as.integer(DSTRBCD1_cond1 %in% 20:29)]
  pairs_out[, had_wind_t1    := as.integer(DSTRBCD1_cond1 == 50)]
  pairs_out[, had_harvest_t1 := as.integer(DSTRBCD1_cond1 == 80)]
}
if ("TRTCD1_cond1" %in% names(pairs_out)) {
  pairs_out[, had_cutting_t1   := as.integer(TRTCD1_cond1 == 10)]
  pairs_out[, had_natregen_t1  := as.integer(TRTCD1_cond1 == 20)]
  pairs_out[, had_artregen_t1  := as.integer(TRTCD1_cond1 == 30)]
  pairs_out[, had_site_prep_t1 := as.integer(TRTCD1_cond1 == 40)]
}

## Time-since variables. Use INVYR1 as the reference year. Clip to
## [0, 50] so that very old / invalid years do not anchor the decay at
## implausible values. Retained NA when DSTRBYR / TRTYR is NA or <=0.

clip_yst <- function(inv_yr, ev_yr) {
  y <- inv_yr - ev_yr
  y[is.na(y) | ev_yr <= 0] <- NA_real_
  pmin(pmax(y, 0), 50)
}

if ("DSTRBYR1_cond1" %in% names(pairs_out)) {
  pairs_out[, years_since_dstrb := clip_yst(INVYR1, DSTRBYR1_cond1)]
}
if ("TRTYR1_cond1" %in% names(pairs_out)) {
  pairs_out[, years_since_trt := clip_yst(INVYR1, TRTYR1_cond1)]
}

## Exponential time-since decay envelopes at three time constants.
## These are ready-made modifier covariates:
##   delta_dist_i = beta_dstrb * dstrb_decay_<lambda>_i * indicator_i
decay_envelope <- function(yst, lambda) {
  ifelse(is.na(yst), 0, exp(-yst / lambda))
}
if ("years_since_dstrb" %in% names(pairs_out)) {
  pairs_out[, dstrb_decay_5yr  := decay_envelope(years_since_dstrb,  5)]
  pairs_out[, dstrb_decay_10yr := decay_envelope(years_since_dstrb, 10)]
  pairs_out[, dstrb_decay_20yr := decay_envelope(years_since_dstrb, 20)]
}
if ("years_since_trt" %in% names(pairs_out)) {
  pairs_out[, trt_decay_5yr  := decay_envelope(years_since_trt,  5)]
  pairs_out[, trt_decay_10yr := decay_envelope(years_since_trt, 10)]
  pairs_out[, trt_decay_20yr := decay_envelope(years_since_trt, 20)]
}

## ---- 3. Coverage / QC summaries -------------------------------------------

cat("\nBinary modifier coverage and prevalence:\n")
for (v in c("is_plantation", "had_fire_t1", "had_insect_t1",
            "had_disease_t1", "had_wind_t1", "had_harvest_t1",
            "had_cutting_t1", "had_natregen_t1", "had_artregen_t1",
            "had_site_prep_t1")) {
  if (v %in% names(pairs_out)) {
    nn  <- sum(!is.na(pairs_out[[v]]))
    mn  <- mean(pairs_out[[v]] == 1L, na.rm = TRUE)
    cat(sprintf("  %-18s n=%s  prev=%.3f\n",
                v, format(nn, big.mark = ","), mn))
  }
}

cat("\nTime-since decay ranges (non-NA rows):\n")
for (v in c("years_since_dstrb", "years_since_trt",
            "dstrb_decay_5yr", "dstrb_decay_10yr", "dstrb_decay_20yr",
            "trt_decay_5yr", "trt_decay_10yr", "trt_decay_20yr")) {
  if (v %in% names(pairs_out)) {
    x <- pairs_out[[v]]
    nn <- sum(!is.na(x))
    rng <- range(x, na.rm = TRUE)
    cat(sprintf("  %-20s n=%s  min=%.3f  max=%.3f\n",
                v, format(nn, big.mark = ","), rng[1], rng[2]))
  }
}

cat("\nDisturbance / treatment prevalence by EPA L1:\n")
if ("EPA_L1_CODE" %in% names(pairs_out)) {
  prev_l1 <- pairs_out[!is.na(is_plantation),
                       .(n_pairs       = .N,
                         pct_plant     = 100 * mean(is_plantation),
                         pct_dstrb_any = 100 * mean(!is.na(DSTRBCD1_cond1) &
                                                     DSTRBCD1_cond1 > 0),
                         pct_trt_any   = 100 * mean(!is.na(TRTCD1_cond1) &
                                                     TRTCD1_cond1 > 0)),
                       by = EPA_L1_CODE][order(EPA_L1_CODE)]
  print(prev_l1)
}

cat(sprintf("\nWriting %s ...\n", opts$out_path))
saveRDS(pairs_out, opts$out_path, compress = "xz")
sz <- file.info(opts$out_path)$size
cat(sprintf("Done. %s rows, %d cols, %.1f MB\n",
            format(nrow(pairs_out), big.mark = ","),
            ncol(pairs_out), sz / 1e6))
