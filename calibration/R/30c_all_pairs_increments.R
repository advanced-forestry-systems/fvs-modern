## =============================================================================
## 30c_all_pairs_increments.R
##
## Expand a long-format per-tree FIA remeasurement table into all pairwise
## remeasurement combinations for fitting DG, HG, HCB, and mortality.
##
## Motivation:
##   David Marshall's CHANGE files currently emit only successive pairs
##   (t_k -> t_{k+1}), so a tree measured 3 times contributes 2 observations.
##   Adapting Salas & Riquelme (2018) creaIncrData we can emit n(n-1)/2
##   pairs, e.g. (t1,t2), (t1,t3), (t2,t3), which is a free sample-size
##   boost for any tree with 3+ measurements and especially valuable for
##   mortality (more exposure-years per tree) and long-interval size effects.
##
##   The Cao/Weiskittel likelihood already scales variance by sqrt(YEARS),
##   so longer-interval pairs are correctly down-weighted. Induced within-
##   tree correlation is absorbed by the nested EPA + tree-level RE structure.
##
## Input:
##   A long-format data.table keyed by `tree_key` with a `MEASYEAR` column
##   and time-varying attributes (DIA, HT, CR, STATUSCD, BA, BAL, BAL_SW,
##   BAL_HW, ... plus any plot / condition covariates). One row per tree per
##   measurement occasion.
##
## Output:
##   A "wide pair" data.table with columns DIA1 / DIA2, HT1 / HT2, CR1 / CR2,
##   STATUS1 / STATUS2, YEARS, plus all .1 / .2 suffixed time-varying
##   columns. Pair identity tracked by tree_key, MEASYEAR1, MEASYEAR2.
##
## Dependencies: gtools::combinations, data.table
##
## Usage:
##   long <- readRDS("calibration/data/conus_tree_long.rds")
##   pairs <- all_pairs_increments(
##     long,
##     tree_key  = "tree_key",
##     time_col  = "MEASYEAR",
##     tv_cols   = c("DIA", "HT", "CR", "STATUSCD",
##                   "BA", "BAL", "BAL_SW", "BAL_HW", "rd_add"),
##     static_cols = c("SPCD", "plot_key", "EPA_L1_CODE",
##                     "EPA_L2_CODE", "EPA_L3_CODE",
##                     "LAT", "LON", "cspi", "SDImax_brms"))
##
## Author: A. Weiskittel, 2026-04-15
## Based on Salas, C. & Riquelme, R. (2018) creaIncrData().
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(gtools)   # for combinations()
})

all_pairs_increments <- function(dat,
                                 tree_key     = "tree_key",
                                 time_col     = "MEASYEAR",
                                 tv_cols      = c("DIA", "HT", "CR", "STATUSCD",
                                                  "BA", "BAL", "BAL_SW", "BAL_HW"),
                                 static_cols  = c("SPCD", "plot_key",
                                                  "EPA_L1_CODE",
                                                  "EPA_L2_CODE",
                                                  "EPA_L3_CODE"),
                                 min_years    = 1,
                                 max_years    = 20,
                                 require_alive_start = TRUE,
                                 verbose      = TRUE) {
  stopifnot(is.data.frame(dat))
  setDT(dat)

  required <- c(tree_key, time_col, tv_cols, static_cols)
  miss <- setdiff(required, names(dat))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  ## Sort once; all downstream work is on subsets.
  setorderv(dat, c(tree_key, time_col))

  ## Per-tree measurement count. Trees with a single measurement are dropped.
  n_meas <- dat[, .N, by = tree_key]
  multi  <- n_meas[N >= 2, get(tree_key)]
  if (verbose) {
    message("Input trees:      ", format(nrow(n_meas), big.mark = ","))
    message("Multi-measure:    ", format(length(multi), big.mark = ","),
            "  (", round(100 * length(multi) / nrow(n_meas), 1), "%)")
    meas_dist <- n_meas[, .N, by = N][order(N)]
    message("Measurement count distribution:")
    print(meas_dist)
  }

  sub <- dat[get(tree_key) %in% multi]

  ## Build pair indices per tree using gtools::combinations on the sorted
  ## row positions within the tree. Vectorized at the tree level via by=.
  sub[, rownum := seq_len(.N), by = tree_key]

  pair_idx <- sub[, {
    if (.N < 2) {
      NULL
    } else {
      cmb <- gtools::combinations(n = .N, r = 2, v = seq_len(.N))
      list(i1 = cmb[, 1], i2 = cmb[, 2])
    }
  }, by = tree_key]

  if (verbose)
    message("Pair rows:        ", format(nrow(pair_idx), big.mark = ","))

  ## Join each pair end against the long table.
  col_suffix <- function(suffix, cols) setNames(cols, paste0(cols, suffix))

  pair_idx[, key_t1 := paste(get(tree_key), i1, sep = "::")]
  pair_idx[, key_t2 := paste(get(tree_key), i2, sep = "::")]

  sub[, key_t := paste(get(tree_key), rownum, sep = "::")]

  keep_cols <- c(time_col, tv_cols)

  a <- sub[, c("key_t", keep_cols), with = FALSE]
  b <- sub[, c("key_t", keep_cols), with = FALSE]

  setnames(a, "key_t", "key_t1")
  setnames(b, "key_t", "key_t2")
  setnames(a, keep_cols, paste0(keep_cols, "1"))
  setnames(b, keep_cols, paste0(keep_cols, "2"))

  pairs <- merge(pair_idx, a, by = "key_t1", all.x = TRUE)
  pairs <- merge(pairs,    b, by = "key_t2", all.x = TRUE)

  ## Attach static covariates (from the t1 row by convention).
  stat <- sub[, c(tree_key, "rownum", static_cols), with = FALSE]
  setnames(stat, "rownum", "i1")
  pairs <- merge(pairs, stat, by = c(tree_key, "i1"), all.x = TRUE)

  ## YEARS, standard naming, status filters.
  pairs[, YEARS := get(paste0(time_col, "2")) - get(paste0(time_col, "1"))]

  setnames(pairs,
           c(paste0(time_col, "1"), paste0(time_col, "2")),
           c("MEASYEAR1", "MEASYEAR2"))

  ## Standard downstream names used by script 31 / 32 / 34.
  tv_renames <- c(
    "DIA1"      = "DBH1",      "DIA2"      = "DBH2",
    "HT1"       = "HT1",       "HT2"       = "HT2",
    "CR1"       = "CR1",       "CR2"       = "CR2",
    "STATUSCD1" = "STATUS1",   "STATUSCD2" = "STATUS2",
    "BA1"       = "BA1",       "BA2"       = "BA2",
    "BAL1"      = "BAL1",      "BAL2"      = "BAL2",
    "BAL_SW1"   = "BAL_SW1",   "BAL_SW2"   = "BAL_SW2",
    "BAL_HW1"   = "BAL_HW1",   "BAL_HW2"   = "BAL_HW2"
  )
  rn <- tv_renames[intersect(names(tv_renames), names(pairs))]
  if (length(rn)) setnames(pairs, names(rn), rn)

  ## Interval filter and basic survival logic.
  before <- nrow(pairs)
  pairs <- pairs[YEARS >= min_years & YEARS <= max_years]
  if (require_alive_start && "STATUS1" %in% names(pairs)) {
    pairs <- pairs[STATUS1 == 1]
  }
  if (verbose)
    message("After filters:    ", format(nrow(pairs), big.mark = ","),
            " of ", format(before, big.mark = ","))

  ## Drop internal bookkeeping.
  pairs[, c("i1", "i2", "key_t1", "key_t2", "rownum") := NULL]

  pairs[]
}

## ---- If run directly, smoke-test on a tiny synthetic example ---------------
if (!interactive() && identical(sys.calls(), list())) {
  set.seed(42)
  nt  <- 50
  tree_key <- rep(sprintf("T%03d", seq_len(nt)), times = sample(1:4, nt, TRUE))
  n   <- length(tree_key)
  demo <- data.table(
    tree_key    = tree_key,
    MEASYEAR    = unlist(lapply(rle(tree_key)$lengths,
                                function(k) 2000 + sort(sample(0:20, k)))),
    DIA         = runif(n, 5, 60),
    HT          = runif(n, 5, 40),
    CR          = runif(n, 0.2, 0.8),
    STATUSCD    = 1L,
    BA          = runif(n, 10, 40),
    BAL         = runif(n,  0, 35),
    BAL_SW      = runif(n,  0, 20),
    BAL_HW      = runif(n,  0, 20),
    SPCD        = sample(c(202, 122, 316), n, TRUE),
    plot_key    = rep("P01", n),
    EPA_L1_CODE = "6",
    EPA_L2_CODE = "6.2",
    EPA_L3_CODE = "6.2.3"
  )
  out <- all_pairs_increments(demo, verbose = TRUE)
  cat("Smoke test produced", nrow(out), "pair rows\n")
  print(head(out, 3))
}
