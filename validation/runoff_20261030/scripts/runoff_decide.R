#!/usr/bin/env Rscript
## runoff_decide.R GATES_BY_ARM NEAR_JOINED OUTDIR
## Pre-registered decision rule (PREREG_runoff_20261030.md section 7), applied per variant:
##  1. Admissible set = arms with admissible == TRUE in catalog_gates.R output on the run-off stands
##     (BA and SDI p99.9 exceedance by 50 yr < 5% with Wilson upper < 10%, Reineke flag rate < 10%,
##     median peak BA >= FIA BA p50).
##  2. Among admissible arms, best = smallest |near-term BA percent bias| on the common stand set.
##  3. Tie = an admissible arm whose paired stand-bootstrap 95% interval of (|pbias_arm| - |pbias_best|)
##     contains 0 (2000 resamples, seed 20261028). Among tied arms the simplest is chosen:
##     usfs_native < default_native < default_gompit < calibrated_native < calibrated_gompit < conus_b1
##     < conus_b1_gompit < conus_b2 < conus_b2_gompit.
##  4. No admissible arm: variant NOT RECOMMENDED for long-horizon use; the near-term best is named for
##     single-cycle use only.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); G <- fread(a[1]); P <- fread(a[2]); OUT <- a[3]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NB <- 2000; set.seed(20261028)
simp <- c(usfs_native = 0, "usfs_native(ne_lib)" = 0, default_native = 10, default_gompit = 11, calibrated_native = 20, calibrated_gompit = 21,
          conus_b1 = 30, conus_b1_gompit = 31, conus_b2 = 40, conus_b2_gompit = 41,
          conus_v2_b1 = 30, conus_v2_b1_gompit = 31, conus_v2_b2 = 40, conus_v2_b2_gompit = 41,
          conus_v6_b1 = 50, conus_v6_b1_gompit = 51, conus_v6_b2 = 60, conus_v6_b2_gompit = 61)
Q <- P[common == TRUE]; rows <- list()
for (vv in sort(unique(Q$variant))) {
  d <- Q[variant == vv]; W <- dcast(d, traj ~ arm, value.var = "BA"); O <- d[, .(BAo = BAo[1]), by = traj][W, on = "traj"]
  arms <- setdiff(names(W), "traj"); n <- nrow(O)
  pb <- function(i) sapply(arms, function(ar) 100 * mean(O[[ar]][i] - O$BAo[i]) / mean(O$BAo[i]))
  est <- pb(seq_len(n)); bs <- replicate(NB, pb(sample.int(n, replace = TRUE))); if (is.null(dim(bs))) bs <- matrix(bs, nrow = 1, dimnames = list(arms, NULL))
  adm <- G[toupper(group) == vv & arm %in% arms, .(arm, admissible = as.logical(admissible))]
  aset <- adm[admissible %in% TRUE, arm]
  ## Amendment A2: v6 projector arms are diagnostic only (stand BA cap suspected miscalibrated); never eligible.
  diag_arms <- grep("^conus_v6_", arms, value = TRUE); aset <- setdiff(aset, diag_arms)
  if (length(aset)) {
    best <- aset[which.min(abs(est[aset]))]
    tied <- aset[sapply(aset, function(ar) { dd <- abs(bs[ar, ]) - abs(bs[best, ]); q <- quantile(dd, c(.025, .975)); q[1] <= 0 & q[2] >= 0 })]
    chosen <- tied[which.min(ifelse(is.na(simp[tied]), 99, simp[tied]))]; path <- "admissible -> near-term BA bias -> simplicity among ties"
  } else {
    elig <- setdiff(arms, diag_arms); best <- elig[which.min(abs(est[elig]))]; tied <- character(); chosen <- NA_character_; path <- "NO ADMISSIBLE ARM: not recommended for long horizon"
  }
  ci <- apply(bs, 1, quantile, c(.025, .975))
  rows[[length(rows) + 1]] <- data.table(variant = vv, n_common = n, arms = paste(arms, collapse = ";"),
    admissible = paste(aset, collapse = ";"), best_by_bias = best, tied = paste(tied, collapse = ";"), chosen = chosen, path = path,
    ba_pbias = paste(sprintf("%s %.1f [%.1f, %.1f]", arms, est, ci[1, ], ci[2, ]), collapse = "; "),
    diagnostic = paste(sprintf("%s (diagnostic: stand BA cap suspected miscalibrated; admissible=%s)", diag_arms,
      adm[match(diag_arms, adm$arm), admissible]), collapse = "; "))
}
D <- rbindlist(rows); fwrite(D, file.path(OUT, "runoff_decision_by_variant.csv")); print(D)
