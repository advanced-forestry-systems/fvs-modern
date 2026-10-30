#!/usr/bin/env Rscript
## runoff_nearterm_obs.R (2026-10-30, v4)
## Observed time-1 and time-2 stand metrics for every plot in the catalog sample frame
## (catalog/fia_variant_assignment.csv), from the CONUS remeasurement pairs RDS. The pairs file is
## tree-level and a plot_key (STATECD-UNITCD-COUNTYCD-PLOT) can hold several remeasurement periods,
## so rows are restricted to the period the catalog seeded from (INVYR1 of the assignment row).
## Units (checked in obsdiag job 16012770): DBH1/DBH2 inches, HT1/HT2 metres, TPH_UNADJ1 trees ha-1
## per tree record (= TPA_UNADJ1 x 2.471). Expansion = TPH_UNADJ1 / CONDPROP_UNADJ_cond1, the same
## condition-acre expansion the catalog applies when seeding (TPA_UNADJ / CONDPROP_UNADJ).
## Time-2 metrics are SURVIVOR metrics (live at time 1 and at time 2); the pairs file carries no
## ingrowth and no contestant regenerates, so this is the like-for-like target.
## Output units: BA m2 ha-1, TPH ha-1, QMD cm, top height m (TPA-weighted mean HT of the 40 largest-DBH
## trees per acre = 98.84 per ha, the FVS ATTOPHT rule), VOLC m3 ha-1 = 0.40 x BA x top height.
## No PLT_CN or coordinates are written.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); PAIRS <- a[1]; ASSIGN <- a[2]; OUT <- a[3]
asg <- fread(ASSIGN, select = c(1, 3, 4, 5, 6, 7, 9)); setnames(asg, c("plot_key","variant","STATECD","COUNTYCD","PLOT","INVYR1","CONDID"))
x <- as.data.table(readRDS(PAIRS))
x <- x[asg[, .(plot_key, INVYR1)], on = .(plot_key, INVYR1), nomatch = NULL]
cat("time-1 tree rows with no time-2 record (STATUS2 NA):", x[!is.finite(STATUS2), .N], "\n")
cat("pairs tree rows for frame plot x INVYR1:", nrow(x), " plots:", uniqueN(x$plot_key), " of", nrow(asg), "\n")
cat("distinct non-missing INVYR2 per plot_key (should be 1):\n"); print(table(x[, uniqueN(INVYR2[is.finite(INVYR2)]), by = plot_key]$V1))
x[, cp := fifelse(is.finite(CONDPROP_UNADJ_cond1) & CONDPROP_UNADJ_cond1 > 0, CONDPROP_UNADJ_cond1, 1)]
x[, e := fifelse(is.finite(TPH_UNADJ1), TPH_UNADJ1, 0)]; x[, ecp := e / cp]
TOPN <- 40 * 2.4710538
topht <- function(d, h, w) { ok <- is.finite(d) & is.finite(h) & is.finite(w) & w > 0; d <- d[ok]; h <- h[ok]; w <- w[ok]
  if (!length(d)) return(NA_real_); o <- order(d, decreasing = TRUE); d <- d[o]; h <- h[o]; w <- w[o]
  cw <- cumsum(w); take <- pmin(w, pmax(TOPN - (cw - w), 0)); if (sum(take) <= 0) return(NA_real_); sum(h * take) / sum(take) }
bam <- function(d) pi / 4 * (d * 0.0254)^2
g <- x[, { l2 <- (STATUS2 == 1) %in% TRUE & is.finite(DBH2)
  .(YEARS = if (any(is.finite(YEARS))) as.numeric(YEARS[is.finite(YEARS)][1]) else NA_real_, INVYR2 = if (any(is.finite(INVYR2))) as.numeric(max(INVYR2, na.rm = TRUE)) else NA_real_, unk_ba_share = sum((ecp * bam(DBH1))[!is.finite(STATUS2)]) / sum(ecp * bam(DBH1)), n_t1 = .N, n_surv = sum(l2), any_t2 = any(is.finite(STATUS2)),
    BA1_file = first(BA1), BA2_file = first(BA2), TPH1_file = first(TPH1), TPH2_file = first(TPH2), HT40_1_file = first(HT40_1),
    BA1_nocp = sum(e * bam(DBH1)), TPH1_nocp = sum(e),
    BA1 = sum(ecp * bam(DBH1)), TPH1 = sum(ecp), BA2 = sum((ecp * bam(DBH2))[l2]), TPH2 = sum(ecp[l2]),
    HTTOP1 = topht(DBH1, HT1, ecp), HTTOP2 = topht(DBH2[l2], HT2[l2], ecp[l2])) }, by = plot_key]
cat("\nagreement with the file's stand columns (median ratio; share within 2 pct):\n")
print(g[is.finite(BA1_file) & BA1_file > 0, .(BA1_cp = median(BA1 / BA1_file), BA1_cp_w2 = mean(abs(BA1 / BA1_file - 1) < .02),
  BA1_nocp = median(BA1_nocp / BA1_file), BA1_nocp_w2 = mean(abs(BA1_nocp / BA1_file - 1) < .02),
  TPH1_cp = median(TPH1 / TPH1_file, na.rm = TRUE), TPH1_nocp = median(TPH1_nocp / TPH1_file, na.rm = TRUE),
  BA2_surv = median(BA2 / BA2_file, na.rm = TRUE), BA2_surv_le1 = mean(BA2 <= BA2_file * 1.001, na.rm = TRUE), BA2_surv_w2 = mean(abs(BA2 / BA2_file - 1) < .02, na.rm = TRUE),
  TPH2_surv = median(TPH2 / TPH2_file, na.rm = TRUE), HTTOP1 = median(HTTOP1 / HT40_1_file, na.rm = TRUE))])
g[, `:=`(QMD1 = sqrt(BA1 / TPH1 / (pi / 40000)), QMD2 = fifelse(TPH2 > 0, sqrt(BA2 / TPH2 / (pi / 40000)), NA_real_))]
g[, `:=`(VOLC1 = 0.40 * BA1 * HTTOP1, VOLC2 = 0.40 * BA2 * HTTOP2)]
g[, interval := fifelse(is.finite(YEARS) & YEARS > 0, YEARS, NA_real_)]
g[, k_cycle := as.integer(round(interval))]
o <- merge(asg, g, by = "plot_key", all.x = TRUE)
o[is.na(interval) & is.finite(INVYR2), `:=`(interval = as.numeric(INVYR2 - INVYR1), k_cycle = as.integer(INVYR2 - INVYR1))]
o[, obs_ok := any_t2 %in% TRUE & is.finite(interval) & interval > 0 & is.finite(BA2) & k_cycle >= 1 & unk_ba_share <= 0.05]
cat("plots excluded for >5 pct of time-1 BA with unknown time-2 fate:", o[any_t2 %in% TRUE & unk_ba_share > 0.05, .N], "\n")
fwrite(o[, .(plot_key, variant, STATECD, COUNTYCD, PLOT, INVYR1, CONDID, INVYR2, interval, k_cycle, obs_ok, n_t1, n_surv,
             BA1, TPH1, QMD1, HTTOP1, VOLC1, BA2, TPH2, QMD2, HTTOP2, VOLC2, BA1_file, BA2_file, TPH1_file, TPH2_file, unk_ba_share)], OUT)
cat("\nwrote", OUT, nrow(o), "plots; obs_ok", sum(o$obs_ok), "\n")
print(o[, .(n = .N, obs_ok = sum(obs_ok), k_med = median(k_cycle, na.rm = TRUE), k_min = min(k_cycle, na.rm = TRUE), k_max = max(k_cycle, na.rm = TRUE)), by = variant][order(variant)])
cat("\nk_cycle distribution:\n"); print(table(o$k_cycle, useNA = "ifany"))
