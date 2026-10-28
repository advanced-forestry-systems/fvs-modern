#!/usr/bin/env Rscript
## Observed BA and SDI percentiles by FVS variant from the tree-level CONUS remeasurement pairs RDS
## (one row per tree; plot-level BA1/TPA1/SDI1/QMD1 repeated per tree, so dedupe on plot_key).
## Writes only percentiles; no coordinates or plot identifiers leave this script.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); x <- as.data.table(readRDS(a[1]))
p <- unique(x[, .(plot_key, variant = toupper(fvs_variant), BA = BA1, TPA = TPA1, QMD = QMD1, SDI = SDI1, RD = RD1)], by = "plot_key")
cat(sprintf("plots %d; BA median %.1f max %.1f; TPA median %.0f; QMD median %.1f; SDI median %.0f\n", nrow(p), median(p$BA,na.rm=T), max(p$BA,na.rm=T), median(p$TPA,na.rm=T), median(p$QMD,na.rm=T), median(p$SDI,na.rm=T)))
## unit check: the file is labelled metric; if BA median > 60 it is ft2/ac and TPA is per acre
if (median(p$BA, na.rm = TRUE) > 60) { cat("imperial detected, converting\n"); p[, BA := BA * 0.2295684][, TPA := TPA * 2.4710538][, QMD := QMD * 2.54] }
p[, SDI_m := ifelse(TPA > 0 & QMD > 0, TPA * (QMD / 25.4)^1.605, NA_real_)]
out <- p[is.finite(BA) & is.finite(SDI_m), .(n = .N, ba_p99 = quantile(BA, .99), ba_p999 = quantile(BA, .999), ba_max = max(BA),
          sdi_p99 = quantile(SDI_m, .99), sdi_p999 = quantile(SDI_m, .999), sdi_max = max(SDI_m), sdi_file_p99 = quantile(SDI, .99, na.rm = TRUE)), by = variant][order(variant)]
fwrite(out, a[2]); print(out, digits = 4)
cat("Units: BA m2 ha-1; SDI_m trees ha-1 at 25.4 cm computed from TPA and QMD; sdi_file_p99 is the file's own SDI column for comparison.\n")
