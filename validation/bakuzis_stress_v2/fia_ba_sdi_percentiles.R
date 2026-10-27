#!/usr/bin/env Rscript
## fia_ba_sdi_percentiles.R  (2026-10-27)
## Observed basal area and SDI percentiles by FVS variant from the CONUS remeasurement pairs, for the
## magnitude test in bakuzis_stress_v2.R (--obs_ba_p999 and --sdimax). Column-name driven: it stops
## with the column list if it cannot find basal area, stems and a variant column, rather than guessing.
## Usage: Rscript fia_ba_sdi_percentiles.R <pairs.rds> <out.csv>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: fia_ba_sdi_percentiles.R <pairs.rds> <out.csv>")
x <- readRDS(args[1]); x <- as.data.table(x); nm <- names(x)
pick <- function(c) { h <- c[c %in% nm]; if (length(h)) h[1] else NA_character_ }
cBA  <- pick(c("BA_M2HA_T1","BAPH_T1","BA_T1","BAPH","BA_M2HA","BA_t1","ba_m2ha_t1","BALIVE_T1"))
cTPH <- pick(c("TPH_T1","TPH","TPA_T1","STEMS_HA_T1","tph_t1","TPA"))
cVAR <- pick(c("VARIANT","variant","FVS_VARIANT","fvs_variant"))
cQMD <- pick(c("QMD_CM_T1","QMD_T1","QMD","qmd_cm_t1"))
if (any(is.na(c(cBA, cTPH, cVAR)))) {
  cat("Could not identify columns. Names present:\n"); print(nm); stop("edit the pick() lists for this file") }
cat(sprintf("using BA=%s stems=%s variant=%s qmd=%s\n", cBA, cTPH, cVAR, ifelse(is.na(cQMD), "(derived)", cQMD)))
d <- x[, .(variant = get(cVAR), BA = as.numeric(get(cBA)), N = as.numeric(get(cTPH)))]
## unit sanity: BA in m2/ha should be < 200 and stems/ha < 50000; if BA looks like ft2/ac and N like TPA, convert
if (median(d$BA, na.rm = TRUE) > 60 && grepl("TPA", cTPH)) { cat("converting imperial to metric\n"); d[, BA := BA * 0.2295684][, N := N * 2.4710538] }
d[, QMD := ifelse(N > 0, sqrt(BA / N * 40000 / pi), NA_real_)]
d[, SDI := ifelse(N > 0 & QMD > 0, N * (QMD / 25.4)^1.605, NA_real_)]
out <- d[is.finite(BA) & is.finite(SDI), .(n = .N, ba_p99 = quantile(BA, 0.99), ba_p999 = quantile(BA, 0.999), ba_max = max(BA),
                                          sdi_p99 = quantile(SDI, 0.99), sdi_p999 = quantile(SDI, 0.999), sdi_max = max(SDI)), by = variant][order(variant)]
fwrite(out, args[2]); print(out)
cat("Units: BA m2 ha-1, SDI trees ha-1 at 25.4 cm. Pass sdi_p999 as --sdimax (group,sdimax) and ba_p999 as --obs_ba_p999.\n")
