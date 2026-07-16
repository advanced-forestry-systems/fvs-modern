#!/usr/bin/env Rscript
# build_emt_td_crosswalk.R  (P11)
# Build a plot_key -> EMT/TD crosswalk so Greg's climate DG (needs EMT) can be
# scored on our CONUS remeasurement pairs, and the climate DG variant becomes
# CONUS-wide. EMT/TD are 1991-2020 ClimateNA normals (per LOCATION, not per
# measurement), so any of a plot's FIA CNs carries the same value; we collapse
# to one EMT/TD per plot_key.
#
# Inputs:
#   ~/fia_data/{STATE}_PLOT.csv        FIA PLOT tables (CN, STATECD, UNITCD, COUNTYCD, PLOT)
#   ~/fvs-modern/config/greg_emt_td_lookup.csv   STAND_CN(=CN), EMT, TD  (1.84M rows)
#   ~/fvs-conus/calibration/data/conus_remeasurement_pairs.rds   plot_key
# Output:
#   .../plot_key_emt_td_crosswalk.csv  + coverage report
# Author: Cowork autopilot 2026-07-03. CN read as character (15-digit).

suppressPackageStartupMessages(library(data.table))
OUT <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/emt_crosswalk"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. FIA PLOT tables -> plot_key -> CN
pfiles <- list.files("~/fia_data", pattern = "^[A-Z]{2}_PLOT\\.csv$",
                     full.names = TRUE)
cat("PLOT files:", length(pfiles), "\n")
plots <- rbindlist(lapply(pfiles, function(f)
  tryCatch(fread(f, select = c("CN","STATECD","UNITCD","COUNTYCD","PLOT"),
                 colClasses = list(character = "CN")), error = function(e) NULL)),
  fill = TRUE)
plots <- plots[!is.na(CN)]
plots[, plot_key := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "-")]
cat("plot rows:", format(nrow(plots), big.mark=","),
    " distinct plot_key:", format(uniqueN(plots$plot_key), big.mark=","), "\n")

# 2. join EMT/TD lookup by CN
lk <- fread("~/fvs-modern/config/greg_emt_td_lookup.csv",
            colClasses = list(character = "STAND_CN"))
setnames(lk, "STAND_CN", "CN")
pj <- merge(plots, lk, by = "CN")
cat("plot rows with EMT:", format(nrow(pj), big.mark=","), "\n")

# 3. collapse to one EMT/TD per plot_key (normals -> mean is exact)
xw <- pj[, .(EMT = mean(EMT, na.rm=TRUE), TD = mean(TD, na.rm=TRUE),
             n_cn = .N), by = plot_key]
fwrite(xw, file.path(OUT, "plot_key_emt_td_crosswalk.csv"))
cat("crosswalk plot_keys:", format(nrow(xw), big.mark=","), "\n")

# 4. coverage against the remeasurement pairs
d <- as.data.table(readRDS("~/fvs-conus/calibration/data/conus_remeasurement_pairs.rds"))
d <- merge(d, xw[, .(plot_key, EMT, TD)], by = "plot_key", all.x = TRUE)
cov <- mean(is.finite(d$EMT))
cat(sprintf("\n=== COVERAGE: %.1f%% of %s remeasurement pairs now have EMT/TD ===\n",
            100*cov, format(nrow(d), big.mark=",")))
cat(sprintf("EMT range [%.1f, %.1f] C; TD range [%.1f, %.1f]\n",
            min(d$EMT,na.rm=TRUE), max(d$EMT,na.rm=TRUE),
            min(d$TD,na.rm=TRUE), max(d$TD,na.rm=TRUE)))
# also emit a pairs-level EMT/TD table for downstream A/B (plot_key,EMT,TD)
fwrite(unique(d[is.finite(EMT), .(plot_key, EMT, TD)]),
       file.path(OUT, "pairs_plot_key_emt_td.csv"))
cat("wrote crosswalk + pairs EMT table to", OUT, "\n")
