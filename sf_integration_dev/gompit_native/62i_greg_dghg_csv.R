#!/usr/bin/env Rscript
# 62i_greg_dghg_csv.R
# Emit Greg's deployed DG and HG coefficient tables as flat per-species CSVs in
# the same GOMPLOAD-style format the native hooks will read (header + SPCD n
# B0.. ). DG uses B0..B6 (7 params); HG uses B0..B8 (9 params, B0 = max_height).
# De-duplicated to one best-fit row per species (isConv, then min RSS).
suppressPackageStartupMessages(library(data.table))
RDS <- "/users/PUOM0008/crsfaaron/fvs_remodeling/rds"
OUT <- "/users/PUOM0008/crsfaaron/fvs-modern/config"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dedup <- function(dt) {
  dt <- copy(dt); setnames(dt, names(dt), names(dt))
  if ("isConv" %in% names(dt)) dt <- dt[isConv == TRUE | is.na(isConv)]
  ord <- c("spcd", intersect("RSS", names(dt)))
  setorderv(dt, ord)
  dt[, .SD[1], by = spcd]
}

dg <- dedup(as.data.table(readRDS(file.path(RDS, "dg_parms.RDS"))))
dgo <- dg[, .(SPCD = as.integer(spcd), n = as.integer(n),
              B0, B1, B2, B3, B4, B5, B6)][order(SPCD)]
dgo <- dgo[is.finite(B0) & is.finite(B1) & is.finite(B2) & is.finite(B3) &
           is.finite(B4) & is.finite(B5) & is.finite(B6)]
fwrite(dgo, file.path(OUT, "greg_dg_coefficients.csv"), quote = FALSE)

hg <- dedup(as.data.table(readRDS(file.path(RDS, "hg_parms.RDS"))))
hgo <- hg[, .(SPCD = as.integer(spcd), n = as.integer(n),
              B0, B1, B2, B3, B4, B5, B6, B7, B8)][order(SPCD)]
hgo <- hgo[is.finite(B0) & is.finite(B1) & is.finite(B2) & is.finite(B3) &
           is.finite(B4) & is.finite(B5) & is.finite(B6) & is.finite(B7) & is.finite(B8)]
fwrite(hgo, file.path(OUT, "greg_hg_coefficients.csv"), quote = FALSE)

cat(sprintf("DG: %d species -> %s/greg_dg_coefficients.csv\n", nrow(dgo), OUT))
cat(sprintf("HG: %d species -> %s/greg_hg_coefficients.csv\n", nrow(hgo), OUT))
cat("DG cols: SPCD,n,B0..B6 ; HG cols: SPCD,n,B0..B8 (B0=max_height)\n")
