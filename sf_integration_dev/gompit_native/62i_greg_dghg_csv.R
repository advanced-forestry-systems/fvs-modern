#!/usr/bin/env Rscript
# 62i_greg_dghg_csv.R
# Emit Greg's deployed DG, HG and crown recession coefficient tables as flat
# per-species CSVs in the GOMPLOAD-style format the native hooks read
# (header + SPCD n B0.. ). DG uses B0..B6 (PC climate form); HG uses B0..B7
# (B0 = max_height, PC climate form); crown recession uses b0..b2.
# De-duplicated to one best-fit row per species (isConv, then min RSS or rmse).
#   Rscript 62i_greg_dghg_csv.R --rds=<fvs_remodeling/rds> --out=<config dir>
suppressPackageStartupMessages(library(data.table))
ga <- function(n,d){ m<-grep(paste0("^--",n,"="),commandArgs(TRUE),value=TRUE); if(length(m)) sub(paste0("^--",n,"="),"",m[1]) else d }
RDS <- ga("rds", "/users/PUOM0008/crsfaaron/fvs_remodeling/rds")
OUT <- ga("out", "/users/PUOM0008/crsfaaron/fvs-modern/config")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dedup <- function(dt, key = "spcd") {
  dt <- copy(dt)
  if ("isConv" %in% names(dt)) dt <- dt[isConv == TRUE | is.na(isConv)]
  ord <- c(key, intersect(c("RSS", "rmse"), names(dt)))
  setorderv(dt, ord)
  dt[, .SD[1], by = key]
}
fin <- function(dt, cols) dt[Reduce(`&`, lapply(cols, function(v) is.finite(dt[[v]])))]

dg <- dedup(as.data.table(readRDS(file.path(RDS, "dg_parms.RDS"))))
dgo <- fin(dg[, .(SPCD = as.integer(spcd), n = as.integer(n),
                  B0, B1, B2, B3, B4, B5, B6)][order(SPCD)], paste0("B", 0:6))
fwrite(dgo, file.path(OUT, "greg_dg_coefficients.csv"), quote = FALSE)

hg <- dedup(as.data.table(readRDS(file.path(RDS, "hg_parms.RDS"))))
hgo <- fin(hg[, .(SPCD = as.integer(spcd), n = as.integer(n),
                  B0, B1, B2, B3, B4, B5, B6, B7)][order(SPCD)], paste0("B", 0:7))
fwrite(hgo, file.path(OUT, "greg_hg_coefficients.csv"), quote = FALSE)

cr <- as.data.table(readRDS(file.path(RDS, "delta_height_to_live_crown.RDS")))
setnames(cr, "SPCD", "spcd", skip_absent = TRUE)
cr <- dedup(cr)
cro <- fin(cr[, .(SPCD = as.integer(spcd), n = as.integer(n), b0, b1, b2)][order(SPCD)], paste0("b", 0:2))
fwrite(cro, file.path(OUT, "greg_crown_change_coefficients.csv"), quote = FALSE)

cat(sprintf("DG: %d species -> %s/greg_dg_coefficients.csv\n", nrow(dgo), OUT))
cat(sprintf("HG: %d species -> %s/greg_hg_coefficients.csv\n", nrow(hgo), OUT))
cat(sprintf("CRW: %d species -> %s/greg_crown_change_coefficients.csv\n", nrow(cro), OUT))
