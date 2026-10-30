#!/usr/bin/env Rscript
## runoff_long_prep.R OUT GROUPS FIRST CDIR
## Engine long table for the run-off scorer: contestant A from catalog/catalog_long_v5.csv (read only,
## never modified) and contestant C from runoff_long_<v>_native.csv files in CDIR. VOL is set to the common
## form-factor volume 0.40 x BAPH x HT for every contestant (the projector's VOL in the scorer is the same
## formula); the engines' native tcuft volume is kept as VOL_NATIVE. FIRST > 0 keeps traj C0001..C<FIRST>.
suppressMessages(library(data.table))
a <- commandArgs(TRUE); OUT <- a[1]; G <- toupper(strsplit(a[2], ",")[[1]]); FIRST <- as.integer(a[3]); CDIR <- a[4]
keep <- c("arm", "group", "site", "density", "traj", "age", "year", "TPH", "BAPH", "QMD", "HT", "VOL", "VOL_NATIVE", "seed_ok")
A <- fread("/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/catalog_long_v5.csv")
A <- A[toupper(group) %in% G]; setnames(A, "VOL", "VOL_NATIVE")
if (FIRST > 0) A <- A[as.integer(sub("^.*:C", "", traj)) <= FIRST]
A[, VOL := 0.40 * BAPH * HT]
cf <- list.files(CDIR, "^runoff_long_[a-z]+_native\\.csv$", full.names = TRUE)
C <- if (length(cf)) rbindlist(lapply(cf, fread), fill = TRUE) else data.table()
if (nrow(C)) { C <- C[toupper(group) %in% G]; C[, traj := paste(group, traj, sep = ":")] }
L <- rbind(A[, ..keep], if (nrow(C)) C[, ..keep] else NULL, fill = TRUE)
fwrite(L, OUT)
print(L[, .(stands = uniqueN(traj), rows = .N), by = .(arm, group)][order(group, arm)])
