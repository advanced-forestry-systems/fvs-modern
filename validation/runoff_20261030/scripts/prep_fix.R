## prep_fix.R (Amendment A3): inputs for the primary (fixed calibrated) and sensitivity (_inproc) rescoring.
suppressMessages(library(data.table))
F <- "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full"; C <- "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog"
dir.create(file.path(F, "score_fix/inputs"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(F, "score_fix_sens/inputs"), recursive = TRUE, showWarnings = FALSE)
old <- Sys.glob(file.path(F, "a_near/runoff_near_*.csv")); old <- old[!grepl("_gate\\.csv$", old)]
A0 <- rbindlist(lapply(old, fread), fill = TRUE)
fwrite(A0[config == "default"], file.path(F, "score_fix/inputs/a_near_default.csv"))
inp <- A0[config == "calibrated"][, arm := paste0(arm, "_inproc")]; fwrite(inp, file.path(F, "score_fix_sens/inputs/a_near_calibrated_inproc.csv"))
fx <- Sys.glob(file.path(F, "engines_fix/a_near/runoff_near_*.csv")); fx <- fx[!grepl("_gate\\.csv$", fx)]
AF <- rbindlist(lapply(fx, fread), fill = TRUE); fwrite(AF, file.path(F, "score_fix/inputs/a_near_calibrated_fix.csv"))
cat("near rows: default", nrow(A0[config == "default"]), " calibrated fixed", nrow(AF), " calibrated inproc", nrow(inp), "\n")
keep <- c("arm", "group", "site", "density", "traj", "age", "year", "TPH", "BAPH", "QMD", "HT", "VOL", "VOL_NATIVE", "seed_ok")
rdcat <- function(f) { d <- fread(f); setnames(d, "VOL", "VOL_NATIVE"); d[, VOL := 0.40 * BAPH * HT]; d[, ..keep] }
A6 <- rdcat(file.path(C, "catalog_long_v6.csv")); A5c <- rdcat(file.path(C, "catalog_long_v5.csv"))[grepl("^calibrated", arm)][, arm := paste0(arm, "_inproc")]
cf <- Sys.glob(file.path(F, "c_long/runoff_long_*_native.csv")); CL <- rbindlist(lapply(cf, fread), fill = TRUE); CL[, traj := paste(group, traj, sep = ":")]
P <- rbind(A6, CL[, ..keep], fill = TRUE); fwrite(P, file.path(F, "score_fix/eng_long.csv")); fwrite(rbind(P, A5c, fill = TRUE), file.path(F, "score_fix_sens/eng_long.csv"))
print(P[, .(stands = uniqueN(traj)), by = arm]); cat("v5 calibrated inproc stands:", uniqueN(A5c$traj), "\n")
