## concat_catalog_v6.R (October 30, 2026): catalog v6 long file and seed gate = runs_v5 default rows (unchanged by
## the PR 135 overlay fix) + runs_v6 calibrated rows (rerun with the fixed code). Thresholds files are reused as is.
suppressMessages(library(data.table))
f5 <- list.files("runs_v5", "^catalog_[a-z]+_(native|gompit)\\.csv$", full.names = TRUE)
f6 <- list.files("runs_v6", "^catalog_[a-z]+_(native|gompit)\\.csv$", full.names = TRUE)
cat("files runs_v5", length(f5), "runs_v6", length(f6), "\n")
d5 <- rbindlist(lapply(f5, fread), fill = TRUE)[config == "default"]
d6 <- rbindlist(lapply(f6, fread), fill = TRUE)
stopifnot(all(d6$config == "calibrated"), !any(as.logical(d6$partial)), !any(as.logical(d5$partial)))
d <- rbind(d5, d6, fill = TRUE); d[, traj := paste(group, traj, sep = ":")]
cat("rows default", nrow(d5), "calibrated", nrow(d6), "arms", paste(sort(unique(d$arm)), collapse = ","), "groups", uniqueN(d$group), "\n")
fwrite(d, "catalog_long_v6.csv")
g5 <- rbindlist(lapply(list.files("runs_v5", "_seedgate\\.csv$", full.names = TRUE), fread), fill = TRUE)
g6 <- rbindlist(lapply(list.files("runs_v6", "_seedgate\\.csv$", full.names = TRUE), fread), fill = TRUE)
sg <- rbind(g5[config == "default"], g6, fill = TRUE)
fwrite(sg, "catalog_seedgate_v6.csv"); print(sg[, .N, by = .(status)])
print(sg[status %in% c("OK", "SEED_FAIL"), .(ok = mean(status == "OK"), n = .N), by = .(variant, arm, config)][order(variant, arm, config)])
c5 <- g5[config == "calibrated", .(variant, arm, plt_cn, s5 = status)]; c6 <- g6[config == "calibrated", .(variant, arm, plt_cn, s6 = status)]
m <- merge(c5, c6, by = c("variant", "arm", "plt_cn"), all = TRUE)
cat("calibrated seed-gate rows v5", nrow(c5), "v6", nrow(c6), "status changed", m[is.na(s5) | is.na(s6) | s5 != s6, .N], "\n")
## no plot identifiers in the copy that leaves Cardinal
fwrite(sg[, !"plt_cn"], "catalog_seedgate_v6_noid.csv")
