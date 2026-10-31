## concat_catalog_v8.R (October 30, 2026): catalog v8 = all four arms rerun with D5 (runs_v8), long file and seed gate.
suppressMessages(library(data.table))
f8 <- list.files("runs_v8", "^catalog_[a-z]+_(native|gompit)\\.csv$", full.names = TRUE); cat("files runs_v8", length(f8), "\n")
d <- rbindlist(lapply(f8, fread), fill = TRUE); stopifnot(!any(as.logical(d$partial))); d[, traj := paste(group, traj, sep = ":")]
cat("rows", nrow(d), "arms", paste(sort(unique(d$arm)), collapse = ","), "groups", uniqueN(d$group), "stands per arm", d[, uniqueN(traj), by = arm][, paste(V1, collapse = ",")], "\n")
fwrite(d, "catalog_long_v8.csv")
sg <- rbindlist(lapply(list.files("runs_v8", "_seedgate\\.csv$", full.names = TRUE), fread), fill = TRUE)
fwrite(sg, "catalog_seedgate_v8.csv"); print(sg[, .N, by = .(status)])
sg[, r := ba0_ft2ac / balive]; cat("year-0 BA / BALIVE (v8): median", round(median(sg$r, na.rm = TRUE), 4), " within 2%:", round(mean(abs(sg$r - 1) <= 0.02, na.rm = TRUE), 4), "\n")
fwrite(sg[, !"plt_cn"], "catalog_seedgate_v8_noid.csv")
