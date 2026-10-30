## d2_compare_v6.R: inproc vs kwonly on the fixed code (catalog d2check_v6), and the default-arm reuse check.
suppressMessages(library(data.table))
D <- "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/d2check_v6"; C <- dirname(D)
rd <- function(m) { d <- rbindlist(lapply(Sys.glob(file.path(D, m, "catalog_*_native.csv")), fread), fill = TRUE); d[, SDI := TPH * (QMD / 25.4)^1.605]; d }
a <- rd("inproc"); b <- rd("kwonly")
m <- merge(a[, .(group, traj, age, BA_i = BAPH, TPH_i = TPH, SDI_i = SDI)], b[, .(group, traj, age, BA_k = BAPH, TPH_k = TPH, SDI_k = SDI)], by = c("group", "traj", "age"), all = TRUE)
m[, rel50 := ifelse(age == 50, abs(BA_i - BA_k) / BA_k, NA_real_)]
r <- m[, .(stands = uniqueN(traj), rows = .N, rows_missing_one_side = sum(!is.finite(BA_i) | !is.finite(BA_k)),
  identical_rows_share = round(mean(abs(BA_i - BA_k) < 1e-9 & abs(TPH_i - TPH_k) < 1e-9, na.rm = TRUE), 4),
  max_abs_dBA = signif(max(abs(BA_i - BA_k), na.rm = TRUE), 4), max_abs_dTPH = signif(max(abs(TPH_i - TPH_k), na.rm = TRUE), 4),
  max_abs_dSDI = signif(max(abs(SDI_i - SDI_k), na.rm = TRUE), 4), max_rel_dBA_yr50 = signif(max(rel50, na.rm = TRUE), 4),
  stands_yr50_dBA_over_1pct = sum(rel50 > 0.01, na.rm = TRUE)), by = group][order(group)]
s <- rbindlist(lapply(Sys.glob(file.path(D, "*", "spsdi_*.csv")), function(f) { d <- fread(f); d[, group := toupper(sub("spsdi_(.*)\\.csv", "\\1", basename(f)))]; d }), fill = TRUE)
z <- s[, .(runs = .N, n_species = max(n_species), max_spsdi_zero = max(n_spsdi_zero)), by = .(group, mode)]
r <- merge(r, dcast(z, group ~ mode, value.var = c("max_spsdi_zero", "n_species")), by = "group", all.x = TRUE)
options(width = 250); print(r); fwrite(r, file.path(D, "d2_compare_by_variant_v6.csv"))
## default arm: fixed code vs runs_v5
n <- rbindlist(lapply(Sys.glob(file.path(D, "default", "catalog_*_native.csv")), fread), fill = TRUE)
o <- rbindlist(lapply(c("ca", "wc"), function(v) fread(file.path(C, "runs_v5", sprintf("catalog_%s_native.csv", v)))), fill = TRUE)[config == "default" & traj %in% n$traj]
k <- c("TPH", "BAPH", "QMD", "HT", "VOL")
dd <- merge(n[, c("group", "traj", "age", k), with = FALSE], o[, c("group", "traj", "age", k), with = FALSE], by = c("group", "traj", "age"), suffixes = c("_fix", "_v5"), all = TRUE)
dr <- dd[, .(stands = uniqueN(traj), rows = .N, missing = sum(!complete.cases(.SD)),
             max_abs_diff = max(sapply(k, function(x) max(abs(get(paste0(x, "_fix")) - get(paste0(x, "_v5"))), na.rm = TRUE)))), by = group]
print(dr); fwrite(dr, file.path(D, "default_reuse_check_v6.csv"))
