## Observed FIA heights (live, measured, unbroken) by species vs Greg's HG asymptote B0, 2026 FIADB build.
suppressPackageStartupMessages(library(data.table))
D <- "/users/PUOM0008/crsfaaron/fvs_remodeling_data/2026-11-04"; O <- "/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/ooda_20261110"
Z <- file.path(D, "TREEdata-11-04-2026.zip")
fs <- grep("_TREEdata.CSV$", unzip(Z, list = TRUE)$Name, value = TRUE)
NE <- c("ME","NH","VT","NY","MA","CT","RI","PA","NJ")
sp <- c(318,371,531,541,833,129,316,802,12,97,95,261)
rd <- function(f) fread(cmd = sprintf("unzip -p %s %s", Z, f),
  select = c("SPCD","STATUSCD","HT","HTCD","ACTUALHT","DIA"), showProgress = FALSE)[STATUSCD == 1 & HTCD == 1 & !is.na(HT) & HT == ACTUALHT][, st := substr(f, 1, 2)]
x <- rbindlist(lapply(fs, rd))
q <- function(v, p) as.numeric(quantile(v, p))
s <- x[, .(n = .N, p95 = q(HT,.95), p99 = q(HT,.99), p999 = q(HT,.999), max = max(HT)), by = .(SPCD, region = fifelse(st %in% NE, "NE", "other"))]
s <- rbind(s, x[, .(n = .N, p95 = q(HT,.95), p99 = q(HT,.99), p999 = q(HT,.999), max = max(HT), region = "CONUS"), by = SPCD])
hg <- fread("/users/PUOM0008/crsfaaron/fvs-modern-wt-greg2026/config/greg_hg_coefficients.csv")[, .(SPCD, B0)]
out <- merge(s[SPCD %in% sp], hg, by = "SPCD")[order(SPCD, region)]
fwrite(out, file.path(O, "fia_obs_height_vs_B0.csv")); print(out[region != "other"])
cat("rows used", nrow(x), "\n")
