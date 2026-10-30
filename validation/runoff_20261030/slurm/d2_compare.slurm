#!/bin/bash
#SBATCH --job-name=ro_d2cmp --account=PUOM0008 --time=00:20:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/d2check/logs/d2cmp_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
D=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/d2check
for f in $D/logs/d2_*_*.out; do echo "$(grep '>>>' $f | cut -d' ' -f3-4) | $(grep -E 'status counts' $f | sed 's/.*status counts: //') $(grep -c -E 'Traceback|failed' $f)err"; done
Rscript -e '
suppressMessages(library(data.table)); D <- "'$D'"
rd <- function(m) { d <- rbindlist(lapply(Sys.glob(file.path(D, m, "runoff_long_*_native.csv")), fread), fill = TRUE); d[, SDI := TPH * (QMD / 25.4)^1.605]; d }
a <- rd("inproc"); b <- rd("kwonly")
m <- merge(a[, .(group, traj, age, BA_i = BAPH, TPH_i = TPH, SDI_i = SDI)], b[, .(group, traj, age, BA_k = BAPH, TPH_k = TPH, SDI_k = SDI)], by = c("group", "traj", "age"), all = TRUE)
m[, rel50 := ifelse(age == 50, abs(BA_i - BA_k) / BA_k, NA_real_)]
r <- m[, .(stands = uniqueN(traj), rows = .N, rows_missing_one_side = sum(!is.finite(BA_i) | !is.finite(BA_k)),
  identical_rows_share = round(mean(abs(BA_i - BA_k) < 1e-9 & abs(TPH_i - TPH_k) < 1e-9, na.rm = TRUE), 4),
  max_abs_dBA = signif(max(abs(BA_i - BA_k), na.rm = TRUE), 4), max_abs_dTPH = signif(max(abs(TPH_i - TPH_k), na.rm = TRUE), 4),
  max_abs_dSDI = signif(max(abs(SDI_i - SDI_k), na.rm = TRUE), 4), max_rel_dBA_yr50 = signif(max(rel50, na.rm = TRUE), 4),
  stands_yr50_dBA_over_1pct = sum(rel50 > 0.01, na.rm = TRUE)), by = group][order(group)]
s <- rbindlist(lapply(Sys.glob(file.path(D, "*", "spsdi_*.csv")), function(f) { d <- fread(f); d[, group := toupper(sub("spsdi_(.*)\\.csv", "\\1", basename(f)))]; d }), fill = TRUE)
z <- s[, .(runs = .N, n_species = max(n_species), max_spsdi_zero = max(n_spsdi_zero), min_spsdi_zero = min(n_spsdi_zero)), by = .(group, mode)]
zz <- dcast(z, group ~ mode, value.var = c("max_spsdi_zero", "n_species"))
r <- merge(r, zz, by = "group", all.x = TRUE); print(r); fwrite(r, file.path(D, "d2_compare_by_variant.csv"))
cat("MATERIAL (any yr-50 BA diff > 1 pct):", paste(r[stands_yr50_dBA_over_1pct > 0, group], collapse = ","), "\n")'
md5sum $D/d2_compare_by_variant.csv $D/d2_wrapper.py
echo D2CMP_DONE
