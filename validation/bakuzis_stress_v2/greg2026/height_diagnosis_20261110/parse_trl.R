suppressPackageStartupMessages(library(data.table))
O <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/bakuzis_greg2026_20261109/treelist_20261110"
P <- "/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/ooda_20261110"
grp <- function(n) c(rep("Spruce-Fir",9),rep("Northern-Hardwood",9),rep("Pine",9),rep("Oak-Pine",9))[n]
cell <- function(n) { k <- (n-1) %% 9; paste(c("Low","Medium","High")[k %/% 3 + 1], c("Low","Medium","High")[k %% 3 + 1], sep="/") }
rd <- function(d) {
  L <- readLines(file.path(d, "bk.trl")); yr <- NA_integer_; out <- vector("list", length(L)); k <- 0L
  for (l in L) {
    m <- regmatches(l, regexec("YEAR:\\s*([0-9]{4})", l))[[1]]
    if (length(m)) { yr <- as.integer(m[2]); next }
    t <- strsplit(trimws(l), "\\s+")[[1]]
    if (length(t) >= 18 && grepl("^[0-9]+$", t[1]) && grepl("^[A-Z]{2}$", t[3])) { k <- k + 1L
      out[[k]] <- list(year = yr, sp = t[3], tpa = as.numeric(t[8]), dbh = as.numeric(t[10]), ht = as.numeric(t[12])) }
  }
  x <- rbindlist(out[seq_len(k)]); tag <- sub("run_([a-z]+)_BK.*", "\\1", basename(d)); n <- as.integer(sub(".*BK", "", basename(d)))
  x[, `:=`(arm = tag, bk = n, group = grp(n), cell = cell(n))]
}
x <- rbindlist(lapply(list.dirs(O, recursive = FALSE), rd)); x[, yr := year - 2000L]
top <- x[tpa > 0][order(arm, bk, yr, -dbh)][, cum := cumsum(tpa), by = .(arm, bk, yr)][cum - tpa < 40]
top[, w := pmin(tpa, 40 - (cum - tpa))]
s1 <- top[, .(topht = sum(w*ht)/sum(w)), by = .(arm, group, cell, yr)]
s2 <- top[, .(share = sum(w)/40, ht = sum(w*ht)/sum(w)), by = .(arm, group, cell, yr, sp)]
fwrite(x, file.path(P, "treelist_long.csv.gz")); fwrite(s1, file.path(P, "topht_by_arm.csv")); fwrite(s2, file.path(P, "topht_species_share.csv"))
print(dcast(s1[cell == "Medium/Medium" & yr %in% c(0,50,100)], group + yr ~ arm, value.var = "topht"))
print(s2[cell == "Medium/Medium" & yr == 100][order(group, arm, -share)][, .(group, arm, sp, share = round(share, 2), ht = round(ht))], nrows = 80)
