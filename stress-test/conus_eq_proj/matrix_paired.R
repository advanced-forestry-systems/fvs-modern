#!/usr/bin/env Rscript
# Stand-paired growth x mortality 2x2 for the CONUS equation arms.
# Fixes the population confound in growth_mort_matrix_2100.csv: b1_gompit pooled only
# n=185,658 stands (array tasks 9/11/13 OOM'd, 2 variants missing) vs n=741,392 for the
# other three arms, so the unpaired pooled means mixed the mortality effect with composition.
# Here the mortality effect is estimated WITHIN stand, on the set of stands present in ALL
# four arms (b1/b2 growth x native/gompit mortality). Two-pass streaming read, memory bounded.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260629)
D   <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
OUT <- D
elog <- file.path(OUT, "error_log.txt")
logerr <- function(where, e) cat(sprintf("[%s] %s: %s\n", Sys.time(), where, conditionMessage(e)),
                                 file = elog, append = TRUE)

arms <- list(
  b1_native = file.path(D, "out_conus_eq",          "*conus_b1_metrics.csv"),
  b2_native = file.path(D, "out_conus_eq",          "*conus_b2_metrics.csv"),
  b2_gompit = file.path(D, "out_conus_eq_gompit",   "*conus_b2_gompit_metrics.csv"),
  b1_gompit = file.path(D, "out_conus_eq_b1_gompit","*conus_b1_gompit_metrics.csv")
)
sel <- c("STAND_CN","STATE","PROJ_YEAR","BA_FT2AC","QMD_IN","TPH","AGB_TONS_AC")
add_sdi <- function(d) {
  d <- d[is.finite(BA_FT2AC) & is.finite(QMD_IN) & is.finite(TPH) & QMD_IN > 0]
  d[, TPA := TPH / 2.471][, SDI := TPA * (QMD_IN / 10)^1.605][]
}

## ---- PASS 1: terminal-year row per stand, per arm ----
terminal <- list()
for (a in names(arms)) {
  fs <- Sys.glob(arms[[a]])
  rows <- rbindlist(lapply(fs, function(f) tryCatch({
    d <- fread(f, select = sel); d <- add_sdi(d)
    d[d[, .I[PROJ_YEAR == max(PROJ_YEAR)], by = STAND_CN]$V1]
  }, error = function(e) { logerr(paste("pass1", f), e); NULL })), fill = TRUE)
  rows <- rows[, .(STATE = STATE[1], BA = BA_FT2AC, SDI = SDI, AGB = AGB_TONS_AC), by = STAND_CN]
  terminal[[a]] <- rows
  cat(sprintf("pass1 %s: %d stands\n", a, nrow(rows)))
  gc()
}
## inner join across all four arms on STAND_CN -> common stand set
common <- Reduce(function(x, y) merge(x, y, by = "STAND_CN"),
                 lapply(names(arms), function(a) {
                   z <- terminal[[a]][, .(STAND_CN, STATE,
                                          BA = BA, SDI = SDI, AGB = AGB)]
                   setnames(z, c("BA","SDI","AGB"), paste0(c("BA_","SDI_","AGB_"), a)); z
                 }))
common <- unique(common, by = "STAND_CN")
n_common <- nrow(common)
cat(sprintf("common stands present in all 4 arms: %d\n", n_common))
fwrite(common, file.path(OUT, "matrix_paired_terminal_wide.csv"))

## ---- 2x2 arm means on the COMMON set ----
arm_mean <- rbindlist(lapply(names(arms), function(a) {
  data.table(arm = a,
             BA_2100  = round(mean(common[[paste0("BA_",  a)]], na.rm = TRUE), 1),
             SDI_2100 = round(mean(common[[paste0("SDI_", a)]], na.rm = TRUE)),
             AGB_2100 = round(mean(common[[paste0("AGB_", a)]], na.rm = TRUE), 1),
             n = n_common)
}))
fwrite(arm_mean, file.path(OUT, "matrix_paired_2100.csv"))
print(arm_mean)

## ---- within-stand paired mortality effect (native - gompit) ----
paired <- function(g) {
  dBA  <- common[[paste0("BA_",  g, "_native")]] - common[[paste0("BA_",  g, "_gompit")]]
  dSDI <- common[[paste0("SDI_", g, "_native")]] - common[[paste0("SDI_", g, "_gompit")]]
  data.table(growth_eq = g,
             dBA_mean = round(mean(dBA, na.rm = TRUE), 1),
             dBA_sd   = round(sd(dBA,  na.rm = TRUE), 1),
             dBA_median = round(median(dBA, na.rm = TRUE), 1),
             BA_pct_reduction = round(100 * mean(dBA, na.rm = TRUE) /
                                      mean(common[[paste0("BA_", g, "_native")]], na.rm = TRUE), 1),
             dSDI_mean = round(mean(dSDI, na.rm = TRUE)),
             SDI_pct_reduction = round(100 * mean(dSDI, na.rm = TRUE) /
                                       mean(common[[paste0("SDI_", g, "_native")]], na.rm = TRUE), 1))
}
# rename terminal cols to the g_mortality scheme expected above
setnames(common,
         old = c("BA_b1_native","SDI_b1_native","BA_b1_gompit","SDI_b1_gompit",
                 "BA_b2_native","SDI_b2_native","BA_b2_gompit","SDI_b2_gompit"),
         new = c("BA_b1_native","SDI_b1_native","BA_b1_gompit","SDI_b1_gompit",
                 "BA_b2_native","SDI_b2_native","BA_b2_gompit","SDI_b2_gompit"),
         skip_absent = TRUE)
paired_tab <- rbindlist(lapply(c("b1","b2"), paired))
fwrite(paired_tab, file.path(OUT, "matrix_paired_effect.csv"))
print(paired_tab)
keep_cn <- common$STAND_CN
gc()

## ---- PASS 2: mean BA trajectory on the COMMON stands only ----
traj <- rbindlist(lapply(names(arms), function(a) {
  fs <- Sys.glob(arms[[a]])
  acc <- rbindlist(lapply(fs, function(f) tryCatch({
    d <- fread(f, select = c("STAND_CN","PROJ_YEAR","BA_FT2AC"))
    d <- d[STAND_CN %in% keep_cn & is.finite(BA_FT2AC)]
    d[, .(sumBA = sum(BA_FT2AC), n = .N), by = PROJ_YEAR]
  }, error = function(e) { logerr(paste("pass2", f), e); NULL })), fill = TRUE)
  acc <- acc[, .(meanBA = sum(sumBA) / sum(n)), by = PROJ_YEAR][, arm := a]
  gc(); acc[]
}))
setorder(traj, arm, PROJ_YEAR)
fwrite(traj, file.path(OUT, "matrix_paired_trajectory.csv"))

## ---- Output C: headless figures ----
cols <- c(b1_native = "#d62728", b2_native = "#e6ab02",
          b1_gompit = "#7b2d8e", b2_gompit = "#1f77b4")
tryCatch({
  png(file.path(OUT, "fig_matrix_paired.png"), width = 1300, height = 720, res = 300)
  par(mar = c(4, 4, 3, 1))
  plot(NA, xlim = range(traj$PROJ_YEAR), ylim = c(0, max(traj$meanBA) * 1.05),
       xlab = "projection year", ylab = "mean BA (ft2/ac)",
       main = sprintf("Stand-paired arms, common stands n = %d", n_common), cex.main = 0.9)
  for (a in names(cols)) { s <- traj[arm == a]
    lines(s$PROJ_YEAR, s$meanBA, col = cols[a], lwd = 2,
          lty = if (grepl("gompit", a)) 1 else 2) }
  legend("topleft", c("b1 native","b2 native","b1 + gompit","b2 + gompit"),
         col = cols[c("b1_native","b2_native","b1_gompit","b2_gompit")],
         lty = c(2, 2, 1, 1), lwd = 2, bty = "n", cex = 0.8)
  dev.off()
}, error = function(e) logerr("fig_traj", e))

tryCatch({
  dB1 <- common$BA_b1_native - common$BA_b1_gompit
  dB2 <- common$BA_b2_native - common$BA_b2_gompit
  png(file.path(OUT, "fig_paired_diff.png"), width = 1300, height = 620, res = 300)
  par(mar = c(4, 4, 2, 1))
  xr <- range(c(dB1, dB2), na.rm = TRUE)
  h1 <- hist(dB1, breaks = 60, plot = FALSE); h2 <- hist(dB2, breaks = 60, plot = FALSE)
  plot(h1, col = "#7b2d8e80", border = NA, xlim = xr,
       ylim = c(0, max(h1$counts, h2$counts)),
       xlab = "within-stand BA reduction, native - gompit (ft2/ac)",
       main = "Per-stand mortality effect", cex.main = 0.9)
  plot(h2, col = "#1f77b480", border = NA, add = TRUE)
  abline(v = 0, lty = 3); abline(v = mean(dB1, na.rm = TRUE), col = "#7b2d8e", lwd = 2)
  abline(v = mean(dB2, na.rm = TRUE), col = "#1f77b4", lwd = 2)
  legend("topright", c("b1 (species-free)","b2 (species-dep)"),
         fill = c("#7b2d8e80","#1f77b480"), border = NA, bty = "n", cex = 0.8)
  dev.off()
}, error = function(e) logerr("fig_diff", e))

cat("matrix_paired done\n")
