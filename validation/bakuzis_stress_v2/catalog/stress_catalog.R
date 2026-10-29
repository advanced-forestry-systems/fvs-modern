#!/usr/bin/env Rscript
## stress_catalog.R (2026-10-29)
## Stress test of the shipped catalog (primary rule) against the analyst choices that could have gone another way.
## Two parts, both computed from the trajectory-level data without new engine runs:
##   A. Rule grid: horizon (50, 100 yr) x exceedance tolerance (3, 5, 10 percent) x Wilson upper bound (10, 15 percent)
##      x basal area percentile (p99.9, p99) x Reineke flag-rate tolerance (10, 15, 25 percent) x minimum evaluated
##      trajectories for the Reineke gate (10, 30). Reports, per arm and variant, the share of the 144 grid cells in which
##      the arm is admissible, and whether the shipped verdict sits in the majority.
##   B. Stand bootstrap: 500 resamples of stands with replacement per arm and variant, primary rule unchanged; reports
##      the probability that the arm is admissible and that the variant has at least one admissible engine arm.
## Inputs: --long (engine trajectories with seed_ok), --traj (scorer by-trajectory file, T3 statuses), --thresholds,
##         --gates (shipped gate table, for the shipped verdict), --out
suppressMessages(library(data.table))
args <- commandArgs(TRUE)
ga <- function(k, d = NA) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v else d }
L <- fread(ga("long")); T <- fread(ga("traj")); TH <- fread(ga("thresholds")); G <- fread(ga("gates")); OUT <- ga("out", "stress"); B <- as.integer(ga("boot", "500"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(20261029)
L <- L[seed_ok %in% c(TRUE, "TRUE", "True") & is.finite(BAPH) & is.finite(TPH) & TPH > 0]
L[, `:=`(group = toupper(group), age = as.numeric(age))]; L[, SDI := TPH * (QMD / 25.4)^1.605]
setnames(TH, "variant", "group"); TH[, group := toupper(group)]
pk <- L[, .(ba_max50 = max(BAPH[age <= 50]), ba_max100 = max(BAPH), sdi_max50 = max(SDI[age <= 50]), sdi_max100 = max(SDI),
            ba50 = if (any(age >= 50)) BAPH[which.min(abs(age - 50))] else NA_real_), by = .(arm, group, traj)]
T3 <- T[grepl("^T3 Reineke", relation), .(arm, group = toupper(group), traj, t3 = fifelse(status == "FLAG", 1L, fifelse(status == "PASS", 0L, NA_integer_)))]
pk <- merge(pk, T3, by = c("arm", "group", "traj"), all.x = TRUE)
pk <- merge(pk, TH, by = "group")
wil <- function(x, n, z = 1.96) { p <- x / n; (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n) }
adm <- function(d, horizon = 50, tol = 0.05, wmax = 0.10, pct = "p999", rk_tol = 0.10, rk_min = 10) {
  n <- nrow(d); if (n == 0) return(NA)
  bth <- if (pct == "p999") d$ba_p999[1] else d$ba_p99[1]
  bam <- if (horizon == 50) d$ba_max50 else d$ba_max100; sdm <- if (horizon == 50) d$sdi_max50 else d$sdi_max100
  xb <- sum(bam > bth); xs <- sum(sdm > d$sdi_p999_ha[1])
  ev <- sum(!is.na(d$t3)); fl <- sum(d$t3 == 1, na.rm = TRUE)
  rk_ok <- if (ev < rk_min) FALSE else (fl / ev) < rk_tol
  (xb / n < tol) & (wil(xb, n) < wmax) & (xs / n < tol) & (wil(xs, n) < wmax) & rk_ok & (median(d$ba50, na.rm = TRUE) >= d$ba_p50[1])
}
## A. rule grid
grid <- CJ(horizon = c(50, 100), tol = c(0.03, 0.05, 0.10), wmax = c(0.10, 0.15), pct = c("p999", "p99"), rk_tol = c(0.10, 0.15, 0.25), rk_min = c(10, 30))
cat("grid cells:", nrow(grid), "\n")
keys <- unique(pk[, .(arm, group)])
gridres <- rbindlist(lapply(seq_len(nrow(keys)), function(i) {
  d <- pk[arm == keys$arm[i] & group == keys$group[i]]
  r <- sapply(seq_len(nrow(grid)), function(j) adm(d, grid$horizon[j], grid$tol[j], grid$wmax[j], grid$pct[j], grid$rk_tol[j], grid$rk_min[j]))
  data.table(arm = keys$arm[i], group = keys$group[i], n = nrow(d), grid_admissible_share = mean(r, na.rm = TRUE),
             shipped = adm(d), admissible_h100 = adm(d, horizon = 100), admissible_p99 = adm(d, pct = "p99"), admissible_tol03 = adm(d, tol = 0.03),
             admissible_rk25 = adm(d, rk_tol = 0.25), admissible_rkmin30 = adm(d, rk_min = 30))
}))
## B. stand bootstrap under the shipped rule
bootres <- rbindlist(lapply(seq_len(nrow(keys)), function(i) {
  d <- pk[arm == keys$arm[i] & group == keys$group[i]]; n <- nrow(d)
  r <- replicate(B, adm(d[sample.int(n, n, replace = TRUE)]))
  data.table(arm = keys$arm[i], group = keys$group[i], boot_admissible_prob = mean(r, na.rm = TRUE))
}))
res <- merge(gridres, bootres, by = c("arm", "group"))
res <- merge(res, G[, .(arm, group, shipped_catalog = admissible)], by = c("arm", "group"), all.x = TRUE)
res[, agrees_with_catalog := shipped == shipped_catalog]
res[, verdict_stability := fifelse(shipped, fifelse(boot_admissible_prob >= 0.9 & grid_admissible_share >= 0.5, "STABLE admissible", "FRAGILE admissible"),
                                   fifelse(boot_admissible_prob <= 0.1 & grid_admissible_share <= 0.5, "STABLE not admissible", "FRAGILE not admissible"))]
setorder(res, group, arm)
fwrite(res, file.path(OUT, "stress_by_arm.csv"))
## variant level: probability of at least one admissible engine arm
vb <- rbindlist(lapply(unique(keys$group), function(gv) {
  arms <- keys[group == gv, arm]; dd <- lapply(arms, function(a) pk[arm == a & group == gv])
  n <- nrow(dd[[1]])
  r <- replicate(B, { idx <- sample.int(n, n, replace = TRUE); any(sapply(dd, function(d) adm(d[idx]))) })
  data.table(group = gv, n = n, shipped_any = any(sapply(dd, adm)), boot_any_admissible_prob = mean(r, na.rm = TRUE),
             grid_any_share = mean(sapply(seq_len(nrow(grid)), function(j) any(sapply(dd, function(d) adm(d, grid$horizon[j], grid$tol[j], grid$wmax[j], grid$pct[j], grid$rk_tol[j], grid$rk_min[j]))))))
}))
setorder(vb, group); fwrite(vb, file.path(OUT, "stress_by_variant.csv"))
cat("disagreements with shipped catalog:", sum(!res$agrees_with_catalog, na.rm = TRUE), "\n"); print(res[, .(group, arm, n, shipped, boot = round(boot_admissible_prob, 2), grid = round(grid_admissible_share, 2), verdict_stability)], nrows = 200); print(vb)
