#!/usr/bin/env Rscript
## catalog_gates.R  (2026-10-28)
## Configuration catalog gates for FVS variants x survival arm x parameter set.
##
## Inputs (repeatable, as --long=<file> and --traj=<file>):
##   --long  long-format trajectories: arm, group, traj, age, TPH, BAPH, QMD (the v2 scorer's
##           *_long_metric.csv, or conus_catalog_runs.py output). Units: ha, m2 ha-1, cm, years.
##   --traj  v2 scorer *_scorecard_by_trajectory.csv (for the RD-gated Reineke verdicts, T3)
##   --thresholds  csv: variant, N, ba_p50, ba_p90, ba_p99, ba_p999, sdi_p999_ha (FIA plot level)
##   --seedgate    optional csv rows with variant, arm (survival), config, status (OK / SEED_FAIL ...)
##   --out   output directory
## Pre-registered rule (2026-10-28, after the plan red team):
##   admissible  <-  share(max BA over 0 to 50 yr > FIA BA p99.9) < 5% AND its Wilson 95% upper bound < 10%
##               AND share(peak SDI over 0 to 50 yr > FIA SDI p99.9) < 5% AND Wilson upper < 10%
##               AND Reineke flag rate among RD-gated evaluated trajectories < 10% (NA if fewer than 10 evaluated)
##               AND median peak BA (0 to 100 yr) >= FIA BA p50 (the over-kill lower bound)
##   Among admissible arms the simplest wins (default < calibrated < conus_b1 < conus_b2; native < gompit
##   only when the near-term stocking evidence does not separate them). The near-term tie-break is applied
##   only where per-variant near-term evidence exists (Northeast cohort, Paper 2 Table 1, basal area bias:
##   native +10.78%, GOMPIT +4.62%), which prefers GOMPIT in NE. Elsewhere the catalog records the
##   admissible set and the simplest admissible arm, and flags the missing near-term evidence.
##   A variant with no admissible arm is NOT RECOMMENDED for long-horizon use; the least-bad arm by
##   50 yr BA p99.9 exceedance is named for information only.
suppressMessages({library(data.table)})
args <- commandArgs(TRUE)
ga <- function(k) sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE))
LONG <- ga("long"); TRAJ <- ga("traj"); TH <- ga("thresholds")[1]; SG <- ga("seedgate"); OUT <- ga("out")[1]
if (is.na(OUT)) OUT <- "catalog_out"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(length(LONG) > 0, !is.na(TH))

rd <- function(f) { d <- fread(f); d[, src := basename(f)]; d }
L <- rbindlist(lapply(LONG, rd), fill = TRUE)
need <- c("arm", "group", "traj", "age", "TPH", "BAPH", "QMD"); miss <- setdiff(need, names(L))
if (length(miss)) stop("long input missing: ", paste(miss, collapse = ","))
L[, `:=`(age = as.numeric(age), TPH = as.numeric(TPH), BAPH = as.numeric(BAPH), QMD = as.numeric(QMD), group = toupper(group))]
L <- L[is.finite(BAPH) & is.finite(TPH) & TPH > 0]
L[, SDI := TPH * (QMD / 25.4)^1.605]     # Reineke SDI, metric (trees ha-1 at QMD 25.4 cm)
## Output red team, October 29 (defect 4): stands that failed the year-0 seed gate are excluded from the gates.
## Engine rows carry seed_ok; projector rows do not and are kept.
if ("seed_ok" %in% names(L)) {
  n_before <- uniqueN(L[, .(arm, group, traj)]); L <- L[is.na(seed_ok) | seed_ok %in% c(TRUE, "TRUE", "True")]
  cat("seed gate filter: trajectories", n_before, "->", uniqueN(L[, .(arm, group, traj)]), "\n")
}
## Output red team (defect 2): a calibrated arm whose trajectories equal the default arm's to the last digit was never
## distinct (the calibrated config did not load or carried nothing for the variant) and must not be counted as evidence.
dist <- L[grepl("^(default|calibrated)_", arm), .(arm, group, traj, age, BAPH, TPH)]
dist[, `:=`(pset = sub("_.*", "", arm), surv = sub(".*_", "", arm))]
dw <- dcast(dist, group + surv + traj + age ~ pset, value.var = c("BAPH", "TPH"))
nd <- dw[, .(max_abs_diff = max(abs(BAPH_calibrated - BAPH_default), abs(TPH_calibrated - TPH_default), na.rm = TRUE)), by = .(group, surv)]
nd[, `:=`(arm = paste0("calibrated_", surv), calibrated_distinct = max_abs_diff > 1e-9)]
cat("calibrated arms not distinct from default:", paste(nd[calibrated_distinct == FALSE, paste(group, surv)], collapse = ", "), "\n")
th <- fread(TH); setnames(th, "variant", "group"); th[, group := toupper(group)]

wilson_hi <- function(x, n, z = 1.96) { p <- x / n; (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n) }

## per trajectory peaks over the two horizons
pk <- L[, .(ba_max50 = max(BAPH[age <= 50]), ba_max100 = max(BAPH), sdi_max50 = max(SDI[age <= 50]), sdi_max100 = max(SDI),
            tph0 = TPH[which.min(age)], tph50 = if (any(age >= 50)) TPH[which.min(abs(age - 50))] else NA_real_,
            ba0 = BAPH[which.min(age)], ba50 = if (any(age >= 50)) BAPH[which.min(abs(age - 50))] else NA_real_,
            n_age = .N, age_max = max(age)), by = .(arm, group, traj)]
pk <- merge(pk, th, by = "group", all.x = TRUE)
g <- pk[, .(n_traj = .N, fia_n_plots = N[1],
            ba_p999 = ba_p999[1], ba_p99 = ba_p99[1], ba_p50 = ba_p50[1], sdi_p999 = sdi_p999_ha[1],
            x_ba999_50 = sum(ba_max50 > ba_p999), x_ba999_100 = sum(ba_max100 > ba_p999),
            x_ba99_50 = sum(ba_max50 > ba_p99), x_ba99_100 = sum(ba_max100 > ba_p99),
            x_sdi_50 = sum(sdi_max50 > sdi_p999_ha), x_sdi_100 = sum(sdi_max100 > sdi_p999_ha),
            med_peak_ba100 = median(ba_max100), med_peak_sdi100 = median(sdi_max100),
            med_ba0 = median(ba0), med_ba50 = median(ba50, na.rm = TRUE), med_stem_ret50 = median(tph50 / tph0, na.rm = TRUE),
            share_peak_below_p50 = mean(ba_max100 < ba_p50)), by = .(arm, group)]
g[, `:=`(share_ba999_50 = x_ba999_50 / n_traj, share_ba999_100 = x_ba999_100 / n_traj,
         share_ba99_50 = x_ba99_50 / n_traj, share_ba99_100 = x_ba99_100 / n_traj,
         share_sdi_50 = x_sdi_50 / n_traj, share_sdi_100 = x_sdi_100 / n_traj)]
g[, `:=`(wil_ba999_50 = wilson_hi(x_ba999_50, n_traj), wil_sdi_50 = wilson_hi(x_sdi_50, n_traj))]

## Reineke from the scorer's by-trajectory file (T3, RD gated)
if (length(TRAJ)) {
  T <- rbindlist(lapply(TRAJ, fread), fill = TRUE)
  T <- T[grepl("^T3 Reineke", relation)]; T[, group := toupper(group)]
  keep <- unique(L[, .(arm, group, traj)]); T <- merge(T, keep, by = c("arm", "group", "traj"))   # same seed-OK trajectory set as the magnitude gates
  r <- T[, .(reineke_eval = sum(status %in% c("PASS", "FLAG")), reineke_flag = sum(status == "FLAG"),
             reineke_na = sum(grepl("NOT APPLICABLE", status)),
             reineke_slope_med = median(as.numeric(value)[status %in% c("PASS", "FLAG")], na.rm = TRUE)), by = .(arm, group)]
  r[, reineke_flag_rate := ifelse(reineke_eval >= 10, reineke_flag / reineke_eval, NA_real_)]
  r[, reineke_eval_share := reineke_eval / (reineke_eval + reineke_na)]
  g <- merge(g, r, by = c("arm", "group"), all.x = TRUE)
} else g[, `:=`(reineke_eval = NA_integer_, reineke_flag = NA_integer_, reineke_flag_rate = NA_real_, reineke_eval_share = NA_real_, reineke_slope_med = NA_real_)]

## seed gate share (engine runs only)
if (length(SG)) {
  S <- rbindlist(lapply(SG, fread), fill = TRUE)
  if (all(c("variant", "arm", "config", "status") %in% names(S))) {
    S[, `:=`(group = toupper(variant), arm2 = paste(config, arm, sep = "_"))]
    s <- S[!is.na(config) & status %in% c("OK", "SEED_FAIL"), .(seed_n = .N, seed_ok_share = mean(status == "OK")), by = .(arm = arm2, group)]
    s2 <- S[!is.na(config), .(run_fail_share = mean(status == "RUN_FAIL")), by = .(arm = arm2, group)]
    s <- merge(s, s2, by = c("arm", "group"), all = TRUE)
    g <- merge(g, s, by = c("arm", "group"), all.x = TRUE)
  }
}

## admissibility
g <- merge(g, nd[, .(arm, group, calibrated_distinct)], by = c("arm", "group"), all.x = TRUE)
g[is.na(calibrated_distinct), calibrated_distinct := TRUE]
g[, reineke_unknown := is.na(reineke_flag_rate)]
## Repairs after the October 29 output red team (logged in implementation-notes.md): the under-stocking bound is the
## median basal area at year 50 (the 100 year peak included year 0 and could not fail); a Reineke gate that cannot
## be tested (fewer than 10 evaluated trajectories) is NOT TESTED and does not admit; a calibrated arm identical to the
## default is not evidence.
g[, admissible := share_ba999_50 < 0.05 & wil_ba999_50 < 0.10 & share_sdi_50 < 0.05 & wil_sdi_50 < 0.10 &
                   !reineke_unknown & reineke_flag_rate < 0.10 & med_ba50 >= ba_p50 & calibrated_distinct]
## Stress-test alternative (pre-registered as the protocol's alternative, not the primary rule): the Reineke gate on
## its Wilson 95 percent upper bound under 25 percent with at least 30 evaluated trajectories (fewer: not gating),
## because 10 to 12 evaluated trajectories make a 10 percent point gate flip on one trajectory.
g[, reineke_flag_wil := ifelse(!is.na(reineke_eval) & reineke_eval > 0, wilson_hi(reineke_flag, reineke_eval), NA_real_)]
g[, admissible_alt := share_ba999_50 < 0.05 & wil_ba999_50 < 0.10 & share_sdi_50 < 0.05 & wil_sdi_50 < 0.10 &
                       !is.na(reineke_eval) & reineke_eval >= 30 & reineke_flag_wil < 0.25 & med_ba50 >= ba_p50 & calibrated_distinct]
g[, why_not := paste0(ifelse(share_ba999_50 >= 0.05 | wil_ba999_50 >= 0.10, "BA>p99.9 by 50 yr; ", ""),
                      ifelse(share_sdi_50 >= 0.05 | wil_sdi_50 >= 0.10, "SDI>p99.9 by 50 yr; ", ""),
                      ifelse(!is.na(reineke_flag_rate) & reineke_flag_rate >= 0.10, "Reineke flags; ", ""),
                      ifelse(reineke_unknown, sprintf("Reineke not testable (%s evaluated); ", ifelse(is.na(reineke_eval), "0", reineke_eval)), ""),
                      ifelse(is.na(med_ba50) | med_ba50 < ba_p50, "under-stocks (median BA at 50 yr < FIA p50); ", ""),
                      ifelse(!calibrated_distinct, "calibrated config did not load (identical to default); ", ""))]
g[, threshold_n := fia_n_plots]

## simplicity order and recommendation
simp <- function(a) { s <- 0
  s <- s + ifelse(grepl("^default", a), 0, ifelse(grepl("^calibrated", a), 1, ifelse(grepl("b1", a), 2, ifelse(grepl("b2", a), 3, 4))))
  s * 10 + ifelse(grepl("gompit", a), 1, 0) }
g[, simplicity := simp(arm)]
nearterm <- c(NE = "gompit")   # Paper 2 Table 1 (Northeast cohort): GOMPIT BA bias +4.62% vs native +10.78%
rec <- g[!is.na(ba_p999), {
  ad <- .SD[admissible == TRUE]
  ad_alt <- .SD[admissible_alt == TRUE]
  sens <- !setequal(ad$arm, ad_alt$arm)
  if (nrow(ad)) {
    pref <- if (!is.na(nearterm[group[1]])) ad[grepl(nearterm[group[1]], arm)] else ad
    if (!nrow(pref)) pref <- ad
    best <- pref[order(simplicity)][1]
    list(recommended = best$arm, status = if (best$n_traj < 100) sprintf("RECOMMENDED (low sample: %d stands)", best$n_traj) else "RECOMMENDED", admissible_arms = paste(ad[order(simplicity), arm], collapse = "; "),
         nearterm_evidence = ifelse(is.na(nearterm[group[1]]), "none per variant", "Paper 2 NE cohort (stocking)"),
         rule_sensitive = sens, admissible_arms_alt = paste(ad_alt[order(simplicity), arm], collapse = "; "))
  } else {
    best <- .SD[order(share_ba999_50, share_sdi_50)][1]
    list(recommended = best$arm, status = if (max(n_traj) < 100) sprintf("INSUFFICIENT SAMPLE (%d stands); least-bad arm named", max(n_traj)) else "NOT RECOMMENDED for long horizons (least-bad arm named)", admissible_arms = "",
         nearterm_evidence = ifelse(is.na(nearterm[group[1]]), "none per variant", "Paper 2 NE cohort (stocking)"),
         rule_sensitive = sens, admissible_arms_alt = paste(ad_alt[order(simplicity), arm], collapse = "; "))
  }
}, by = group]
noth <- g[is.na(ba_p999), .(recommended = "", status = "NOT EVALUATED (no FIA thresholds for this variant)", admissible_arms = "",
                            nearterm_evidence = "none per variant", rule_sensitive = FALSE, admissible_arms_alt = ""), by = group]
rec <- rbind(rec, noth, fill = TRUE)
setorder(g, group, simplicity)
fwrite(g, file.path(OUT, "catalog_gates_by_arm.csv"))
fwrite(rec, file.path(OUT, "catalog_recommendation_by_variant.csv"))
cat("arms x variants:", nrow(g), " variants:", uniqueN(g$group), "\n")
print(rec)
print(g[, .(group, arm, n_traj, share_ba999_50 = round(share_ba999_50, 3), share_sdi_50 = round(share_sdi_50, 3),
            reineke_flag_rate = round(reineke_flag_rate, 2), med_peak_ba100 = round(med_peak_ba100, 1), ba_p50 = round(ba_p50, 1), admissible)], nrows = 200)
