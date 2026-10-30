## compare_gates_v5_v6.R (October 30, 2026): catalog v6 (calibrated arms rerun with the PR 135 overlay fix)
## against catalog v5, per variant and arm, and every verdict change.
suppressMessages(library(data.table))
a5 <- fread("gates_v5/catalog_gates_by_arm.csv"); a6 <- fread("gates_v6/catalog_gates_by_arm.csv")
r5 <- fread("gates_v5/catalog_recommendation_by_variant.csv"); r6 <- fread("gates_v6/catalog_recommendation_by_variant.csv")
s5 <- fread("stress_v5/stress_by_arm.csv"); s6 <- fread("stress_v6/stress_by_arm.csv")
keep <- c("arm", "group", "admissible", "admissible_alt", "calibrated_distinct", "share_ba999_50", "share_sdi_50", "reineke_flag_rate", "med_ba50", "why_not")
m <- merge(a5[, ..keep], a6[, ..keep], by = c("arm", "group"), suffixes = c("_v5", "_v6"), all = TRUE)
st <- merge(s5[, .(arm, group, stress_v5 = verdict_stability, boot_v5 = boot_admissible_prob)],
            s6[, .(arm, group, stress_v6 = verdict_stability, boot_v6 = boot_admissible_prob)], by = c("arm", "group"), all = TRUE)
m <- merge(m, st, by = c("arm", "group"), all = TRUE)
m <- m[!grepl("^conus_", arm)]
num <- c("share_ba999_50", "share_sdi_50", "reineke_flag_rate", "med_ba50")
for (k in num) m[, paste0("d_", k) := get(paste0(k, "_v6")) - get(paste0(k, "_v5"))]
setorder(m, group, arm)
fwrite(m, "catalog_v5_v6_comparison_by_arm.csv")
rr <- merge(r5[, .(group, recommended_v5 = recommended, status_v5 = status, admissible_arms_v5 = admissible_arms, admissible_arms_alt_v5 = admissible_arms_alt)],
            r6[, .(group, recommended_v6 = recommended, status_v6 = status, admissible_arms_v6 = admissible_arms, admissible_arms_alt_v6 = admissible_arms_alt)], by = "group", all = TRUE)
fwrite(rr, "catalog_v5_v6_recommendation.csv")
ch <- rbind(
  m[admissible_v5 != admissible_v6, .(group, arm, field = "admissible", v5 = as.character(admissible_v5), v6 = as.character(admissible_v6))],
  m[admissible_alt_v5 != admissible_alt_v6, .(group, arm, field = "admissible_alt", v5 = as.character(admissible_alt_v5), v6 = as.character(admissible_alt_v6))],
  m[calibrated_distinct_v5 != calibrated_distinct_v6, .(group, arm, field = "calibrated_distinct", v5 = as.character(calibrated_distinct_v5), v6 = as.character(calibrated_distinct_v6))],
  m[stress_v5 != stress_v6, .(group, arm, field = "stress_verdict", v5 = stress_v5, v6 = stress_v6)],
  rr[recommended_v5 != recommended_v6, .(group, arm = "(variant)", field = "recommended", v5 = recommended_v5, v6 = recommended_v6)],
  rr[status_v5 != status_v6, .(group, arm = "(variant)", field = "status", v5 = status_v5, v6 = status_v6)],
  rr[admissible_arms_v5 != admissible_arms_v6, .(group, arm = "(variant)", field = "admissible_arms", v5 = admissible_arms_v5, v6 = admissible_arms_v6)])
setorder(ch, group, arm, field)
fwrite(ch, "catalog_v5_v6_verdict_changes.csv")
options(width = 250)
cat("arms compared", nrow(m), "| verdict changes", nrow(ch), "\n"); print(ch)
cat("\nlargest year-50 median BA changes (calibrated arms):\n")
print(m[grepl("^calibrated", arm)][order(-abs(d_med_ba50))][1:12, .(group, arm, med_ba50_v5, med_ba50_v6, d_med_ba50, share_ba999_50_v5, share_ba999_50_v6, share_sdi_50_v5, share_sdi_50_v6, reineke_flag_rate_v5, reineke_flag_rate_v6)], digits = 3)
cat("\ndefault arms identical v5 vs v6 (all numeric gate fields):", m[grepl("^default", arm), all(abs(d_med_ba50) < 1e-9 & abs(d_share_ba999_50) < 1e-9 & abs(d_share_sdi_50) < 1e-9, na.rm = TRUE)], "\n")
