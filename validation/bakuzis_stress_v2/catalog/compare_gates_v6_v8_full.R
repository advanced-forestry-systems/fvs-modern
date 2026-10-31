## compare_gates_v6_v8_full.R (October 30, 2026): catalog v8 (D5 plot-size fix, PR 144) against catalog v6, every
## verdict field under the primary and the alternative rule, the stress labels, and the gate metrics.
## Supersedes catalog_v6_v8_verdict_changes.csv from compare_gates_v6_v8.R, which compared admissible only.
## Run from fix_b_20261030/.
suppressMessages(library(data.table))
ne <- function(a, b) { a <- as.character(a); b <- as.character(b); (is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b) }
a6 <- fread("../catalog/gates_v6/catalog_gates_by_arm.csv"); a8 <- fread("gates_v8/catalog_gates_by_arm.csv")
r6 <- fread("../catalog/gates_v6/catalog_recommendation_by_variant.csv"); r8 <- fread("gates_v8/catalog_recommendation_by_variant.csv")
s6 <- fread("../catalog/stress_v6/stress_by_arm.csv"); s8 <- fread("stress_v8/stress_by_arm.csv")
v6 <- fread("../catalog/stress_v6/stress_by_variant.csv"); v8 <- fread("stress_v8/stress_by_variant.csv")
keep <- c("arm", "group", "admissible", "admissible_alt", "calibrated_distinct", "share_ba999_50", "share_sdi_50", "reineke_flag_rate", "reineke_flag_wil", "reineke_eval", "med_ba50", "why_not")
m <- merge(a6[, ..keep], a8[, ..keep], by = c("arm", "group"), suffixes = c("_v6", "_v8"), all = TRUE)
m <- merge(m, merge(s6[, .(arm, group, stress_v6 = verdict_stability, boot_v6 = boot_admissible_prob, grid_v6 = grid_admissible_share)],
                    s8[, .(arm, group, stress_v8 = verdict_stability, boot_v8 = boot_admissible_prob, grid_v8 = grid_admissible_share)], by = c("arm", "group"), all = TRUE),
           by = c("arm", "group"), all = TRUE)
for (k in c("share_ba999_50", "share_sdi_50", "reineke_flag_rate", "med_ba50")) m[, paste0("d_", k) := get(paste0(k, "_v8")) - get(paste0(k, "_v6"))]
setorder(m, group, arm); fwrite(m, "catalog_v6_v8_comparison_by_arm_full.csv")
rr <- merge(r6[, .(group, recommended_v6 = recommended, status_v6 = status, admissible_arms_v6 = admissible_arms, admissible_arms_alt_v6 = admissible_arms_alt, rule_sensitive_v6 = rule_sensitive)],
            r8[, .(group, recommended_v8 = recommended, status_v8 = status, admissible_arms_v8 = admissible_arms, admissible_arms_alt_v8 = admissible_arms_alt, rule_sensitive_v8 = rule_sensitive)], by = "group", all = TRUE)
rr <- merge(rr, merge(v6[, .(group, boot_any_v6 = boot_any_admissible_prob, grid_any_v6 = grid_any_share)], v8[, .(group, boot_any_v8 = boot_any_admissible_prob, grid_any_v8 = grid_any_share)], by = "group", all = TRUE), by = "group", all = TRUE)
fwrite(rr, "catalog_v6_v8_recommendation_full.csv")
e <- m[!grepl("^conus_", arm)]
ch <- rbind(
  e[ne(admissible_v6, admissible_v8), .(group, arm, field = "admissible", v6 = as.character(admissible_v6), v8 = as.character(admissible_v8))],
  m[ne(admissible_alt_v6, admissible_alt_v8), .(group, arm, field = "admissible_alt", v6 = as.character(admissible_alt_v6), v8 = as.character(admissible_alt_v8))],
  e[ne(calibrated_distinct_v6, calibrated_distinct_v8), .(group, arm, field = "calibrated_distinct", v6 = as.character(calibrated_distinct_v6), v8 = as.character(calibrated_distinct_v8))],
  e[ne(stress_v6, stress_v8), .(group, arm, field = "stress_verdict", v6 = stress_v6, v8 = stress_v8)],
  rr[ne(recommended_v6, recommended_v8), .(group, arm = "(variant)", field = "recommended", v6 = recommended_v6, v8 = recommended_v8)],
  rr[ne(status_v6, status_v8), .(group, arm = "(variant)", field = "status", v6 = status_v6, v8 = status_v8)],
  rr[ne(admissible_arms_v6, admissible_arms_v8), .(group, arm = "(variant)", field = "admissible_arms", v6 = admissible_arms_v6, v8 = admissible_arms_v8)],
  rr[ne(admissible_arms_alt_v6, admissible_arms_alt_v8), .(group, arm = "(variant)", field = "admissible_arms_alt", v6 = admissible_arms_alt_v6, v8 = admissible_arms_alt_v8)],
  rr[ne(rule_sensitive_v6, rule_sensitive_v8), .(group, arm = "(variant)", field = "rule_sensitive", v6 = as.character(rule_sensitive_v6), v8 = as.character(rule_sensitive_v8))])
setorder(ch, group, arm, field); fwrite(ch, "catalog_v6_v8_verdict_changes_full.csv")
options(width = 250); cat("engine arms compared", nrow(e), "| field changes", nrow(ch), "\n"); print(ch)
adm <- e[admissible_v6 == TRUE & admissible_v8 == TRUE]
cat("largest |d med_ba50| among admissible arms:", round(max(abs(adm$d_med_ba50)), 3), adm[which.max(abs(d_med_ba50)), paste(group, arm)], "\n")
cat("largest |d med_ba50| any engine arm:", round(max(abs(e$d_med_ba50), na.rm = TRUE), 3), e[which.max(abs(d_med_ba50)), paste(group, arm)], "\n")
cat("stability counts v8:", paste(names(table(e$stress_v8)), table(e$stress_v8), collapse = "; "), "\n")
