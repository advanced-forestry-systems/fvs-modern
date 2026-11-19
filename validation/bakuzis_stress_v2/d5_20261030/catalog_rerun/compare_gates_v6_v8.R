## compare_gates_v6_v8.R (October 30, 2026): catalog v8 (D5 plot-size expansion, all arms rerun) against catalog v6.
suppressMessages(library(data.table)); C <- "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog"
a6 <- fread(file.path(C, "gates_v6/catalog_gates_by_arm.csv")); a8 <- fread("gates_v8/catalog_gates_by_arm.csv")
r6 <- fread(file.path(C, "gates_v6/catalog_recommendation_by_variant.csv")); r8 <- fread("gates_v8/catalog_recommendation_by_variant.csv")
keep <- c("arm", "group", "admissible", "share_ba999_50", "share_sdi_50", "reineke_flag_rate", "med_ba50", "why_not")
m <- merge(a6[, ..keep], a8[, ..keep], by = c("arm", "group"), suffixes = c("_v6", "_v8"), all = TRUE)[!grepl("^conus_", arm)]
for (k in c("share_ba999_50", "share_sdi_50", "reineke_flag_rate", "med_ba50")) m[, paste0("d_", k) := get(paste0(k, "_v8")) - get(paste0(k, "_v6"))]
setorder(m, group, arm); fwrite(m, "catalog_v6_v8_comparison_by_arm.csv")
rr <- merge(r6[, .(group, recommended_v6 = recommended, status_v6 = status, admissible_arms_v6 = admissible_arms)],
            r8[, .(group, recommended_v8 = recommended, status_v8 = status, admissible_arms_v8 = admissible_arms)], by = "group", all = TRUE)
fwrite(rr, "catalog_v6_v8_recommendation.csv")
ch <- rbind(m[admissible_v6 != admissible_v8, .(group, arm, field = "admissible", v6 = as.character(admissible_v6), v8 = as.character(admissible_v8))],
            rr[recommended_v6 != recommended_v8, .(group, arm = "(variant)", field = "recommended", v6 = recommended_v6, v8 = recommended_v8)],
            rr[admissible_arms_v6 != admissible_arms_v8, .(group, arm = "(variant)", field = "admissible_arms", v6 = admissible_arms_v6, v8 = admissible_arms_v8)])
fwrite(ch, "catalog_v6_v8_verdict_changes.csv"); options(width = 250)
cat("verdict changes:", nrow(ch), "\n"); print(ch)
print(m[, .(max_abs_d_ba999_50 = max(abs(d_share_ba999_50), na.rm = TRUE), max_abs_d_reineke = max(abs(d_reineke_flag_rate), na.rm = TRUE), max_abs_d_med_ba50 = max(abs(d_med_ba50), na.rm = TRUE))])
