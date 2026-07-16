#!/usr/bin/env Rscript
# 51e_compare_gompit_modifiers.R
# 3-way PSIS-LOO comparison of the gompit disturbance modifier legs fitted on the
# SAME species-dependent survival base (eta_base with z_sp): common (global alpha
# only) vs traitmed (alpha + W*gamma) vs speciesdep (alpha + W*gamma + z_sp).
# Reads modifier_loo.rds from each leg's output dir and runs loo::loo_compare.
suppressPackageStartupMessages({ library(loo); library(data.table) })
base <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/mort"
legs <- c(common   = file.path(base, "gompit_modifier_common",   "modifier_loo.rds"),
          traitmed = file.path(base, "gompit_modifier_traitmed", "modifier_loo.rds"),
          speciesdep = file.path(base, "gompit_modifier_speciesdep", "modifier_loo.rds"))
present <- legs[file.exists(legs)]
if (length(present) < 2) { cat("Need >=2 legs; found:", names(present), "\n"); quit(status = 0) }
loos <- lapply(present, readRDS)
cmp  <- loo_compare(loos)
cat("=== loo_compare (top row = best ELPD) ===\n"); print(cmp)
out <- data.table(model = rownames(cmp), as.data.table(cmp))
fwrite(out, file.path(base, "gompit_modifier_3way_elpd.csv"))
# Verdict: a leg is justified over the simpler one only if elpd_diff exceeds ~2*se_diff.
best <- rownames(cmp)[1]
cat(sprintf("\nBEST: %s\n", best))
for (i in 2:nrow(cmp)) {
  cat(sprintf("  %s vs best: elpd_diff=%.1f  se_diff=%.1f  %s\n",
              rownames(cmp)[i], cmp[i, "elpd_diff"], cmp[i, "se_diff"],
              ifelse(abs(cmp[i,"elpd_diff"]) > 2*cmp[i,"se_diff"], "DISTINGUISHED", "within noise")))
}
cat("\nNote: per-disturbance leg choice also follows sigma_sp / gamma_alpha CIs from the\n",
    "speciesdep summary (insect carries the species-dependent signal; others ~ common).\n")
