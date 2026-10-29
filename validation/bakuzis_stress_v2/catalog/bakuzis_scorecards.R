#!/usr/bin/env Rscript
## bakuzis_scorecards.R (2026-10-29)
## One Bakuzis scorecard per FVS variant from the v2 scorer outputs on the catalog trajectories (four fvs-modern arms on
## the common seed-passing stand sample, with within-variant site and density terciles) plus the four CONUS projector arms.
## Inputs: --traj (scorecard_by_trajectory.csv), --gates (catalog_gates_by_arm.csv), --thresholds, --out (directory),
##         --step (projection step, years), --qmd_rate (cm yr-1 used by the scorer), --label
## Output: <out>/<VARIANT>_bakuzis_scorecard.md, <out>/scorecard_index.md, <out>/scorecard_cells.csv
suppressMessages(library(data.table))
args <- commandArgs(TRUE)
ga <- function(k, d = NA) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v else d }
T <- fread(ga("traj")); G <- fread(ga("gates")); TH <- fread(ga("thresholds")); OUT <- ga("out", "scorecards"); STEP <- ga("step", "5"); QR <- ga("qmd_rate", "0.10"); LAB <- ga("label", "catalog v3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
T[, group := toupper(group)]
T[, rel := sub("^(T[0-9]+b?) .*", "\\1", relation)]
T[, short := fifelse(status %in% c("PASS", "FLAG"), status, fifelse(grepl("NOT APPLICABLE", status), "NA", fifelse(grepl("NOT EVALUATED", status), "NE", status)))]
cells <- T[, .(evaluated = sum(short %in% c("PASS", "FLAG")), flagged = sum(short == "FLAG"), not_applicable = sum(short == "NA"), not_evaluated = sum(short == "NE")),
           by = .(group, arm, rel, relation)]
cells[, verdict := fifelse(evaluated == 0, "NOT EVALUATED", fifelse(flagged == 0, "PASS", fifelse(flagged == evaluated, "FLAG (all)", "FLAG")))]
fwrite(cells, file.path(OUT, "scorecard_cells.csv"))
arm_order <- c("default_native", "default_gompit", "calibrated_native", "calibrated_gompit", "conus_b1", "conus_b1_gompit", "conus_b2", "conus_b2_gompit")
arm_lab <- c(default_native = "default, native", default_gompit = "default, GOMPIT", calibrated_native = "calibrated, native", calibrated_gompit = "calibrated, GOMPIT",
             conus_b1 = "CONUS b1, native", conus_b1_gompit = "CONUS b1, GOMPIT", conus_b2 = "CONUS b2, native", conus_b2_gompit = "CONUS b2, GOMPIT")
rel_order <- c("T0", "T1", "T2", "T3", "T3b", "T4", "T5", "T6", "T7", "T8", "T9")
rel_name <- c(T0 = "Initial-step shock", T1 = "Monotone height and QMD", T2 = "Site ordering (height, volume)", T3 = "Reineke slope, within stand (RD gated)",
              T3b = "Reineke frontier (tau 0.95)", T4 = "Magnitude vs FIA envelope (BA p99, p99.9; SDI)", T5 = "Eichhorn (volume at matched height across sites)",
              T6 = "Sukachev (per-tree volume vs density)", T7 = "Yield-density (final BA vs initial density)", T8 = "Basal area plateau", T9 = "BA and volume inversions")
## Component attribution per relation, from references/bakuzis-method.md, written once and reused.
attrib <- c(
  T0 = "first-cycle mortality on dense stands: the mortality routine (native SDI-based or GOMPIT crown-closure) removes a large share of stems in the first step; check the density term and the cycle-1 calibration against the observed first remeasurement",
  T1 = "diameter or height increment turning negative: the increment equations or the summary accounting; a monotone failure is a defect, not a calibration question",
  T2 = "height submodel site term too weak or inverted (site curves cross), or volume ordering lost downstream of self-thinning; on real FIA stands the site classes are terciles of SICOND, so species and age mix are confounded with site and a flag here is a design signal before it is a model one",
  T3 = "mortality and self-thinning calibration: SDImax, the Reineke slope used by the density-dependent mortality, or GOMPIT's crown-closure response; a shallow slope means stems are not removed as stands densify",
  T3b = "the maximum-density line across the sample: SDImax level and slope; a pass here with T3 flags means the cap shapes the densest stands while typical stands do not self-thin",
  T4 = "the growth and mortality balance against the observed envelope: basal area or SDI beyond anything measured in the fitting data points at the mortality submodel (native) or the BA cap; GOMPIT typically restores it",
  T5 = "the volume, height and basal area chain: a site-dependent volume-height relation points at the form factor or the BA-height-site interaction; on real stands with tercile site classes, species mix is confounded with site",
  T6 = "mortality's site and density response: per-tree volume should fall with density at matched age; on real stands the density terciles differ in age and structure, so read this with T7",
  T7 = "self-thinning calibration: a lower final BA from a higher initial density means the density-dependent mortality over-corrects (or the stands were older and already declining); check SDImax and the mortality multiplier",
  T8 = "a pinned plateau (CV under 3 percent) means a hard cap (BAMAX or SDImax) is binding and the stand is held there rather than self-thinning along a trajectory",
  T9 = "single-step drops above 10 percent in BA or volume: mortality pulses (native SDI mortality firing in one cycle) or a cap engaging; GOMPIT arms carry more of these because mortality arrives earlier"
)
groups <- sort(unique(G[!is.na(ba_p999), group]))
index <- c(sprintf("# Bakuzis scorecards by FVS variant (%s)", LAB), "",
           sprintf("Projection step %s yr, Reineke screen qmd_rate %s cm yr-1, relative-density gate 0.55, thresholds from plot-level FIA fitting data by variant. One file per variant; `scorecard_cells.csv` has every cell.", STEP, QR), "",
           "| Variant | Stands (engine) | Recommended | Status | Relations flagged on the recommended arm |", "|---|---|---|---|---|")
rec <- fread(sub("catalog_gates_by_arm.csv", "catalog_recommendation_by_variant.csv", ga("gates")))
for (gv in groups) {
  g <- G[group == gv]; th <- TH[toupper(variant) == gv]
  cv <- cells[group == gv]
  arms <- intersect(arm_order, unique(c(g$arm, cv$arm)))
  r <- rec[group == gv]
  n_eng <- max(g[!grepl("^conus", arm), n_traj], 0)
  lines <- c(sprintf("# %s: Bakuzis scorecard (%s)", gv, LAB), "",
             sprintf("Engine arms: %d seed-passing FIA stands, 100 yr at a %s yr step, no regeneration, posterior draws omitted. Projector arms: 1,500 stands each from the October 27 rescore (no dominant height, so the height, site, Eichhorn and Sukachev relations are NOT EVALUATED there). Site and density classes for the engine arms are within-variant terciles of SICOND and BALIVE. Thresholds: BA p99 %.1f and p99.9 %.1f m2 ha-1, SDI p99.9 %.0f trees ha-1 from %d fitting plots%s. Reineke screen qmd_rate %s cm yr-1, relative-density gate 0.55.",
                     n_eng, STEP, th$ba_p99, th$ba_p999, th$sdi_p999_ha, th$N, ifelse(th$N < 1000, " (fewer than 1,000 plots, so the p99.9 is the sample maximum)", ""), QR), "",
             sprintf("Catalog verdict (primary rule): **%s**; recommended arm %s; admissible: %s%s.", r$status, ifelse(nchar(r$recommended), r$recommended, "none"),
                     ifelse(nchar(r$admissible_arms), r$admissible_arms, "none"), ifelse(isTRUE(r$rule_sensitive), "; rule-sensitive (alternative Reineke gate admits a different set)", "")), "",
             "## Relation table (evaluated / flagged; NE = not evaluated, NA = not applicable below the relative-density gate)", "")
  hdr <- paste0("| Relation | ", paste(arm_lab[arms], collapse = " | "), " |"); sep <- paste0("|---|", paste(rep("---", length(arms)), collapse = "|"), "|")
  lines <- c(lines, hdr, sep)
  for (rl in rel_order) {
    rr <- cv[rel == rl]; if (!nrow(rr)) next
    rels <- unique(rr$relation)
    for (rn in rels) {
      row <- sapply(arms, function(a) { x <- rr[arm == a & relation == rn]; if (!nrow(x)) "NE" else if (x$evaluated == 0) sprintf("NE (%d NA)", x$not_applicable) else sprintf("%d / **%d**", x$evaluated, x$flagged) })
      lines <- c(lines, paste0("| ", rn, " | ", paste(row, collapse = " | "), " |"))
    }
  }
  ## magnitude block from the gate table
  lines <- c(lines, "", "## Magnitude and gates (engine and projector arms)", "", "| Arm | n | BA > p99.9 by 50 yr (%) | Wilson upper (%) | SDI > p99.9 by 50 yr (%) | Reineke flag rate (evaluated) | median slope | median BA at 50 yr (m2 ha-1) vs p50 | stems retained at 50 yr | seed gate OK | admissible | admissible (alt) |", "|---|---|---|---|---|---|---|---|---|---|---|---|")
  for (a in arms) { x <- g[arm == a]; if (!nrow(x)) next
    lines <- c(lines, sprintf("| %s | %d | %.1f | %.1f | %.1f | %s | %s | %.1f vs %.1f | %.2f | %s | %s | %s |", arm_lab[a], x$n_traj, 100 * x$share_ba999_50, 100 * x$wil_ba999_50, 100 * x$share_sdi_50,
                              ifelse(is.na(x$reineke_flag_rate), sprintf("not testable (%s)", ifelse(is.na(x$reineke_eval), 0, x$reineke_eval)), sprintf("%.2f (%d)", x$reineke_flag_rate, x$reineke_eval)),
                              ifelse(is.na(x$reineke_slope_med), "", sprintf("%.2f", x$reineke_slope_med)), x$med_ba50, x$ba_p50, x$med_stem_ret50,
                              ifelse(is.na(x$seed_ok_share), "", sprintf("%.2f", x$seed_ok_share)), ifelse(isTRUE(x$admissible), "yes", "no"), ifelse(isTRUE(x$admissible_alt), "yes", "no"))) }
  ## attribution: relations flagged on the recommended (or least-bad engine) arm
  ra <- if (nchar(r$recommended) && !grepl("^conus", r$recommended)) r$recommended else g[!grepl("^conus", arm)][order(share_ba999_50, share_sdi_50)][1, arm]
  fl <- cv[arm == ra & verdict != "PASS" & verdict != "NOT EVALUATED"][order(match(rel, rel_order))]
  lines <- c(lines, "", sprintf("## Reading and component attribution (arm %s)", arm_lab[ra]), "")
  if (!nrow(fl)) lines <- c(lines, "No relation is flagged on this arm among those evaluated. That is a failure to falsify, not a confirmation of accuracy; the relations that could not be evaluated are listed above as NE.")
  for (i in seq_len(nrow(fl))) {
    x <- fl[i]
    lines <- c(lines, sprintf("- **%s** (%s): %d of %d evaluated flagged. %s.", rel_name[x$rel], x$relation, x$flagged, x$evaluated, attrib[x$rel]))
  }
  strong <- fl[rel %in% c("T3", "T5")]
  lines <- c(lines, "", if (nrow(strong)) sprintf("Strong relations (Eichhorn, Reineke) flagged: %s. These outrank the others; start there.", paste(unique(strong$rel), collapse = ", ")) else "Neither strong relation (Eichhorn, Reineke) is flagged on this arm among evaluated trajectories.",
             "", "Figures: the scorer's matrix panels for this variant are `*_<arm>_<VARIANT>_bakuzis_matrix.png` and `*_<arm>_selfthinning.png` beside this file.", "")
  writeLines(lines, file.path(OUT, sprintf("%s_bakuzis_scorecard.md", gv)))
  index <- c(index, sprintf("| %s | %d | %s | %s | %s |", gv, n_eng, ifelse(nchar(r$recommended), r$recommended, "none"), r$status, ifelse(nrow(fl), paste(fl$rel, collapse = ", "), "none")))
}
writeLines(index, file.path(OUT, "scorecard_index.md"))
cat("wrote", length(groups), "scorecards to", OUT, "\n")
