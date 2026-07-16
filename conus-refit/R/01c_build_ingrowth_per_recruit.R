##=============================================================================
## 01c_build_ingrowth_per_recruit.R
##
## Builds the foundational per-recruit dataset required for a four-stage
## ingrowth model:
##   1. P(occurrence)   — any recruits per plot
##   2. N | occurrence  — count given recruits exist
##   3. DBH distribution — size of each recruit
##   4. Species composition — SPCD of each recruit
##
## Output is one row per ingrowth event (tree) with:
##   PLT_CN, SPCD, DBH_cm, TPA, state, COMPONENT_label
##
## Source: TREE_GRM_COMPONENT (DIA_END is the recruit's first measured
## diameter; SUBPTYP_BEGIN/MIDPT/END flags identify subplot type) joined
## to TREE.csv (CN -> TRE_CN) to attach SPCD.
##
## Aggregations are NOT done here. The downstream model code builds:
##   - plot-level (PLT_CN, n_recruits, has_recruit) for stages 1+2
##   - recruit-level (DBH_cm) for stage 3
##   - (PLT_CN, SPCD) counts for stage 4
##=============================================================================

suppressPackageStartupMessages({ library(data.table) })

FIA_DIR  <- "/users/PUOM0008/crsfaaron/fia_data"
ALT_DIR  <- "/users/PUOM0008/crsfaaron/fvs-conus/data/raw_fia"
OUT_FILE <- "/users/PUOM0008/crsfaaron/fvs-conus/data/ingrowth_per_recruit.rds"

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

grm_files <- c(
  list.files(FIA_DIR, pattern = "TREE_GRM_COMPONENT\\.csv$", full.names = TRUE),
  list.files(ALT_DIR, pattern = "TREE_GRM_COMPONENT\\.csv$", full.names = TRUE)
)
grm_files <- grm_files[!grepl("MIDPT", grm_files)]
get_state <- function(f) sub("_TREE_GRM_COMPONENT\\.csv$", "",
                              sub("^TREE_GRM_COMPONENT_", "", basename(f)))
grm_dt <- data.table(file = grm_files, state = vapply(grm_files, get_state, character(1)))
grm_dt <- grm_dt[!duplicated(state)]
cat(sprintf("Discovered %d state GRM files\n", nrow(grm_dt)))

INGROWTH_LABELS <- c("INGROWTH", "REVERSION1", "REVERSION2")
GRM_COLS <- c("TRE_CN", "PLT_CN", "STATECD", "DIA_END",
              "SUBP_COMPONENT_AL_FOREST", "SUBP_TPAGROW_UNADJ_AL_FOREST")

per_state <- list()
for (i in seq_len(nrow(grm_dt))) {
  st <- grm_dt$state[i]
  grm_f <- grm_dt$file[i]
  tree_f <- c(file.path(FIA_DIR, paste0(st, "_TREE.csv")),
              file.path(ALT_DIR, paste0(st, "_TREE.csv")))
  tree_f <- tree_f[file.exists(tree_f)][1]
  if (is.na(tree_f)) { message("  skip ", st, ": no TREE.csv"); next }

  grm <- tryCatch(fread(grm_f, select = GRM_COLS, showProgress = FALSE),
                  error = function(e) { message("  GRM fail ", st); NULL })
  if (is.null(grm) || nrow(grm) == 0) next
  grm <- grm[SUBP_COMPONENT_AL_FOREST %in% INGROWTH_LABELS]
  if (nrow(grm) == 0) next

  tree <- tryCatch(fread(tree_f, select = c("CN", "SPCD"), showProgress = FALSE),
                   error = function(e) { message("  TREE fail ", st); NULL })
  if (is.null(tree) || nrow(tree) == 0) { message("  skip ", st, ": no rows"); next }
  if (!all(c("CN", "SPCD") %in% names(tree))) {
    message("  skip ", st, ": TREE.csv missing CN or SPCD column (partial download)")
    next
  }

  setkey(tree, CN)
  rec <- tree[grm, on = c(CN = "TRE_CN"), nomatch = 0]
  # Convert DIA_END from inches to cm to match model conventions
  rec[, DBH_cm := DIA_END * 2.54]
  rec[, state := st]
  rec[, comp_label := SUBP_COMPONENT_AL_FOREST]
  cat(sprintf("  %s: %s recruits, %s plots, %s species, DBH_cm median=%.1f\n",
              st, format(nrow(rec), big.mark = ","),
              format(uniqueN(rec$PLT_CN), big.mark = ","),
              format(uniqueN(rec$SPCD), big.mark = ","),
              median(rec$DBH_cm, na.rm = TRUE)))
  per_state[[st]] <- rec[, .(PLT_CN, SPCD, DBH_cm,
                              tpa = SUBP_TPAGROW_UNADJ_AL_FOREST,
                              state, comp_label)]
}

all_rec <- rbindlist(per_state, use.names = TRUE, fill = TRUE)
cat(sprintf("\n=== CONUS-wide totals ===\n"))
cat(sprintf("  Per-recruit rows:   %s\n",
            format(nrow(all_rec), big.mark = ",")))
cat(sprintf("  Unique plots:       %s\n",
            format(uniqueN(all_rec$PLT_CN), big.mark = ",")))
cat(sprintf("  Unique species:     %s\n",
            format(uniqueN(all_rec$SPCD), big.mark = ",")))
cat(sprintf("  DBH_cm summary:\n"))
print(summary(all_rec$DBH_cm))

# Recruit DBH distribution overview by species (top 20)
cat("\nTop 20 species by recruit count + DBH stats:\n")
print(all_rec[!is.na(DBH_cm), .(
    n_recruits = .N,
    dbh_p5  = quantile(DBH_cm, 0.05),
    dbh_p50 = quantile(DBH_cm, 0.50),
    dbh_p95 = quantile(DBH_cm, 0.95)
  ), by = SPCD][order(-n_recruits)][1:20])

saveRDS(all_rec, OUT_FILE)
cat(sprintf("\nWrote: %s (%s rows)\n", OUT_FILE,
            format(nrow(all_rec), big.mark = ",")))
