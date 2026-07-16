# =============================================================================
# 01b_build_ingrowth_counts.R
#
# Builds plot-level counts of FIA ingrowth events from TREE_GRM_COMPONENT.
# Replaces the TREESTATUS-based approximation in 35_fit_ingrowth_negbinom.R
# which yields zero events because the remeasurement_pairs source only
# contains t1-live trees by construction.
#
# Strategy:
#   1. Auto-discover all *_TREE_GRM_COMPONENT.csv files in ~/fia_data/
#   2. Read SUBP_COMPONENT_AL_FOREST column (canonical FIA ingrowth flag for
#      subplot, all-live, on forest land).
#   3. Group by PLT_CN, count rows where flag is INGROWTH (and optionally
#      REVERSION1/REVERSION2 if reversions should also be counted).
#   4. Save as data/ingrowth_counts_by_plot.rds.
#
# The downstream ingrowth model script (35_fit_ingrowth_negbinom.R) then
# left-joins n_recruits from this file instead of computing from
# TREESTATUS1/TREESTATUS2.
#
# Usage:
#   Rscript --vanilla R/01b_build_ingrowth_counts.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

FIA_DIR  <- '/users/PUOM0008/crsfaaron/fia_data'
OUT_FILE <- '/users/PUOM0008/crsfaaron/fvs-conus/data/ingrowth_counts_by_plot.rds'

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

files <- list.files(FIA_DIR,
                    pattern = 'TREE_GRM_COMPONENT\\.csv$',
                    full.names = TRUE)
files <- files[!grepl('MIDPT', files)]
cat(sprintf('Discovered %d TREE_GRM_COMPONENT files\n', length(files)))
print(basename(files))

INGROWTH_LABELS <- c('INGROWTH', 'REVERSION1', 'REVERSION2')

read_one <- function(f) {
  st <- sub('_TREE_GRM_COMPONENT\\.csv$', '',
            sub('^TREE_GRM_COMPONENT_', '', basename(f)))
  cols <- c('PLT_CN', 'TRE_CN', 'STATECD',
            'SUBP_COMPONENT_AL_FOREST', 'SUBP_TPAGROW_UNADJ_AL_FOREST')
  d <- tryCatch(
    fread(f, select = cols, showProgress = FALSE),
    error = function(e) {
      message('  skip ', basename(f), ': ', conditionMessage(e))
      NULL
    }
  )
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d[, state := st]
  d
}

cat('Loading ...\n')
all_grm <- rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)
cat(sprintf('Total rows loaded: %s\n', format(nrow(all_grm), big.mark = ',')))

cat('Flag counts:\n')
print(all_grm[, .N, by = SUBP_COMPONENT_AL_FOREST][order(-N)])

ingrowth <- all_grm[SUBP_COMPONENT_AL_FOREST %in% INGROWTH_LABELS]
cat(sprintf('Ingrowth rows: %s across %s plots\n',
            format(nrow(ingrowth), big.mark = ','),
            format(length(unique(ingrowth$PLT_CN)), big.mark = ',')))

counts <- ingrowth[, .(n_recruits = .N,
                       n_recruits_tpa = sum(SUBP_TPAGROW_UNADJ_AL_FOREST,
                                            na.rm = TRUE)),
                   by = PLT_CN]
cat(sprintf('Plots with at least 1 recruit: %s\n',
            format(nrow(counts), big.mark = ',')))
cat('Per-plot summary:\n')
print(summary(counts$n_recruits))

saveRDS(counts, OUT_FILE)
cat(sprintf('Wrote: %s\n', OUT_FILE))
