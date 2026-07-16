#!/usr/bin/env Rscript
# CONUS-wide assessment: compare the three equation methods (Greg, species-free b1,
# species-dependent b2) on AGB, BA, QMD, TPH, CCH across all variants and the
# projection horizon. Memory-careful: read one metrics file at a time, aggregate
# to (METHOD, VARIANT, PROJ_YEAR) medians, drop, gc. Headless figures at 300 dpi.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(2026)
elog <- function(e) cat(conditionMessage(e), "\n", file = "error_log.txt", append = TRUE)
P   <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
OUT <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/method_assessment"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

methods <- list(
  greg = list(dir = file.path(P, "out_conus_eq_greg"),        pat = "_greg_metrics.csv$"),
  b1   = list(dir = file.path(P, "out_conus_eq_b1_gompit"),   pat = "_conus_b1_gompit_metrics.csv$"),
  b2   = list(dir = file.path(P, "out_conus_eq_b2_gompit"),   pat = "_conus_b2_gompit_metrics.csv$")
)
mvars <- c("AGB_TONS_AC","BA_FT2AC","QMD_IN","TPH","CCH_MEAN")

agg_one <- function(f, method) {
  tryCatch({
    d <- fread(f, showProgress = FALSE)
    if (!all(c("VARIANT","PROJ_YEAR") %in% names(d))) return(NULL)
    have <- intersect(mvars, names(d))
    a <- d[, lapply(.SD, function(x) as.numeric(median(x, na.rm = TRUE))),
           by = .(VARIANT, PROJ_YEAR), .SDcols = have]
    a[, METHOD := method]; rm(d); gc(verbose = FALSE); a
  }, error = function(e) { elog(e); NULL })
}

res <- list()
for (m in names(methods)) {
  dir <- methods[[m]]$dir
  if (!dir.exists(dir)) { cat("skip", m, "(no dir yet)\n"); next }
  fs <- list.files(dir, pattern = methods[[m]]$pat, full.names = TRUE)
  cat(m, ":", length(fs), "variant files\n")
  for (f in fs) res[[length(res) + 1]] <- agg_one(f, m)
}
S <- rbindlist(res, fill = TRUE)
if (!nrow(S)) { cat("no data\n"); quit(status = 0) }
fwrite(S, file.path(OUT, "method_comparison_by_variant_year.csv"))

## CONUS-pooled trajectory: median across variants per METHOD x PROJ_YEAR
pooled <- S[, lapply(.SD, median, na.rm = TRUE), by = .(METHOD, PROJ_YEAR),
            .SDcols = intersect(mvars, names(S))]
fwrite(pooled, file.path(OUT, "method_comparison_conus_pooled.csv"))
cat("=== CONUS-pooled medians by method (final projection year) ===\n")
print(pooled[PROJ_YEAR == max(PROJ_YEAR)])

## figures: pooled trajectories, one panel per metric
mlt <- melt(pooled, id.vars = c("METHOD","PROJ_YEAR"),
            measure.vars = intersect(mvars, names(pooled)), variable.name = "metric")
p <- ggplot(mlt, aes(PROJ_YEAR, value, color = METHOD)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "CONUS-wide method comparison (median across 19 variants)",
       x = "Projection year", y = NULL, color = "Method") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
tryCatch({
  png(file.path(OUT, "method_comparison_conus_pooled.png"), width = 2400, height = 1500, res = 300)
  print(p); dev.off()
}, error = function(e) elog(e))
cat("DONE: wrote summary CSVs + method_comparison_conus_pooled.png to", OUT, "\n")
cat("methods present:", paste(unique(S$METHOD), collapse = ", "), "\n")
