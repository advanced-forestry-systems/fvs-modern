#!/usr/bin/env Rscript
# Scoring step for the GOMPIT mortality-cap Bakuzis validation.
# Parses the FVS "SUMMARY STATISTICS" table from each stand .out file,
# extracts basal area (BA, ft2/acre) at simulation ages 100/200/300 yr,
# converts to m2/ha, and writes a per-(config, stand, year) results CSV.
#
# Summary table is fixed-width; after the header the data rows lead with a
# 4-digit calendar YEAR then AGE, TREES, BA, ... (whitespace-separated,
# all numeric). Field 1 = YEAR, 2 = AGE, 3 = TREES, 4 = BA.

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "/fs/scratch/PUOM0008/crsfaaron/rescore_mortcap"

FT2AC_TO_M2HA <- 0.2295684        # 1 ft^2/acre = 0.2295684 m^2/ha
TARGET_AGES  <- c(100, 200, 300)  # horizons relative to simulation start
configs <- c("native", "allon_mortcap")

parse_out <- function(f) {
  ln <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  if (!length(ln)) return(NULL)
  i <- grep("SUMMARY STATISTICS", ln)
  if (!length(i)) return(NULL)
  sub <- ln[i[1]:length(ln)]
  # keep only data rows: leading 4-digit year then whitespace
  rows <- grep("^[0-9]{4}[[:space:]]", sub, value = TRUE)
  if (!length(rows)) return(NULL)
  recs <- lapply(rows, function(r) {
    fld <- strsplit(trimws(r), "[[:space:]]+")[[1]]
    if (length(fld) < 4) return(NULL)
    suppressWarnings(data.frame(
      year_cal = as.integer(fld[1]),
      age      = as.integer(fld[2]),
      trees    = as.numeric(fld[3]),
      ba_ft2ac = as.numeric(fld[4]),
      stringsAsFactors = FALSE
    ))
  })
  recs <- do.call(rbind, recs[!vapply(recs, is.null, logical(1))])
  recs
}

out <- list()
for (cfg in configs) {
  wdir <- file.path(BASE, "work", cfg)
  outs <- list.files(wdir, pattern = "\\.out$", full.names = TRUE)
  cat(sprintf("config=%s  .out files=%d\n", cfg, length(outs)))
  for (f in outs) {
    stand <- sub("\\.out$", "", basename(f))
    d <- parse_out(f)
    if (is.null(d)) { cat(sprintf("  no summary table: %s\n", stand)); next }
    d <- d[d$age %in% TARGET_AGES, , drop = FALSE]
    if (!nrow(d)) next
    out[[length(out) + 1]] <- data.frame(
      config   = cfg,
      stand    = stand,
      year     = d$age,               # horizon (yr since sim start)
      year_cal = d$year_cal,
      BA       = round(d$ba_ft2ac * FT2AC_TO_M2HA, 4),  # m2/ha
      BA_ft2ac = d$ba_ft2ac,
      trees_ac = d$trees,
      stringsAsFactors = FALSE
    )
  }
}

res <- if (length(out)) do.call(rbind, out) else
  data.frame(config = character(), stand = character(), year = integer(),
             year_cal = integer(), BA = numeric(), BA_ft2ac = numeric(),
             trees_ac = numeric())

csv <- file.path(BASE, "bakuzis_mortcap_BA.csv")
write.csv(res, csv, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(res), csv))

# Median BA by (config, year), the native-vs-cap comparison of interest
if (nrow(res)) {
  agg <- aggregate(BA ~ config + year, data = res, FUN = median)
  agg <- agg[order(agg$year, agg$config), ]
  cat("\nMedian BA (m2/ha) by config x horizon:\n")
  print(agg, row.names = FALSE)
  wide <- reshape(agg, idvar = "year", timevar = "config", direction = "wide")
  msum <- file.path(BASE, "bakuzis_mortcap_median_summary.csv")
  write.csv(wide, msum, row.names = FALSE)
  cat(sprintf("Wrote median summary to %s\n", msum))
}
