# profile the conus remeasurement data for forest-type and ecoregion columns
fs <- c(Sys.glob("~/fvs-conus/calibration/data/conus_remeasurement_pairs*.rds"),
        Sys.glob("~/fvs-conus/calibration/data/*remeas*.rds"),
        Sys.glob("~/fvs-conus/data/*remeas*.rds"))
fs <- unique(fs[file.exists(fs)])
cat("candidate data files:\n"); print(fs)
if(length(fs)){
  d <- readRDS(fs[1]); cat("\nusing:", fs[1], " rows=", nrow(d), "\n")
  hits <- grep("FORTYP|FLDTYP|FORTYPCD|FORTYPGRP|TYPE|ECO|L1|L2|L3|^FT|forest|eco", names(d),
               ignore.case=TRUE, value=TRUE)
  cat("forest-type / ecoregion columns present:\n"); print(hits)
  for(cc in intersect(c("FORTYPCD","FLDTYPCD","FORTYPGRP","FORTYPGCD"), names(d)))
    cat(cc, ": n_levels=", length(unique(d[[cc]])), " example=", paste(head(unique(d[[cc]]),6),collapse=","), "\n")
}
cat("DONE_PROFILE\n")
