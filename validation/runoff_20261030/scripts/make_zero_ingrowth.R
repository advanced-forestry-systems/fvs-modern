## Zero-ingrowth lookup for the fvs-conus projector (contestant B): same structure as the banked
## fvs-conus output/comparisons_overstory/intermediate/ingrowth_lookup.rds, every rate set to 0, so the
## projector's FIX 2 recruitment is off, matching the engines' no-regeneration protocol.
src <- "/users/PUOM0008/crsfaaron/fvs-conus/output/comparisons_overstory/intermediate/ingrowth_lookup.rds"
x <- readRDS(src); z <- lapply(x, function(e) { if (is.list(e)) { e$med_ann_TPA <- 0; e$med_ann_BA <- 0 }; e })
names(z) <- names(x); saveRDS(z, commandArgs(TRUE)[1])
cat("entries:", length(z), " names:", paste(head(names(z), 30), collapse = ","), "\n")
cat("source rates (rec/ac/yr) by entry:\n"); print(sapply(x, function(e) if (is.list(e)) e$med_ann_TPA else NA))
