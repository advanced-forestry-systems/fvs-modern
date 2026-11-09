suppressPackageStartupMessages(library(data.table))
O <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/bakuzis_greg2026_20261109"
rd <- function(t) fread(sprintf("%s/engine_%s/engine_bakuzis_100yr_%s.csv", O, t, t))
map <- c("Spruce-Fir"="gregME_sf","Northern-Hardwood"="gregME_nh","Pine"="gregME_pine","Oak-Pine"="gregME_op")
cmp <- rbindlist(lapply(names(map), function(g) rd(map[[g]])[species_group == g]))
cmp[, tag := "gregME_grp"]
dir.create(file.path(O, "engine_gregME_grp"), showWarnings = FALSE)
fwrite(cmp, file.path(O, "engine_gregME_grp/engine_bakuzis_100yr_gregME_grp.csv"))
cat("rows", nrow(cmp), "\n")
