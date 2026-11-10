suppressPackageStartupMessages(library(data.table))
O <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/bakuzis_greg2026_20261109"
rd <- function(t) fread(sprintf("%s/engine_%s/engine_bakuzis_100yr_%s.csv", O, t, t))
P <- Sys.getenv("PFX", "gregME")
map <- c("Spruce-Fir"=paste0(P,"_sf"),"Northern-Hardwood"=paste0(P,"_nh"),"Pine"=paste0(P,"_pine"),"Oak-Pine"=paste0(P,"_op"))
cmp <- rbindlist(lapply(names(map), function(g) rd(map[[g]])[species_group == g]))
cmp[, tag := paste0(P,"_grp")]
dir.create(file.path(O, paste0("engine_",P,"_grp")), showWarnings = FALSE)
fwrite(cmp, file.path(O, sprintf("engine_%s_grp/engine_bakuzis_100yr_%s_grp.csv",P,P)))
cat("rows", nrow(cmp), "\n")
