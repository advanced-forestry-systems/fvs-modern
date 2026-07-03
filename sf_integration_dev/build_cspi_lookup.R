## Build a PLT_CN -> CSPI v4 lookup at the calibration plots.
## Coords: de-fuzzed locator candidate (rank 1) where PLT_CN_cond1 == PLOT_UID,
## else the calibration public LAT/LON. Predictors: 33 climate (v4 stack) + 11
## soil/terrain/canopy (aligned_1km), point-extracted at each plot.
suppressMessages({library(data.table); library(bit64); library(terra); library(ranger)})
CAL <- "/users/PUOM0008/crsfaaron/fvs-modern/calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"
LOC <- "/fs/scratch/PUOM0008/crsfaaron/fia_locator/national/conus_ranked_plot_candidates.csv"
C   <- "/fs/scratch/PUOM0008/crsfaaron/cspi_v3"
A   <- file.path(C,"aligned_1km")
OUT <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

dat <- readRDS(CAL); setDT(dat)
plots <- unique(dat[, .(PLT_CN_cond1, LAT, LON)])
plots <- plots[is.finite(LAT) & is.finite(LON)]
cat("unique calibration plots with coords:", nrow(plots), "\n")

## de-fuzz via locator rank1
loc <- fread(LOC, select=c("PLOT_UID","candidate_lat","candidate_lon","rank"))
loc <- loc[rank==1, .(PLT_CN_cond1=PLOT_UID, dlat=candidate_lat, dlon=candidate_lon)]
plots <- merge(plots, loc, by="PLT_CN_cond1", all.x=TRUE)
nmatch <- plots[is.finite(dlat), .N]
cat("locator de-fuzz matched:", nmatch, "of", nrow(plots), sprintf("(%.1f%%)\n", 100*nmatch/nrow(plots)))
plots[, use_lat := fifelse(is.finite(dlat), dlat, LAT)]
plots[, use_lon := fifelse(is.finite(dlon), dlon, LON)]

## ---- extract predictors at coords ----
m   <- readRDS(file.path(C,"cspi_v4_rf_model.rds"))
ver <- fread(file.path(C,"v4_climate_grid_verification.csv"))
need <- m$forest$independent.variable.names
pts4326 <- vect(as.data.frame(plots[, .(use_lon, use_lat)]), geom=c("use_lon","use_lat"), crs="EPSG:4326")

# climate (33) from v4 stack
rc <- rast(file.path(C,"v4_climate_33_1km.tif"))
clim <- terra::extract(rc, project(pts4326, crs(rc)), ID=FALSE); setDT(clim)
for (cn in intersect(need, ver$CNA)) {
  bm <- ver[CNA==cn, best_match][1]
  if (!is.null(bm) && bm %in% names(clim)) plots[, (cn) := clim[[bm]]]
  else if (cn %in% names(clim)) plots[, (cn) := clim[[cn]]]
}
# soil/terrain/canopy (11) from aligned_1km; map h_loss_raw -> h_loss_10y
afiles <- c(h_tc2000="h_tc2000", h_loss_10y="h_loss_raw", sand_0_5="sand_0_5",
            soc_0_5="soc_0_5", bdod_0_5="bdod_0_5", cec_0_5="cec_0_5",
            nitrogen_0_5="nitrogen_0_5", phh2o_0_5="phh2o_0_5",
            elev="elev", slope="slope", aspect="aspect")
for (nm in names(afiles)) {
  f <- file.path(A, paste0(afiles[[nm]], ".tif"))
  if (file.exists(f)) {
    rr <- rast(f); v <- terra::extract(rr, project(pts4326, crs(rr)), ID=FALSE)
    plots[, (nm) := v[[1]]]
  }
}
have <- intersect(need, names(plots)); miss <- setdiff(need, names(plots))
cat("predictors have", length(have), "of", length(need), "; missing:", paste(head(miss,20),collapse=", "), "\n")

if (length(miss)==0) {
  nd <- as.data.frame(plots[, ..need]); ok <- complete.cases(nd)
  plots[, cspi_v4 := NA_real_]
  plots[ok, cspi_v4 := predict(m, data=nd[ok,,drop=FALSE])$predictions]
  cat("predicted cspi_v4 for", sum(ok), "plots\n")
  cat("cspi_v4 summary:\n"); print(summary(plots$cspi_v4))
  out <- plots[, .(PLT_CN_cond1=as.character(PLT_CN_cond1), cspi_v4,
                   use_lat, use_lon, defuzzed=is.finite(dlat))]
  fwrite(out, file.path(OUT,"cspi_v4_at_calib_plots.csv"))
  cat("wrote", file.path(OUT,"cspi_v4_at_calib_plots.csv"), "\n")
} else cat("SKIP predict (missing predictors)\n")
cat("DONE\n")
