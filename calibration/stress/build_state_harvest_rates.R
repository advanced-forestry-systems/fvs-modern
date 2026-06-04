#!/usr/bin/env Rscript
# build_state_harvest_rates.R
# Per-state annual harvest + disturbance rates for the FVS "managed (harvest)"
# PERSEUS scenario, sampled at the actual FIA plot locations so each state's
# rate is the harvest/disturbance regime its own plots sit in (data-driven,
# region/forest-type specific via plot geography, repeatable).
#
# Harvest driver: conus_expected_ba_removed_annual (occurrence x intensity,
#   the expected FRACTION of standing basal area removed per year), NAD83 Albers.
#   Also samples p_partial / p_clearcut / p_stand_replacement for context.
# Disturbance driver: p_disturbance_2022.tif (annual P(disturbance)), WGS84.
#
# Output: state_harvest_rates.csv  (state, n, harvest_frac_yr, p_partial,
#   p_clearcut, p_standrepl, disturbance_frac_yr)
suppressMessages({library(terra); library(data.table)})

SD   <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
MAPS <- "/fs/scratch/PUOM0008/crsfaaron/conus_hcs/data/analytic/maps_conus_canonical"
DIST <- "/fs/scratch/PUOM0008/crsfaaron/TREEMAP_outputs_wgs84/p_disturbance_2022.tif"
OUT  <- file.path(SD, "state_harvest_rates.csv")

# ---- plot points (lat/lon/state) from all standinit files ----
sis <- list.files(file.path(SD, "standinit_by_variant"),
                  pattern = "^standinit_.*\\.csv$", full.names = TRUE)
pts <- rbindlist(lapply(sis, function(f)
  fread(f, select = c("LATITUDE", "LONGITUDE", "STATE"))), fill = TRUE)
pts <- pts[is.finite(LATITUDE) & is.finite(LONGITUDE) &
           LONGITUDE < -50 & LATITUDE > 17]
# cap per state for speed (rate is a mean; 8000 pts/state is ample)
set.seed(1)
pts <- pts[, .SD[sample(.N, min(.N, 8000))], by = STATE]
cat(sprintf("%d plot points across %d states\n", nrow(pts), uniqueN(pts$STATE)))

FIPS <- c("1"="AL","2"="AK","4"="AZ","5"="AR","6"="CA","8"="CO","9"="CT","10"="DE",
 "12"="FL","13"="GA","16"="ID","17"="IL","18"="IN","19"="IA","20"="KS","21"="KY",
 "22"="LA","23"="ME","24"="MD","25"="MA","26"="MI","27"="MN","28"="MS","29"="MO",
 "30"="MT","31"="NE","32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC",
 "38"="ND","39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD",
 "47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV","55"="WI","56"="WY")
pts[, ST := FIPS[as.character(as.integer(STATE))]]
pts <- pts[!is.na(ST)]

vpt <- vect(as.data.frame(pts), geom = c("LONGITUDE", "LATITUDE"),
            crs = "EPSG:4326")

samp <- function(path, albers) {
  r <- rast(path)
  v <- if (albers) project(vpt, crs(r)) else vpt
  as.numeric(terra::extract(r, v, ID = FALSE)[, 1])
}

cat("sampling harvest rasters (Albers)...\n")
pts[, harvest_frac := samp(file.path(MAPS,
  "conus_expected_ba_removed_annual_240m_2024_conus.tif"), TRUE)]
pts[, p_partial   := samp(file.path(MAPS,
  "conus_p_partial_annual_240m_2024_conus.tif"), TRUE)]
pts[, p_clearcut  := samp(file.path(MAPS,
  "conus_p_clearcut_annual_240m_2024_conus.tif"), TRUE)]
pts[, p_standrepl := samp(file.path(MAPS,
  "conus_p_stand_replacement_annual_240m_2024_conus.tif"), TRUE)]
cat("sampling disturbance raster (WGS84)...\n")
pts[, dist_frac   := samp(DIST, FALSE)]

agg <- pts[, .(
  n                  = .N,
  harvest_frac_yr    = mean(harvest_frac, na.rm = TRUE),
  p_partial          = mean(p_partial,    na.rm = TRUE),
  p_clearcut         = mean(p_clearcut,   na.rm = TRUE),
  p_standrepl        = mean(p_standrepl,  na.rm = TRUE),
  disturbance_frac_yr= mean(dist_frac,    na.rm = TRUE)
), by = ST][order(ST)]
fwrite(agg, OUT)
cat("wrote", OUT, "\n")
print(agg, nrows = 60)
