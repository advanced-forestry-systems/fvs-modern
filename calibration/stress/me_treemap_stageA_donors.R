## me_treemap_stageA_donors.R
## ME TreeMap pilot, Stage A: clip TreeMap2022 to Maine, extract the donor FIA
## plot set + per-plot pixel area. Output me_treemap_donors.csv with
##   TM_ID, PLT_CN, n_pix, area_ha
## so Stage B runs FVS on exactly those donor plots and Stage C paints/compares.
##
## TreeMap2022_CONUS.tif value = TM_ID; VAT maps TM_ID -> PLT_CN (= FIA plot CN
## = campaign STAND_CN). 30 m pixels -> 0.09 ha each.
##
## Usage: Rscript me_treemap_stageA_donors.R <tm_tif> <vat_dbf> <me_gpkg> <out_csv>

suppressWarnings(suppressMessages({library(terra); library(foreign)}))
a <- commandArgs(trailingOnly = TRUE)
tm_tif  <- a[1]; vat_dbf <- a[2]; me_gpkg <- a[3]; out_csv <- a[4]

r  <- rast(tm_tif)
me <- vect(me_gpkg)
me <- project(me, crs(r))                      # boundary -> raster CRS (Albers)

cat("cropping TreeMap to Maine extent...\n")
rc <- crop(r, me)                              # windowed read, ME extent only
rc <- mask(rc, me)                             # keep only inside ME

cat("tabulating donor plot pixel counts...\n")
levels(rc) <- NULL                             # strip RAT -> raw integer TM_ID
f <- freq(rc)                                  # value (TM_ID) -> count
f <- f[!is.na(f$value) & f$value > 0, c("value", "count")]
names(f) <- c("TM_ID", "n_pix")
f$TM_ID <- as.integer(f$TM_ID)

vat <- read.dbf(vat_dbf, as.is = TRUE)         # TM_ID -> PLT_CN
nm <- toupper(names(vat)); names(vat) <- nm
idcol  <- if ("TM_ID" %in% nm) "TM_ID" else if ("VALUE" %in% nm) "VALUE" else nm[1]
vat <- vat[, c(idcol, "PLT_CN")]
names(vat) <- c("TM_ID", "PLT_CN")
vat$TM_ID <- as.integer(vat$TM_ID)

d <- merge(f, vat, by = "TM_ID", all.x = TRUE)
d$area_ha <- d$n_pix * 0.09                     # 30 m pixel
d$PLT_CN  <- format(d$PLT_CN, scientific = FALSE, trim = TRUE)
d <- d[order(-d$n_pix), ]

write.csv(d, out_csv, row.names = FALSE)
cat(sprintf("ME donors: %d unique TM_ID, %d unique PLT_CN, %.0f forested ha, ",
            nrow(d), length(unique(d$PLT_CN)), sum(d$area_ha)))
cat(sprintf("%d total pixels -> %s\n", sum(d$n_pix), out_csv))
