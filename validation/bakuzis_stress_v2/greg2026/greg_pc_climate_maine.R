## Greg's ClimateNA PC inputs as 1991-2020 normals: DD_0 and DD_18 are degree-days BELOW 0 and 18 C
## (ClimateNA convention, as in Greg's climate_si script), TD = MWMT - MCMT, PPT_sm = Jun-Aug precip. From the
## Daymet daily cache built for the Maine FIA covariates (public coordinates, lapse-adjusted
## to FIA field elevation with the same monthly lapse rates).
suppressPackageStartupMessages(library(data.table))
B <- "/fs/scratch/PUOM0008/crsfaaron/maine_cov_2026/climate"
O <- "/users/PUOM0008/crsfaaron/fvs_remodeling_data/impact/climate"
lap <- fread(file.path(B, "climate_daymet/daymet_lapse_rates_degC_per_km.csv"))
pl <- fread(file.path(B, "ME_PLOT.csv"), select = c("CN","LAT","LON","ELEV","MEASYEAR","PLOT_STATUS_CD","STATECD","UNITCD","COUNTYCD","PLOT"))
pl <- pl[!is.na(LAT) & PLOT_STATUS_CD == 1]
pl[, loc_id := sprintf("%.5f_%.5f", LAT, LON)]
co <- fread(file.path(B, "ME_COND.csv"), select = c("PLT_CN","CONDID","FORTYPCD","CONDPROP_UNADJ","COND_STATUS_CD"))
co <- co[COND_STATUS_CD == 1][order(PLT_CN, -CONDPROP_UNADJ)][, .SD[1], by = PLT_CN]
pl <- merge(pl, co[, .(CN = PLT_CN, FORTYPCD)], by = "CN", all.x = TRUE)
pl <- pl[order(STATECD,UNITCD,COUNTYCD,PLOT,-MEASYEAR)][, .SD[1], by = .(STATECD,UNITCD,COUNTYCD,PLOT)]
mo <- as.integer(format(as.Date(0:364, origin = "2001-01-01"), "%m"))
res <- rbindlist(lapply(unique(pl$loc_id), function(id) {
  f <- file.path(B, "climate_daymet/cache", paste0(id, ".rds")); if (!file.exists(f)) return(NULL)
  d <- as.data.table(readRDS(f))[year >= 1991 & year <= 2020]
  d[, month := mo[pmin(yday, 365)]]
  data.table(loc_id = id, alt = d$altitude[1], d[, .(year, month, tmax, tmin, prcp)])[,
    .(tmax = tmax, tmin = tmin, prcp = prcp, year = year, month = month, alt = alt)][, loc_id := id][]
}), fill = TRUE)
el <- unique(pl[, .(loc_id, elev_m = ELEV * 0.3048)])[, .SD[1], by = loc_id]
res <- merge(res, el, by = "loc_id"); res <- merge(res, lap[, .(month, lapse_tmax, lapse_tmin)], by = "month")
res[, dz := (elev_m - alt) / 1000]
res[, tmean := ((tmax + lapse_tmax * dz) + (tmin + lapse_tmin * dz)) / 2]
yr <- res[, .(DD0 = sum(pmax(0 - tmean, 0)), DD_18 = sum(pmax(18 - tmean, 0)),
              PPT_sm = sum(prcp[month %in% 6:8])), by = .(loc_id, year)][,
          .(DD0 = mean(DD0), DD_18 = mean(DD_18), PPT_sm = mean(PPT_sm)), by = loc_id]
mt <- res[, .(t = mean(tmean)), by = .(loc_id, month)][, .(TD = max(t) - min(t)), by = loc_id]
cl <- merge(yr, mt, by = "loc_id")
cl[, PC1 := 0.252*DD0 + 0.002*TD - 0.035*PPT_sm + 0.967*DD_18]
cl[, PC2 := 0.882*DD0 + 0.015*TD + 0.420*PPT_sm - 0.215*DD_18]
pc <- merge(pl[, .(loc_id, FORTYPCD)], cl, by = "loc_id")
pc[, group := fcase(FORTYPCD %between% c(120,129), "Spruce-Fir", FORTYPCD %between% c(800,809), "Northern-Hardwood",
                    FORTYPCD %between% c(100,105), "Pine", FORTYPCD %between% c(400,409) | FORTYPCD %between% c(500,520), "Oak-Pine",
                    default = "Other")]
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
sm <- rbind(pc[, .(group = "All Maine forest", n = .N, stat = c("p10","median","p90"),
                   DD0 = q(DD0, c(.1,.5,.9)), TD = q(TD, c(.1,.5,.9)), PPT_sm = q(PPT_sm, c(.1,.5,.9)), DD_18 = q(DD_18, c(.1,.5,.9)))],
            pc[group != "Other", .(n = .N, stat = "median", DD0 = median(DD0), TD = median(TD), PPT_sm = median(PPT_sm), DD_18 = median(DD_18)), by = group])
sm[, PC1 := 0.252*DD0 + 0.002*TD - 0.035*PPT_sm + 0.967*DD_18][, PC2 := 0.882*DD0 + 0.015*TD + 0.420*PPT_sm - 0.215*DD_18]
fwrite(sm, file.path(O, "maine_greg_pc_climate_summary_1991-2020.csv"))
print(sm)
cat("DONE\n")
