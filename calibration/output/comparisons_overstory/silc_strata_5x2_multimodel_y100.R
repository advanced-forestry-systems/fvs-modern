#!/usr/bin/env Rscript
# silc_strata_5x2_multimodel_y100.R
# Year-100 multi-model comparison (AGM def vs AGM mc vs FVS-NE def vs FVS-NE cal)
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"

def <- read.csv(file.path(od, "silc_strata_5x2_AGM_trajectories.csv"), stringsAsFactors=FALSE)
mc  <- read.csv(file.path(od, "silc_strata_5x2_AGM_MORTCAL_trajectories.csv"), stringsAsFactors=FALSE)
fvs <- read.csv(file.path(od, "silc_strata_100yr_fvsne_results.csv"), stringsAsFactors=FALSE)
sm  <- read.csv(file.path(od, "silc_strata_5x2_mapping.csv"), stringsAsFactors=FALSE)

# Year-100 (Year=2123)
def100 <- def[def$Year == 2123, c("forest_type","density_class","BA","NetCords")]
names(def100)[3:4] <- c("BA_agm_def","Cords_agm_def")

mc100 <- mc[mc$Year == 2123, c("forest_type","density_class","BA","NetCords")]
names(mc100)[3:4] <- c("BA_agm_mc","Cords_agm_mc")

fvs100 <- fvs[fvs$Year == 2123, ]
fvs100$Cords <- fvs100$MCuFt / 79
fvs100 <- merge(fvs100, sm[, c("stand_id","forest_type","density_class")], by.x="StandID", by.y="stand_id")
fvs_agg <- aggregate(cbind(BA, Cords) ~ forest_type + density_class + config, data=fvs100, FUN=mean)
fvs_def <- fvs_agg[fvs_agg$config=="default", c("forest_type","density_class","BA","Cords")]
names(fvs_def)[3:4] <- c("BA_fvs_def","Cords_fvs_def")
fvs_cal <- fvs_agg[fvs_agg$config=="calibrated", c("forest_type","density_class","BA","Cords")]
names(fvs_cal)[3:4] <- c("BA_fvs_cal","Cords_fvs_cal")

y100 <- Reduce(function(a, b) merge(a, b, by=c("forest_type","density_class"), all=TRUE),
               list(def100, mc100, fvs_def, fvs_cal))
y100 <- y100[order(y100$forest_type, y100$density_class), ]
print(y100, row.names=FALSE, digits=3)
write.csv(y100, file.path(od, "silc_strata_5x2_year100_multimodel.csv"), row.names=FALSE)
