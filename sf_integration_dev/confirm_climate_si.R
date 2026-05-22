## Confirm whether the calibration data climate_si is the broken climate-RF (csi)
## by correlating it with FIA SICOND and with cspi_v2/v3 at the same plots.
suppressMessages(library(data.table))
dat <- as.data.table(readRDS("calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"))
nm <- names(dat)
cat("candidate cols:", paste(grep("si$|_si|cspi|csi|PLT|^CN$|PLOT", nm, value=TRUE, ignore.case=TRUE), collapse=", "), "\n")
pid <- intersect(c("PLT_CN","PLT_CN1","PLOT_CN","CN","PLOT"), nm)[1]
sivar <- intersect(c("climate_si","CSPI","csi"), nm)[1]
cat("using plot id:", pid, " | productivity col:", sivar, "\n")
if (is.na(pid) || is.na(sivar)) { cat("MISSING pid or sivar; abort\n"); quit(save="no", status=0) }
d1 <- unique(dat[, c(pid, sivar), with=FALSE]); setnames(d1, c("PLT_CN","prod"))
d1[, PLT_CN := format(PLT_CN, scientific=FALSE, trim=TRUE)]
v2 <- fread("/fs/scratch/PUOM0008/crsfaaron/cspi_v2/cspi_v2_at_plots.csv", colClasses=list(character="PLT_CN"))
m <- merge(d1, v2[, .(PLT_CN, SICOND, cspi_v2)], by="PLT_CN")
m <- m[is.finite(SICOND) & SICOND>0 & is.finite(prod)]
cat("joined plots:", nrow(m), "\n")
co <- function(y,x){ok<-is.finite(x)&is.finite(y); round(cor(y[ok],x[ok]),3)}
cat(sprintf("cor(calib %s, FIA SICOND) = %.3f\n", sivar, co(m$SICOND, m$prod)))
cat(sprintf("cor(calib %s, cspi_v2)    = %.3f\n", sivar, co(m$cspi_v2, m$prod)))
cat(sprintf("cor(cspi_v2, FIA SICOND)  = %.3f\n", co(m$SICOND, m$cspi_v2)))
cat(sprintf("range calib %s: %s\n", sivar, paste(round(range(m$prod),2), collapse=" .. ")))
fwrite(m, "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration/climate_si_check.csv")
cat("DONE\n")
