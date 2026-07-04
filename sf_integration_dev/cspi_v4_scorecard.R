## Add CSPI v4 to the site-productivity scorecard (ranger model). Predict v4 at
## the head-to-head plots (SICOND truth + cspi_v2/v3/csi + v3 covars), then
## correlate every layer against FIA SICOND.
suppressMessages({library(data.table); library(terra); library(ranger)})
C <- "/fs/scratch/PUOM0008/crsfaaron/cspi_v3"
d  <- as.data.table(readRDS(file.path(C,"head_to_head/head_to_head_dat.rds")))
m  <- readRDS(file.path(C,"cspi_v4_rf_model.rds"))
ver<- fread(file.path(C,"v4_climate_grid_verification.csv"))
need <- m$forest$independent.variable.names
cat("plots:", nrow(d), "| class:", paste(class(m),collapse=","), "| needs", length(need), "predictors\n")

r <- rast(file.path(C,"v4_climate_33_1km.tif"))
pts <- vect(as.data.frame(d[, .(LON_F, LAT_F)]), geom=c("LON_F","LAT_F"), crs="EPSG:4326")
pts <- project(pts, crs(r))
clim <- terra::extract(r, pts, ID=FALSE); setDT(clim)
for (cn in intersect(need, ver$CNA)) {
  bm <- ver[CNA==cn, best_match][1]
  if (!is.null(bm) && bm %in% names(clim)) { d[, (cn) := clim[[bm]]] }
  else if (cn %in% names(clim)) { d[, (cn) := clim[[cn]]] }
}
have <- intersect(need, names(d)); miss <- setdiff(need, names(d))
cat("have", length(have), "of", length(need), "predictors\n")
if (length(miss) > 0) cat("missing:", paste(head(miss,25), collapse=", "), "\n")
if (length(miss) == 0) {
  nd <- as.data.frame(d[, ..need])
  ok <- complete.cases(nd)
  d[, cspi_v4 := NA_real_]
  d[ok, cspi_v4 := predict(m, data=nd[ok,,drop=FALSE])$predictions]
  cat("predicted cspi_v4 for", sum(ok), "complete-case plots\n")
} else cat("SKIP v4 predict (missing predictors)\n")

co <- function(y,x){ok<-is.finite(x)&is.finite(y); if(sum(ok)<30) return(NA); round(cor(y[ok],x[ok]),3)}
r2 <- function(y,x){ok<-is.finite(x)&is.finite(y); if(sum(ok)<30) return(NA); round(summary(lm(y[ok]~x[ok]))$r.squared,3)}
layers <- intersect(c("csi","cspi_v2","cspi_v3","cspi_v4"), names(d))
out <- rbindlist(lapply(layers, function(L) data.table(layer=L,
        cor_SICOND=co(d$SICOND,d[[L]]), R2_SICOND=r2(d$SICOND,d[[L]]),
        n=sum(is.finite(d$SICOND)&is.finite(d[[L]])))))
cat("\n=== site-productivity scorecard vs FIA SICOND (head-to-head plots) ===\n"); print(out)
fwrite(out, file.path(C,"head_to_head/scorecard_vs_sicond_with_v4.csv"))
cat("\nDONE\n")
