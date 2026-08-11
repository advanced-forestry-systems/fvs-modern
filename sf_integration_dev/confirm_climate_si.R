## Within calibration data: compare climate_si (used by models) vs the unused
## cspi column, and against any site-index truth column present. Also diagnose
## the key formats for the external join.
suppressMessages(library(data.table))
dat <- as.data.table(readRDS("calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds"))
nm <- names(dat)
cat("=== key format diagnosis ===\n")
cat("PLT_CN_cond1 head:", paste(head(as.character(dat$PLT_CN_cond1),4),collapse=" | "), "\n")
v2 <- fread("/fs/scratch/PUOM0008/crsfaaron/cspi_v2/cspi_v2_at_plots.csv", colClasses=list(character="PLT_CN"))
cat("v2 PLT_CN head:", paste(head(v2$PLT_CN,4),collapse=" | "), "\n\n")

prods <- intersect(c("climate_si","cspi","CSPI","csi"), nm)
cat("=== productivity cols present:", paste(prods,collapse=", "), "===\n")
for (p in prods) {
  x <- dat[[p]]; cat(sprintf("  %-12s NA%%=%.1f  range=[%s]  sd=%.3f\n", p,
    100*mean(!is.finite(x)), paste(round(range(x[is.finite(x)]),2),collapse=", "), sd(x[is.finite(x)])))
}
if (all(c("climate_si","cspi") %in% nm)) {
  ok <- is.finite(dat$climate_si) & is.finite(dat$cspi)
  cat(sprintf("\n  cor(climate_si, cspi) = %.3f  (n=%d)\n", cor(dat$climate_si[ok], dat$cspi[ok]), sum(ok)))
}
## look for any site index truth column
sicols <- grep("SICOND|site_index|^SI_|_SI$|SI_m|SDImax|SICOND_cond", nm, value=TRUE, ignore.case=TRUE)
cat("\n=== candidate site-index truth cols:", paste(head(sicols,10),collapse=", "), "===\n")
for (s in head(sicols,6)) {
  for (p in prods) {
    ok <- is.finite(dat[[s]]) & is.finite(dat[[p]])
    if (sum(ok) > 50) cat(sprintf("  cor(%s, %s) = %.3f (n=%d)\n", p, s, cor(dat[[p]][ok], dat[[s]][ok]), sum(ok)))
  }
}
cat("DONE\n")
