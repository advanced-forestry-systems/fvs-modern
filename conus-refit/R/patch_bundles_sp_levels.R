## Quick patch: add sp_levels to residual bundles that lack them.
suppressPackageStartupMessages(library(data.table))

bundles <- list(
  list(bundle = "calibration/output/conus/dg_kue/v8/dg_kuehne_v8_smoke_residuals.rds",
       meta   = "calibration/output/conus/dg_kue/v8/dg_kuehne_v8_smoke_meta.rds"),
  list(bundle = "calibration/output/conus/cr/speciesfree/cr_sf_smoke_residuals.rds",
       meta   = "calibration/output/conus/cr/speciesfree/cr_sf_smoke_meta.rds"),
  list(bundle = "calibration/output/conus/hcb/speciesfree/hcb_sf_smoke_residuals.rds",
       meta   = "calibration/output/conus/hcb/speciesfree/hcb_sf_smoke_meta.rds")
)

for (b in bundles) {
  cat("Patching:", b$bundle, "\n")
  if (!file.exists(b$bundle) || !file.exists(b$meta)) {
    cat("  missing, skip\n"); next
  }
  bundle <- readRDS(b$bundle)
  meta <- readRDS(b$meta)
  bundle$sp_levels <- meta$sp_levels
  saveRDS(bundle, b$bundle)
  cat("  added sp_levels (n =", length(bundle$sp_levels), ")\n")
}
cat("\nDone.\n")
