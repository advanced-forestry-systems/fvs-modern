# Export population SDIMAX posterior draws (exp(Intercept) at QMD=10in) as a
# CSV sidecar next to stand_density_samples.rds for the Python density cap.
args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else
  "/users/PUOM0008/crsfaaron/fvs-conus/output/variants/ne/stand_density_samples.rds"
s <- as.data.frame(readRDS(rds))
col <- if ("Intercept" %in% names(s)) "Intercept" else "b_Intercept"
sdimax <- exp(s[[col]])
sdimax <- sdimax[is.finite(sdimax) & sdimax > 0]
out <- paste0(rds, ".sdimax.csv")
write.table(sdimax, out, row.names = FALSE, col.names = FALSE)
cat("wrote", length(sdimax), "SDIMAX draws to", out, "\n")
cat("median SDIMAX:", round(median(sdimax), 1),
    " q05:", round(quantile(sdimax, .05), 1),
    " q95:", round(quantile(sdimax, .95), 1), "\n")
