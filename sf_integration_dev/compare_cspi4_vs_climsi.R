## HT-DBH: does CSPI v4 beat climate_si as the productivity covariate?
suppressMessages({library(data.table); library(loo)})
B <- "/users/PUOM0008/crsfaaron/fvs-conus/output/conus/htdbh"
climsi_loo <- Sys.glob(file.path(B,"unified_v2split","*_loo.rds"))
cspi4_loo  <- Sys.glob(file.path(B,"unified_cspi4_v2split","*_loo.rds"))
climsi_sum <- Sys.glob(file.path(B,"unified_v2split","*_summary.csv"))
cspi4_sum  <- Sys.glob(file.path(B,"unified_cspi4_v2split","*_summary.csv"))
cat("found climsi_loo:", length(climsi_loo), " cspi4_loo:", length(cspi4_loo), "\n")
if (length(climsi_loo) && length(cspi4_loo)) {
  a <- readRDS(climsi_loo[1]); b <- readRDS(cspi4_loo[1])
  cat("N obs: climate_si=", dim(a$pointwise)[1], " cspi4=", dim(b$pointwise)[1], "\n")
  cat("\n=== LOO compare (positive favors top row) ===\n")
  print(loo_compare(list(cspi_v4=b, climate_si=a)))
}
ac <- function(f){ s<-fread(f); s[variable=="a_cspi", .(mean,q5,q95)] }
if (length(climsi_sum) && length(cspi4_sum)) {
  cat("\n=== a_cspi (productivity coefficient) ===\n")
  cat("climate_si:", paste(unlist(ac(climsi_sum[1])),collapse=" "), "\n")
  cat("cspi_v4   :", paste(unlist(ac(cspi4_sum[1])),collapse=" "), "\n")
}
cat("DONE\n")
