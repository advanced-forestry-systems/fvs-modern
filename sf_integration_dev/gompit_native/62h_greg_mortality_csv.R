#!/usr/bin/env Rscript
# 62h_greg_mortality_csv.R
# Emit greg_mortality_coefficients.csv in the exact format base/gompmort.f90
# (GOMPLOAD) reads: one header line, then rows of  SPCD n b0 b1 b2 b3 b4
# (list-directed, comma-separated OK). FIA SPCD is mapped onto FVS species by the
# engine via FIAJSP. This is the one missing file that lets Greg's cch-gompit
# mortality (eta = b0 + b1*(cr+0.01)^b2 + b3*cch^b4) run natively inside the FVS
# growth loop, activated by FVS_GOMPIT=1 + FVS_GOMPIT_COEF=<this file>.
suppressPackageStartupMessages(library(data.table))
ga <- function(n,d){ m<-grep(paste0("^--",n,"="),commandArgs(TRUE),value=TRUE); if(length(m)) sub(paste0("^--",n,"="),"",m[1]) else d }
SRC <- ga("src", "/users/PUOM0008/crsfaaron/fvs_remodeling/rds/mort_parm_base_rate_cr_cch.RDS")
OUT <- ga("out", "/users/PUOM0008/crsfaaron/fvs-modern/config/greg_mortality_coefficients.csv")
mo <- as.data.table(readRDS(SRC))
stopifnot(all(c("SPCD","n","b0","b1","b2","b3","b4") %in% names(mo)))
mo[, SPCD := as.integer(SPCD)]
# de-duplicate: Greg's RDS carries multiple fits for some species; keep the best
# (lowest nll) single row per SPCD so the engine's FIAJSP map is unambiguous.
has_nll <- "nll" %in% names(mo)
if (has_nll) setorder(mo, SPCD, nll) else setorder(mo, SPCD)
dd <- mo[, .SD[1], by = SPCD]
out <- dd[, .(SPCD, n = as.integer(n), b0, b1, b2, b3, b4)][order(SPCD)]
ok <- out[is.finite(b0) & is.finite(b1) & is.finite(b2) & is.finite(b3) & is.finite(b4)]
fwrite(ok, OUT, quote = FALSE)
cat(sprintf("wrote %s : %d unique species (from %d rows); dedup by %s\n",
            OUT, nrow(ok), nrow(mo), if (has_nll) "min nll" else "first"))
cat("activate: export FVS_GOMPIT=1; export FVS_GOMPIT_COEF=", OUT, "\n", sep = "")
