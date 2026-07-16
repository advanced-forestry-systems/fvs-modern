#!/usr/bin/env Rscript
# 62j_greg_crown_csv.R
# Emit Greg's CONUS crown-change coefficients as a flat per-species CSV in the
# same GOMPLOAD-style format as 62i (header + SPCD n B0..). Companion to
# 62i_greg_dghg_csv.R and 62h_greg_mortality_csv.R; closes the crown gap that
# 62g still notes as "Greg repo has no fitted crown-change equation".
#
# Source: gregjohnsonbiometrics/fvs_remodeling commit 8310b2a (1 Jul 2026),
#   rds/delta_height_to_live_crown.RDS  (72 species; cols SPCD,CommonName,n,b0,b1,b2,b3,b4,rmse)
#   scripts/crown_change/Crown_Change_Equations_for_CONUS.qmd  (est_htlc)
#
# Equation (change in height to live crown over a period):
#   dHTLC = (CHT - CHTLC) * (1 - exp(B0 + B1*dHT + B2*dCCH))
#   CHT=total height, CHTLC=height to live crown, dHT=height growth,
#   dCCH=change in crown competition (CCH). b3,b4 are fit to 0 and dropped.
# New HTLC = CHTLC + dHTLC; CR = 1 - HTLC/HT.
suppressPackageStartupMessages(library(data.table))
RDS <- Sys.getenv("FVS_REMODELING_RDS", "/users/PUOM0008/crsfaaron/fvs_remodeling/rds")
OUT <- Sys.getenv("FVS_MODERN_CONFIG", "/users/PUOM0008/crsfaaron/fvs-modern/config")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cr <- as.data.table(readRDS(file.path(RDS, "delta_height_to_live_crown.RDS")))
cro <- cr[, .(SPCD = as.integer(SPCD), n = as.integer(n),
              B0 = b0, B1 = b1, B2 = b2)][order(SPCD)]
cro <- cro[is.finite(B0) & is.finite(B1) & is.finite(B2)]
fwrite(cro, file.path(OUT, "greg_crown_coefficients.csv"), quote = FALSE)
cat(sprintf("CROWN (delta HTLC): %d species -> %s/greg_crown_coefficients.csv\n",
            nrow(cro), OUT))
cat("cols: SPCD,n,B0,B1,B2  form: dHTLC=(CHT-CHTLC)*(1-exp(B0+B1*dHT+B2*dCCH))\n")

# ---- 62g follow-up (documented, not applied here) ---------------------------
# In 62g_greg_to_variant_json.R replace the crown stand-in:
#   crown = list(source = "fvs_conus_cr_recession",
#                note = "Greg repo has no fitted crown-change equation ...")
# with:
#   crown = list(source = "greg_delta_htlc",
#                coefficients = "greg_crown_coefficients.csv",
#                form = "dHTLC=(CHT-CHTLC)*(1-exp(B0+B1*dHT+B2*dCCH))",
#                note = "Greg CONUS crown-change (fvs_remodeling 8310b2a, 2026-07-01)")
# then re-run 62g to land categories_conus_greg.crown, and A/B vs the kernel.
# Landing is gated on the greg-arm draft branch rebuild + A/B.
