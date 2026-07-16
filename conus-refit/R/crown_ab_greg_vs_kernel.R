#!/usr/bin/env Rscript
# crown_ab_greg_vs_kernel.R
# A/B: Greg CONUS crown-change (est_htlc) vs the fvs-conus CR-recession kernel
# (Hann & Hanus, script 35) at predicting OBSERVED delta-HCB on held-out FIA
# remeasurement pairs. Decision input for the greg-arm crown swap (P9 / PR #90).
#
# Faithfulness:
#  - Kernel predictors are built by the REAL prep function (prepare_cr_data,
#    sourced from 35_fit_crown_recession.R) so transforms match the fit exactly.
#  - Kernel eta uses fitted r0..r6 + per-species RE; spatial L1/L2/L3 REs set to
#    their mean 0 (marginal prediction). Greg has no spatial RE, so this is a
#    fair like-for-like marginal comparison.
#  - Greg unit convention (fvs_remodeling may be imperial; our pairs are metric)
#    is RESOLVED EMPIRICALLY: Greg is run with height inputs in meters AND in
#    feet; the convention whose predicted delta-HCB magnitude matches observed
#    is the intended one. A dCCH=0 arm bounds competition-term sensitivity.
# Author: Cowork autopilot, 2026-07-03. Seed 20260703.

suppressPackageStartupMessages({ library(data.table) })
set.seed(20260703)
M2FT <- 3.28084
OUT  <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/crown_ab"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FS   <- "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus"

# ---- kernel prep (reuse real code; guard stops the fitting driver) ----------
Sys.setenv(FVS_CONUS_SKIP_DRIVER = "1")
source("R/35_fit_crown_recession.R")     # defines prepare_cr_data()

d <- readRDS("calibration/data/conus_remeasurement_pairs.rds"); setDT(d)
# rd_add is a relative-density-to-SDImax ratio in (0,2) (prep filters on that range).
if (!"rd_add" %in% names(d)) d[, rd_add := if ("rd_sdimax" %in% names(d)) rd_sdimax else RD1]
# Pairs carry only EPA_L3_CODE. Spatial L1/L2/L3 REs are marginalized to 0 in the
# kernel prediction below, so L1/L2 need only be non-NA to satisfy the prep filter.
if (!"EPA_L1_CODE" %in% names(d)) d[, EPA_L1_CODE := EPA_L3_CODE]
if (!"EPA_L2_CODE" %in% names(d)) d[, EPA_L2_CODE := EPA_L3_CODE]
# held-out sample: kernel was fit on ~401k of 8.22M; a fresh 400k is ~95% unseen
if (nrow(d) > 400000L) d <- d[sample(.N, 400000L)]

d[, SPCD_orig := SPCD]   # preserve pre-pooling species for a fair Greg join
prep <- prepare_cr_data(d, site_var = "cspi", min_sp_n = 5000L)
dd <- prep$data
cat("prepped obs:", nrow(dd), "\n")

# ---- kernel prediction -------------------------------------------------------
r  <- fread(file.path(FS, "cr_recession_cspi_traits1_fixed_summary.csv"))
rv <- setNames(r$mean, r$variable)
zs <- fread(file.path(FS, "cr_recession_cspi_traits1_species_intercepts.csv"))  # variable,mean,...,SPCD
dd <- merge(dd, zs[, .(SPCD, z_sp = mean)], by = "SPCD", all.x = TRUE)
dd[is.na(z_sp), z_sp := 0]
dd[, eta := rv["r0"] + z_sp +
      rv["r1"]*ln_cr + rv["r2"]*CR1 + rv["r3"]*sqrt_ba +
      rv["r4"]*ln_bal_ba + rv["r5"]*cr_over_rd + rv["r6"]*ln_cspi_shift]
dd[, dhcb_kernel := max_dhcb / (1 + exp(eta))]        # (CL_S+dH) * inv_logit(-eta)

# ---- Greg prediction (both unit conventions + dCCH=0 sensitivity) -----------
g <- fread(file.path("/users/PUOM0008/crsfaaron/fvs-modern/config",
                     "greg_crown_coefficients.csv"))   # SPCD,n,B0,B1,B2
medB <- g[, .(B0 = median(B0), B1 = median(B1), B2 = median(B2))]
dd <- merge(dd, g[, .(SPCD_orig = SPCD, gB0 = B0, gB1 = B1, gB2 = B2)], by = "SPCD_orig", all.x = TRUE)
matched <- dd[!is.na(gB0), .N] / nrow(dd)
dd[is.na(gB0), `:=`(gB0 = medB$B0, gB1 = medB$B1, gB2 = medB$B2)]
dd[, dcch := CCH2 - CCH1]
prop <- function(b0,b1,b2,dht,dcch) pmin(pmax(1 - exp(b0 + b1*dht + b2*dcch), 0), 1)
dd[, dhcb_greg_m   := cl_s * prop(gB0, gB1, gB2, delta_h,         dcch)]  # dHT in meters
dd[, dhcb_greg_ft  := cl_s * prop(gB0, gB1, gB2, delta_h*M2FT,    dcch)]  # dHT in feet
dd[, dhcb_greg_ft0 := cl_s * prop(gB0, gB1, gB2, delta_h*M2FT,    0)]     # dCCH=0 sensitivity

# ---- metrics -----------------------------------------------------------------
obs <- dd$delta_hcb
met <- function(p) c(RMSE = sqrt(mean((p-obs)^2)), MBE = mean(p-obs),
                     MAE = mean(abs(p-obs)), r = suppressWarnings(cor(p, obs)))
res <- rbind(
  kernel        = met(dd$dhcb_kernel),
  greg_meters   = met(dd$dhcb_greg_m),
  greg_feet     = met(dd$dhcb_greg_ft),
  greg_feet_nocch = met(dd$dhcb_greg_ft0))
res <- data.table(model = rownames(res), res)
cat("\n=== overall (observed dHCB, meters) ===\n"); print(res)
cat(sprintf("\nmean observed dHCB: %.3f m | greg species matched: %.1f%%\n",
            mean(obs), 100*matched))

# pick the Greg unit convention with smaller overall |bias| (empirical unit resolution)
greg_col <- if (abs(mean(dd$dhcb_greg_ft - obs)) <= abs(mean(dd$dhcb_greg_m - obs)))
              "dhcb_greg_ft" else "dhcb_greg_m"
cat("\nchosen Greg unit convention (min |bias|):", greg_col, "\n")
# by softwood/hardwood (SPCD < 300 ~ softwood in FIA)
dd[, grp := ifelse(SPCD > 0 & SPCD < 300, "softwood", "hardwood")]
dd[, greg_best := get(greg_col)]
by_grp <- dd[, .(n=.N, obs=mean(delta_hcb),
                 kernel_RMSE=sqrt(mean((dhcb_kernel-delta_hcb)^2)),
                 greg_RMSE=sqrt(mean((greg_best-delta_hcb)^2))),
             by = grp]
cat("\n=== RMSE by group ===\n"); print(by_grp)

fwrite(res, file.path(OUT, "crown_ab_overall_metrics.csv"))
fwrite(by_grp, file.path(OUT, "crown_ab_by_group.csv"))
saveRDS(dd[, .(SPCD, grp, delta_hcb, dhcb_kernel, dhcb_greg_m, dhcb_greg_ft)],
        file.path(OUT, "crown_ab_predictions.rds"))

# obs-vs-pred figure (thumbnail-safe)
png(file.path(OUT, "crown_ab_obs_vs_pred.png"), width = 1200, height = 500, res = 110)
op <- par(mfrow = c(1,2))
for (nm in c("dhcb_kernel","dhcb_greg_ft")) {
  s <- dd[sample(.N, min(.N, 20000))]
  plot(s$delta_hcb, s[[nm]], pch=".", col="#3355aa55",
       xlab="observed dHCB (m)", ylab=paste("pred", nm),
       main=nm, xlim=c(0,15), ylim=c(0,15)); abline(0,1,col="red")
}
par(op); dev.off()
cat("\nwrote:", OUT, "\n")
