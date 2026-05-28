#!/usr/bin/env python3
# Build cardinal_acadgy_siteindex_v27b.R = v24 (12.3.8 ingrowth-fix FIA
# baseline) but compare three site productivity metrics already in vdat:
#   csi_baseline    ClimateSI_ft as CSI (matches v24, sanity)
#   sicond_replace  SICOND (FIA condition site index in feet)
#   fvssi_replace   FVS_SITE_INDEX (FVS-translated site index in feet)
# All three are read directly from validation_data_acd_post.csv. BGI from the
# ME_BGI_V1.tif raster is a follow-on (needs plot lat/lon and raster::extract).
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_ingrowthfix_v24.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_siteindex_v27b.R"
s = open(src).read()

s = s.replace("[v24]", "[v27b]")
s = s.replace("acadgy_ingrowthfix_v24_results.csv", "acadgy_siteindex_v27b_results.csv")

# Inject SICOND / FVS_SITE_INDEX lookups before run block
anchor = "stand_init <- unique(data.frame(STAND = samp$PLT_CN_t1_str, CSI = samp$ClimateSI_m, ELEV = samp$ELEV_m, stringsAsFactors = FALSE))"
assert anchor in s
inject = (
    anchor + "\n"
    "# v27b: parallel stand_init dataframes for each site productivity metric.\n"
    "# CSI is already in meters (ClimateSI_ft * 0.3048). SICOND and FVS_SITE_INDEX\n"
    "# are in feet in vdat; convert to meters for parity with the CSI variable\n"
    "# that AcadianGYOneStand consumes.\n"
    "stand_init_csi  <- stand_init\n"
    "stand_init_sicond <- data.frame(\n"
    "  STAND = samp$PLT_CN_t1_str,\n"
    "  CSI   = ifelse(!is.na(samp$SICOND) & samp$SICOND > 0, samp$SICOND * 0.3048, samp$ClimateSI_m),\n"
    "  ELEV  = samp$ELEV_m, stringsAsFactors=FALSE)\n"
    "stand_init_sicond <- unique(stand_init_sicond)\n"
    "stand_init_fvssi <- data.frame(\n"
    "  STAND = samp$PLT_CN_t1_str,\n"
    "  CSI   = ifelse(!is.na(samp$FVS_SITE_INDEX) & samp$FVS_SITE_INDEX > 0, samp$FVS_SITE_INDEX * 0.3048, samp$ClimateSI_m),\n"
    "  ELEV  = samp$ELEV_m, stringsAsFactors=FALSE)\n"
    "stand_init_fvssi <- unique(stand_init_fvssi)\n"
    "cat(sprintf(\"[v27b] CSI mean=%.2f m | SICOND mean=%.2f m | FVS_SITE_INDEX mean=%.2f m\\n\",\n"
    "    mean(stand_init_csi$CSI, na.rm=TRUE),\n"
    "    mean(stand_init_sicond$CSI, na.rm=TRUE),\n"
    "    mean(stand_init_fvssi$CSI, na.rm=TRUE)))\n"
)
s = s.replace(anchor, inject, 1)

# Replace the run block: switch stand_init between three configs
import re
new_block = (
    'rows <- list()\n'
    '# v27b: three configs switching stand_init between three productivity\n'
    '# metrics. All else at 12.3.8 production posture (MORTCAL on, CutPoint=0).\n'
    '\n'
    'stand_init <- stand_init_csi\n'
    'cat("[v27b] csi_baseline    (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, CSI = ClimateSI_ft)\\n")\n'
    'rows[["a"]] <- summ("csi_baseline",    run_cfg(TRUE, "Y", 0))\n'
    '\n'
    'stand_init <- stand_init_sicond\n'
    'cat("[v27b] sicond_replace  (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, CSI = SICOND)\\n")\n'
    'rows[["b"]] <- summ("sicond_replace",  run_cfg(TRUE, "Y", 0))\n'
    '\n'
    'stand_init <- stand_init_fvssi\n'
    'cat("[v27b] fvssi_replace   (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, CSI = FVS_SITE_INDEX)\\n")\n'
    'rows[["c"]] <- summ("fvssi_replace",   run_cfg(TRUE, "Y", 0))\n'
    '\n'
    'res <- dplyr::bind_rows(rows)\n'
    'write.csv(res, file.path(OUT_DIR, "acadgy_siteindex_v27b_results.csv"), row.names=FALSE)\n'
    'cat("\\n=== FIA: do alternative site indices (SICOND, FVS_SITE_INDEX) close BA vs ClimateSI? ===\\n")\n'
    'print(format(res, digits=4))'
)
pat = re.compile(r'rows <- list\(\)\s*\ncat\("\[v27b\] baseline.*?print\(format\(res, digits=4\)\)',
                 re.DOTALL)
m = pat.search(s)
if m:
    s = pat.sub(new_block, s, count=1)
else:
    s = s + "\n" + new_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "stand_init_csi_ok", "stand_init_csi" in s,
      "sicond_ok", "sicond_replace" in s,
      "fvssi_ok", "fvssi_replace" in s)
