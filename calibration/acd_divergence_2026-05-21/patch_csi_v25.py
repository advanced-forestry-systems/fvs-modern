#!/usr/bin/env python3
# Build cardinal_acadgy_csisensitivity_v25.R = v24 (12.3.8 ingrowth-fix FIA
# baseline) with CSI scaled by 0.6, 0.8, 1.0, 1.2, 1.5 in production posture
# (calibrated_on_cp0 equivalent). Tests whether the structural +15 percent
# FIA BA residual responds to climate sensitivity. If BA bias drops
# substantially with lower CSI weighting, the Kuehne et al. dDBH climate
# term is over-weighted in the bridge. If it's flat, CSI is not the lever
# and we move to mortality functional form (or accept the ceiling).
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_ingrowthfix_v24.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_csisensitivity_v25.R"
s = open(src).read()

s = s.replace("[v24]", "[v25]")
s = s.replace("acadgy_ingrowthfix_v24_results.csv", "acadgy_csisensitivity_v25_results.csv")

# Replace the two-config run block with a 5-config CSI sweep
import re
# Find the rows<-list() to print block
m = re.search(r'rows <- list\(\)\s*\ncat\("\[v25\] baseline.*?print\(format\(res, digits=4\)\)',
              s, re.DOTALL)
old_block = m.group(0) if m else None
new_block = (
    'rows <- list()\n'
    '# v25: CSI sensitivity sweep at the FIA production posture (MORTCAL on,\n'
    '# CutPoint=0 expected-value ingrowth). Scale stand_init$CSI by 0.6, 0.8,\n'
    '# 1.0, 1.2, 1.5 and measure BA bias response.\n'
    'csi_base <- stand_init$CSI\n'
    'for (scale in c(0.6, 0.8, 1.0, 1.2, 1.5)) {\n'
    '  stand_init$CSI <- csi_base * scale\n'
    '  tag <- sprintf("csi_x%.1f", scale)\n'
    '  cat(sprintf("[v25] %s (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, CSI * %.1f)\\n", tag, scale))\n'
    '  rows[[tag]] <- summ(tag, run_cfg(TRUE, "Y", 0))\n'
    '}\n'
    'stand_init$CSI <- csi_base\n'
    'res <- dplyr::bind_rows(rows)\n'
    'write.csv(res, file.path(OUT_DIR, "acadgy_csisensitivity_v25_results.csv"), row.names=FALSE)\n'
    'cat("\\n=== FIA CSI sensitivity: does climate weighting close the structural BA residual? ===\\n")\n'
    'print(format(res, digits=4))'
)
if old_block:
    s = s.replace(old_block, new_block, 1)
else:
    # Fallback: append
    s = s + "\n" + new_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "v25_ok", "[v25]" in s,
      "csi_sweep_ok", "csi_base * scale" in s)
