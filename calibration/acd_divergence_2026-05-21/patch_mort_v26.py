#!/usr/bin/env python3
# Build cardinal_acadgy_mortmult_v26.R = v24 (12.3.8 ingrowth-fix FIA baseline)
# with global mort.mult scaled by 1.0, 1.3, 1.6, 2.0, 2.5 in production posture.
# If higher mortality closes BA, the Glover/Hool form is under-killing trees and
# the lever is mortality rate. If flat or non-monotonic, the residual sits in
# dDBH base rate, not mortality. Same 100-plot sample as v24/v25, same n_years.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_ingrowthfix_v24.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_mortmult_v26.R"
s = open(src).read()

s = s.replace("[v24]", "[v26]")
s = s.replace("acadgy_ingrowthfix_v24_results.csv", "acadgy_mortmult_v26_results.csv")

# Replace the run block with a 5-config mort.mult sweep applied at tree-input
# stage. p1y currently doesn't know about mort.mult overrides, so set on
# base_init before each scale's run.
import re
m = re.search(r'rows <- list\(\)\s*\ncat\("\[v26\] baseline.*?print\(format\(res, digits=4\)\)',
              s, re.DOTALL)
old_block = m.group(0) if m else None
new_block = (
    'rows <- list()\n'
    '# v26: mort.mult sensitivity scan. Apply a global scale to base_init$mort.mult\n'
    '# before each run, then restore. Production posture (MORTCAL on, CutPoint=0\n'
    '# expected-value ingrowth) so we isolate the mortality-rate lever on the\n'
    '# residual BA.\n'
    'mm_base <- base_init$mort.mult\n'
    'for (scale in c(1.0, 1.3, 1.6, 2.0, 2.5)) {\n'
    '  base_init$mort.mult <- mm_base * scale\n'
    '  tag <- sprintf("mort_x%.1f", scale)\n'
    '  cat(sprintf("[v26] %s (12.3.8, CutPoint 0 EV ingrowth, MORTCAL on, mort.mult * %.1f)\\n", tag, scale))\n'
    '  rows[[tag]] <- summ(tag, run_cfg(TRUE, "Y", 0))\n'
    '}\n'
    'base_init$mort.mult <- mm_base\n'
    'res <- dplyr::bind_rows(rows)\n'
    'write.csv(res, file.path(OUT_DIR, "acadgy_mortmult_v26_results.csv"), row.names=FALSE)\n'
    'cat("\\n=== FIA mortality scan: does higher mort.mult close the structural BA residual? ===\\n")\n'
    'print(format(res, digits=4))'
)
if old_block:
    s = s.replace(old_block, new_block, 1)
else:
    s = s + "\n" + new_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "v26_ok", "[v26]" in s,
      "mort_sweep_ok", "mm_base * scale" in s)
