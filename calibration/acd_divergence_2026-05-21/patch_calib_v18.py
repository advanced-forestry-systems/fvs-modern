#!/usr/bin/env python3
# Build cardinal_acadgy_calib_v18.R = v17 + source 12.3.8 + a 4th config that
# layers CutPoint=0 expected-value ingrowth on top of calibrated_on. v17 was
# the source of the "diameter calibration worsens stand BA" finding, but that
# was measured on 12.3.6 where recruits were silently dropped. v18 tests
# whether layering ingrowth on top of calibrated diameter growth closes BA
# below the +13.3% residual.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_calib_v17.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_calib_v18.R"
s = open(src).read()

# Model swap
s = s.replace('source(file.path(PROJECT_ROOT, "seven_islands/AcadianGY_12.3.6.r"))',
              'source("/users/PUOM0008/crsfaaron/AcadianGY_12.3.8.r")', 1)
s = s.replace("[v17]", "[v18]")
s = s.replace("acadgy_calib_v17_results.csv", "acadgy_calib_v18_results.csv")

# Extend p1y and run_cfg to accept cutpoint
old_p1y_sig = "p1y <- function(trees, mortcal, calib=FALSE) {"
new_p1y_sig = "p1y <- function(trees, mortcal, calib=FALSE, cutpoint=NULL) {"
assert old_p1y_sig in s
s = s.replace(old_p1y_sig, new_p1y_sig, 1)

old_ops_assign = "    ops <- ops0\n    if (mortcal)"
new_ops_assign = "    ops <- ops0\n    if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n    if (mortcal)"
assert old_ops_assign in s
s = s.replace(old_ops_assign, new_ops_assign, 1)

old_p1y_call = "    cur <- p1y(cur, mortcal, calib)"
new_p1y_call = "    cur <- p1y(cur, mortcal, calib, cutpoint)"
assert old_p1y_call in s
s = s.replace(old_p1y_call, new_p1y_call, 1)

old_runcfg_sig = "run_cfg <- function(mortcal, calib=FALSE) {\n  cur <- base_init;"
new_runcfg_sig = "run_cfg <- function(mortcal, calib=FALSE, cutpoint=NULL) {\n  cur <- base_init;"
assert old_runcfg_sig in s
s = s.replace(old_runcfg_sig, new_runcfg_sig, 1)

# Replace the rows assembly with five configs: three sanity checks against v17
# + one new config with ingrowth + one ingrowth-only baseline for comparison
import re
old_run_block_pat = re.compile(
    r'cat\("\[v18\] canonical_off.*?print\(format\(res, digits=4\)\)',
    re.DOTALL
)
new_run_block = (
    'cat("[v18] canonical_off    (12.3.8, no calib, MORTCAL off, default CP)\\n");   rows[["a"]] <- summ("canonical_off",    run_cfg(FALSE, FALSE, NULL))\n'
    'cat("[v18] calibrated_off   (12.3.8, dia calib,  MORTCAL off, default CP)\\n"); rows[["b"]] <- summ("calibrated_off",   run_cfg(FALSE, TRUE,  NULL))\n'
    'cat("[v18] calibrated_on    (12.3.8, dia calib + MORTCAL,    default CP)\\n");  rows[["c"]] <- summ("calibrated_on",    run_cfg(TRUE,  TRUE,  NULL))\n'
    'cat("[v18] calibrated_on_cp0(12.3.8, dia calib + MORTCAL + EV ingrowth)\\n");   rows[["d"]] <- summ("calibrated_on_cp0",run_cfg(TRUE,  TRUE,  0))\n'
    'cat("[v18] ingrowth_only    (12.3.8, no calib,  MORTCAL off, EV ingrowth)\\n"); rows[["e"]] <- summ("ingrowth_only",    run_cfg(FALSE, FALSE, 0))\n'
    'res <- dplyr::bind_rows(rows); write.csv(res, file.path(OUT_DIR, "acadgy_calib_v18_results.csv"), row.names=FALSE)\n'
    'cat("\\n=== FIA: does ingrowth on top of diameter calibration close stand BA? ===\\n")\n'
    'print(format(res, digits=4))'
)
matched = old_run_block_pat.search(s)
if matched:
    s = old_run_block_pat.sub(new_run_block, s, count=1)
else:
    # Fallback: append at end
    s = s + "\n" + new_run_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "model_ok", "AcadianGY_12.3.8.r" in s,
      "cp0_ok", "calibrated_on_cp0" in s,
      "ingrowth_only_ok", "ingrowth_only" in s,
      "p1y_sig_ok", "p1y <- function(trees, mortcal, calib=FALSE, cutpoint=NULL)" in s)
