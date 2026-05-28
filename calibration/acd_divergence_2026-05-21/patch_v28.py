#!/usr/bin/env python3
# Build cardinal_acadgy_csiscale_v28.R = v24 sourcing AcadianGY_12.3.9.r with
# three configs that exercise ops$CSI_SCALE. Default (no scale set) must
# byte-match v24 to confirm backward compat. CSI_SCALE=0.7 is the recommended
# production setting from v25.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_ingrowthfix_v24.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_csiscale_v28.R"
s = open(src).read()

s = s.replace("AcadianGY_12.3.8.r", "AcadianGY_12.3.9.r")
s = s.replace("[v24]", "[v28]")
s = s.replace("acadgy_ingrowthfix_v24_results.csv", "acadgy_csiscale_v28_results.csv")

# Extend run_cfg to pass CSI_SCALE through ops, and use it in the new configs
old_runcfg = ("p1y <- function(trees, mortcal, ingrowth=\"Y\", cutpoint=NULL) {\n"
              "  pc <- list()\n"
              "  for (sid in unique(trees$STAND)) {\n"
              "    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]\n"
              "    if (nrow(sub) == 0) next\n"
              "    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n"
              "    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }")
new_runcfg = ("p1y <- function(trees, mortcal, ingrowth=\"Y\", cutpoint=NULL, csi_scale=NULL) {\n"
              "  pc <- list()\n"
              "  for (sid in unique(trees$STAND)) {\n"
              "    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]\n"
              "    if (nrow(sub) == 0) next\n"
              "    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n"
              "    if (!is.null(csi_scale)) ops$CSI_SCALE <- csi_scale\n"
              "    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }")
assert old_runcfg in s, "p1y signature anchor not found"
s = s.replace(old_runcfg, new_runcfg, 1)

# Threading through run_cfg signature too
s = s.replace(
    'run_cfg <- function(mortcal, ingrowth="Y", cutpoint=NULL) {',
    'run_cfg <- function(mortcal, ingrowth="Y", cutpoint=NULL, csi_scale=NULL) {', 1)
s = s.replace(
    'cur <- p1y(cur, mortcal, ingrowth, cutpoint); if (is.null(cur)) break',
    'cur <- p1y(cur, mortcal, ingrowth, cutpoint, csi_scale); if (is.null(cur)) break', 1)

# Replace the run block with three configs
import re
pat = re.compile(r'rows <- list\(\)\s*\ncat\("\[v28\] baseline.*?print\(format\(res, digits=4\)\)',
                 re.DOTALL)
new_block = (
    'rows <- list()\n'
    'cat("[v28] csi_scale_1.0  (12.3.9, MORTCAL on, CutPoint 0, CSI_SCALE not set; must match v24)\\n")\n'
    'rows[["a"]] <- summ("csi_scale_1.0",  run_cfg(TRUE, "Y", 0, NULL))\n'
    'cat("[v28] csi_scale_0.7  (12.3.9, MORTCAL on, CutPoint 0, CSI_SCALE = 0.7)\\n")\n'
    'rows[["b"]] <- summ("csi_scale_0.7",  run_cfg(TRUE, "Y", 0, 0.7))\n'
    'cat("[v28] csi_scale_0.5  (12.3.9, MORTCAL on, CutPoint 0, CSI_SCALE = 0.5)\\n")\n'
    'rows[["c"]] <- summ("csi_scale_0.5",  run_cfg(TRUE, "Y", 0, 0.5))\n'
    'res <- dplyr::bind_rows(rows)\n'
    'write.csv(res, file.path(OUT_DIR, "acadgy_csiscale_v28_results.csv"), row.names=FALSE)\n'
    'cat("\\n=== FIA: validate 12.3.9 ops$CSI_SCALE knob against v24 baseline and v25 sweep ===\\n")\n'
    'print(format(res, digits=4))'
)
if pat.search(s):
    s = pat.sub(new_block, s, count=1)
else:
    s = s + "\n" + new_block + "\n"

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "model_ok", "AcadianGY_12.3.9.r" in s,
      "csi_scale_threaded", "csi_scale" in s,
      "configs_ok", "csi_scale_1.0" in s and "csi_scale_0.7" in s)
