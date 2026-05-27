#!/usr/bin/env python3
# Build cardinal_magplot_insource_v17.R = v16 + (a) source 12.3.8 instead of
# 12.3.6, (b) add a third config insource_on_cp0 that forces CutPoint = 0
# (expected-value ingrowth). The decisive test: with the carry-through bugs
# fixed and CutPoint = 0 actually producing recruits per cycle, does MORTCAL
# still over-correct Canadian MAGPlot stands, or do recruits compensate?
src = "/users/PUOM0008/crsfaaron/magplot/verify/cardinal_magplot_insource_v16.R"
dst = "/users/PUOM0008/crsfaaron/magplot/verify/cardinal_magplot_insource_v17.R"
s = open(src).read()

# Model swap
s = s.replace('source(file.path(PROJECT_ROOT, "AcadianGY_12.3.6.r"))',
              'source(file.path(PROJECT_ROOT, "AcadianGY_12.3.8.r"))', 1)
s = s.replace("[v16-magplot]", "[v17-magplot]")
s = s.replace("insource_v16", "insource_v17")

# Extend run_cfg to accept a CutPoint argument
old_runcfg_sig = "run_cfg <- function(mortcal) {\n  ops <- list(verbose=FALSE, INGROWTH=\"Y\", MinDBH=3.0)\n  if (mortcal) { ops$MORTCAL <- TRUE; ops$MORTCAL_INTERVAL <- 5 }"
new_runcfg_sig = "run_cfg <- function(mortcal, cutpoint=NULL) {\n  ops <- list(verbose=FALSE, INGROWTH=\"Y\", MinDBH=3.0)\n  if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n  if (mortcal) { ops$MORTCAL <- TRUE; ops$MORTCAL_INTERVAL <- 5 }"
assert old_runcfg_sig in s, "run_cfg signature not found"
s = s.replace(old_runcfg_sig, new_runcfg_sig, 1)

# Extend the run block to include the third config
import re
# Match the rows assembly + write block. The v16 wrote two configs.
m = re.search(r"rows <- list\(\)\s*\n.*?cat\(\"\\n=== MAGPlot.*?\n", s, re.DOTALL)
if m:
    old_block = m.group(0)
else:
    # fallback: find the cat lines
    old_block = None

# Find existing run lines and inject a third
new_lines = (
"rows <- list()\n"
"cat(\"[v17-magplot] canonical_off  (12.3.8, MORTCAL off, default CutPoint)\\n\")\n"
"rows[[\"a\"]] <- summ(\"canonical_off\",   run_cfg(FALSE, NULL))\n"
"cat(\"[v17-magplot] insource_on    (12.3.8, MORTCAL on, default CutPoint)\\n\")\n"
"rows[[\"b\"]] <- summ(\"insource_on\",     run_cfg(TRUE,  NULL))\n"
"cat(\"[v17-magplot] insource_on_cp0(12.3.8, MORTCAL on, CutPoint=0 EV ingrowth)\\n\")\n"
"rows[[\"c\"]] <- summ(\"insource_on_cp0\", run_cfg(TRUE,  0))\n"
"cat(\"[v17-magplot] ingrowth_only  (12.3.8, MORTCAL off, CutPoint=0 EV ingrowth)\\n\")\n"
"rows[[\"d\"]] <- summ(\"ingrowth_only\",   run_cfg(FALSE, 0))\n"
"res <- dplyr::bind_rows(rows)\n"
"write.csv(res, file.path(OUT_DIR, \"magplot_insource_v17_results.csv\"), row.names=FALSE)\n"
"cat(\"\\n=== MAGPlot NB: does MORTCAL still over-correct once ingrowth flows? ===\\n\")\n"
"print(format(res, digits=4))\n"
)

# Replace any prior rows<-list() block with the new one (looking for the original two-row pattern)
pat = re.compile(r"rows <- list\(\).*?print\(format\(res[^\)]*\)\)\s*\n?", re.DOTALL)
s2 = pat.sub(new_lines, s, count=1)
if s2 == s:
    # Append the new block at the end if the pattern was structured differently
    s2 = s + "\n" + new_lines

open(dst, "w").write(s2)
print("wrote", dst, "len", len(s2),
      "model_ok", "AcadianGY_12.3.8.r" in s2,
      "cp0_ok", "insource_on_cp0" in s2,
      "ingrowth_only_ok", "ingrowth_only" in s2)
