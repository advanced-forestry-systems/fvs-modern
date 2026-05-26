#!/usr/bin/env python3
# v22: same harness shape as v21 but (a) sources the instrumented 12.3.7-dbg,
# (b) drops to 10 plots, (c) prints per-cycle TPA/added per stand from p1y, and
# (d) keeps n_years=3 so we can see across-cycle behavior without flooding logs.
# Goal: catch the exact step where v21's CutPoint=0 run loses the ingrowth that
# the standalone probe proves IS produced on cycle 1.
src = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_ingrowthfix_v21.R"
dst = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/cardinal_acadgy_ingrowthfix_v22.R"
s = open(src).read()

s = s.replace("/users/PUOM0008/crsfaaron/AcadianGY_12.3.7.r",
              "/users/PUOM0008/crsfaaron/AcadianGY_12.3.7_dbg.r", 1)
s = s.replace("N_PLOTS <- 100L", "N_PLOTS <- 10L", 1)
s = s.replace("n_years <- 10L", "n_years <- 3L", 1)
s = s.replace("acadgy_ingrowthfix_v21_results.csv",
              "acadgy_ingrowthfix_v22_results.csv", 1)
s = s.replace("[v21] baseline", "[v22] baseline", 1)
s = s.replace("[v21] ingrowthfix", "[v22] ingrowthfix", 1)

# Add per-cycle per-stand TPA prints inside p1y so we can see exactly when/where
# recruits appear and whether they survive the next cycle.
old = ("  for (sid in unique(trees$STAND)) {\n"
       "    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]\n"
       "    if (nrow(sub) == 0) next\n"
       "    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n"
       "    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }\n"
       "    pr <- tryCatch(AcadianGYOneStand(sub, stand = st, ops = ops), error = function(e) NULL)\n"
       "    if (!is.null(pr)) pc[[sid]] <- pr\n"
       "  }\n")
new = ("  for (sid in unique(trees$STAND)) {\n"
       "    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]\n"
       "    if (nrow(sub) == 0) next\n"
       "    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint\n"
       "    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }\n"
       "    cat(sprintf(\"[V22 call] sid=%s nrow_in=%d cp=%s mort=%s\\n\", sid, nrow(sub), as.character(ops$CutPoint %||% \"NA\"), as.character(isTRUE(ops$MORTCAL))))\n"
       "    pr <- tryCatch(AcadianGYOneStand(sub, stand = st, ops = ops), error = function(e) { cat(\"[V22 err] \", conditionMessage(e), \"\\n\"); NULL })\n"
       "    if (!is.null(pr)) { cat(sprintf(\"[V22 ret] sid=%s nrow_out=%d delta=%d\\n\", sid, nrow(pr), nrow(pr)-nrow(sub))); pc[[sid]] <- pr }\n"
       "  }\n")
assert old in s, "p1y loop body not found"
s = s.replace(old, new, 1)

# %||% helper (R lacks it natively unless rlang); add a tiny one at the top
s = s.replace("suppressMessages({ library(dplyr); library(plyr); library(purrr) })",
              "suppressMessages({ library(dplyr); library(plyr); library(purrr) })\n`%||%` <- function(a,b) if (is.null(a)) b else a", 1)

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "dbg_ok", "12.3.7_dbg.r" in s,
      "v22_ok", "[V22 call]" in s,
      "nplots_ok", "N_PLOTS <- 10L" in s,
      "nyears_ok", "n_years <- 3L" in s)
