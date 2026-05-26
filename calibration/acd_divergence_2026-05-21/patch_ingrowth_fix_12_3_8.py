#!/usr/bin/env python3
# Build AcadianGY_12.3.8.r = 12.3.7 + STAND/PLOT carry-through for ingrowth.
# The 12.3.7 block sets dDBH.mult/dHt.mult/mort.mult/max.dbh/max.height so
# recruits grow; this adds STAND and PLOT inheritance from the survivors so
# multi-cycle harnesses don't fragment them off into a phantom stand="1" group
# that the next cycle's per-stand dispatcher cannot find.
src = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.7.r"
dst = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.8.r"
s = open(src).read()

anchor = (
"    ingrow$dDBH.mult  <- 1\n"
"    ingrow$dHt.mult   <- 1\n"
"    ingrow$mort.mult  <- 1\n"
"    ingrow$max.dbh    <- ifelse(is.na(.cap$max.dbh[.mi]),    200, .cap$max.dbh[.mi])\n"
"    ingrow$max.height <- ifelse(is.na(.cap$max.height[.mi]),  60, .cap$max.height[.mi])\n"
"  }\n"
)
assert anchor in s, "12.3.7 ingrow neutral-multiplier block not found"
fix = (
"    ingrow$dDBH.mult  <- 1\n"
"    ingrow$dHt.mult   <- 1\n"
"    ingrow$mort.mult  <- 1\n"
"    ingrow$max.dbh    <- ifelse(is.na(.cap$max.dbh[.mi]),    200, .cap$max.dbh[.mi])\n"
"    ingrow$max.height <- ifelse(is.na(.cap$max.height[.mi]),  60, .cap$max.height[.mi])\n"
"    # ---- 12.3.8: STAND/PLOT inheritance ----------------------------------\n"
"    # AcadianGYOneStand operates on one stand at a time, so survivors share a\n"
"    # single STAND value. ING.TreeList defaults recruits to STAND=1/PLOT=1,\n"
"    # which fragments them off into a phantom stand on the next cycle of any\n"
"    # multi-stand harness. Force inheritance so recruits stay with the parent.\n"
"    if (\"STAND\" %in% names(tree) && length(unique(tree$STAND)) == 1) {\n"
"      ingrow$STAND <- unique(tree$STAND)\n"
"    }\n"
"    if (\"PLOT\" %in% names(tree) && length(unique(tree$PLOT)) == 1) {\n"
"      ingrow$PLOT <- unique(tree$PLOT)\n"
"    }\n"
"  }\n"
)
s = s.replace(anchor, fix, 1)

s = s.replace('AcadianVersionTag = "AcadianV12.3.7"',
              'AcadianVersionTag = "AcadianV12.3.8"', 1)

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "tag_ok", 'AcadianV12.3.8' in s,
      "stand_inherit_ok", "12.3.8: STAND/PLOT inheritance" in s)
