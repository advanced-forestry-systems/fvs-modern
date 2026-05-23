#!/usr/bin/env python3
# Build AcadianGY_12.3.7.r = 12.3.6 + the #127 ingrowth carry-through fix.
# ING.TreeList recruits trees without dDBH.mult/dHt.mult/mort.mult/max.dbh/
# max.height, so on the next cycle dDBH*NA -> NA -> coalesced to 0 (recruits
# never grow). Set neutral multipliers + species size caps on `ingrow` before it
# is bound in, so recruited trees grow into the stand.
src = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.6.r"
dst = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.7.r"
s = open(src).read()

anchor = "  tree = dplyr::bind_rows(tree, ingrow)"
assert anchor in s, "bind_rows(tree, ingrow) anchor not found"
fix = (
"  # ---- #127 ingrowth carry-through fix --------------------------------------\n"
"  # ING.TreeList recruits lack the calibration-multiplier and size-cap columns,\n"
"  # so without this they never grow (dDBH * NA -> coalesced to 0). Give recruits\n"
"  # neutral multipliers and species size caps inherited from the survivors.\n"
"  if (!is.null(ingrow) && nrow(ingrow) > 0) {\n"
"    .cap <- unique(tree[, c(\"SP\",\"max.dbh\",\"max.height\")])\n"
"    .cap <- .cap[!duplicated(.cap$SP), ]\n"
"    .mi  <- match(ingrow$SP, .cap$SP)\n"
"    ingrow$dDBH.mult  <- 1\n"
"    ingrow$dHt.mult   <- 1\n"
"    ingrow$mort.mult  <- 1\n"
"    ingrow$max.dbh    <- ifelse(is.na(.cap$max.dbh[.mi]),    200, .cap$max.dbh[.mi])\n"
"    ingrow$max.height <- ifelse(is.na(.cap$max.height[.mi]),  60, .cap$max.height[.mi])\n"
"  }\n"
+ anchor)
s = s.replace(anchor, fix, 1)

s = s.replace('AcadianVersionTag = "AcadianV12.3.6"', 'AcadianVersionTag = "AcadianV12.3.7"', 1)
# the 12.3.6 file may carry the 12.3.5 tag (uploaded) or 12.3.6 (tag-fixed); handle both
if 'AcadianVersionTag = "AcadianV12.3.7"' not in s:
    s = s.replace('AcadianVersionTag = "AcadianV12.3.5"', 'AcadianVersionTag = "AcadianV12.3.7"', 1)

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "tag_ok", 'AcadianVersionTag = "AcadianV12.3.7"' in s,
      "fix_ok", "#127 ingrowth carry-through fix" in s)
