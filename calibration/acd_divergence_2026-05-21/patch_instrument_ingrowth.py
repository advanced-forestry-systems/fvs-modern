#!/usr/bin/env python3
# Build AcadianGY_12.3.7_dbg.r = 12.3.7 + a cat() inside the ingrowth block that
# prints, per stand call, the NA counts on every covariate that feeds
# Ingrowth.FUN plus the IPH range and recruited row count. This pins down which
# input goes NA on FIA stands and silently zeroes recruitment.
src = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.7.r"
dst = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.7_dbg.r"
s = open(src).read()

# Anchor: the closing line of the mapply(Ingrowth.FUN, ...) assignment to
# Sum.temp$IPH. We need a unique substring; the cyclen)) close + newline is the
# end of that statement in 12.3.6/12.3.7.
import re
m = re.search(r"Sum\.temp\$IPH\s*=\s*as\.numeric\(mapply\(Ingrowth\.FUN[\s\S]*?cyclen\)\)\s*\n", s)
assert m, "Ingrowth.FUN mapply anchor not found"
anchor = m.group(0)

probe = (
"  # ---- #127 part 2 INSTRUMENTATION -----------------------------------------\n"
"  .ing_cols <- intersect(c(\"BAPH\",\"tph\",\"qmd\",\"pHW.ba\",\"pHW.tph\",\"CSI\",\"ELEV\"), names(Sum.temp))\n"
"  .ing_na <- sapply(.ing_cols, function(cn) sum(is.na(Sum.temp[[cn]])))\n"
"  .iph_na <- sum(is.na(Sum.temp$IPH))\n"
"  .iph_rng <- if (any(!is.na(Sum.temp$IPH))) paste(round(range(Sum.temp$IPH, na.rm=TRUE),4), collapse=\"..\") else \"all-NA\"\n"
"  cat(sprintf(\"[INGRWDBG] nrow=%d  IPH NAs=%d  IPH range=%s  covNA=%s\\n\",\n"
"    nrow(Sum.temp), .iph_na, .iph_rng,\n"
"    paste(sprintf(\"%s:%d\", .ing_cols, .ing_na), collapse=\",\")))\n"
)
s = s.replace(anchor, anchor + probe, 1)

# Mark the tag so unit tests can see they're on the dbg build
s = s.replace('AcadianVersionTag = "AcadianV12.3.7"',
              'AcadianVersionTag = "AcadianV12.3.7-dbg"', 1)

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "tag_ok", 'AcadianV12.3.7-dbg' in s,
      "probe_ok", "[INGRWDBG]" in s)
