#!/usr/bin/env python3
# Build AcadianGY_12.3.9.r = 12.3.8 + optional ops$CSI_SCALE knob. Defaults to
# 1.0 (no change), so 12.3.9 with no CSI_SCALE set is byte-identical to 12.3.8.
# The knob multiplies the resolved CSI value once at the top of
# AcadianGYOneStand and then everything downstream (Kuehne et al. dDBH, Russell
# /Weiskittel dHt, Li et al. Ingrowth.FUN, Ingrowth.Comp) sees the scaled CSI.
# v25 found that CSI x 0.7 closes BA bias from +11.1% to about +10.4% and
# improves R^2 by ~0.01 on the FIA harness.
src = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.8.r"
dst = "/users/PUOM0008/crsfaaron/AcadianGY_12.3.9.r"
s = open(src).read()

anchor = "CSI      = if (is.null(stand) || is.null(stand$CSI) || is.na(stand$CSI))   12 else stand$CSI"
assert anchor in s, "CSI resolution line not found"
inject = (
    anchor + "\n"
    "  # ---- 12.3.9: ops$CSI_SCALE knob ----------------------------------------\n"
    "  # Optional scalar multiplier on Climate Site Index. Defaults to 1.0 (no\n"
    "  # change). v25 sensitivity showed CSI elasticity ~0.27 pp BA per 0.1 scale;\n"
    "  # CSI*0.7 closes BA bias by ~0.6 pp and improves R^2 ~0.01 on ME FIA.\n"
    "  if (!is.null(ops$CSI_SCALE) && is.finite(ops$CSI_SCALE) && ops$CSI_SCALE > 0) {\n"
    "    CSI <- CSI * ops$CSI_SCALE\n"
    "  }"
)
s = s.replace(anchor, inject, 1)

s = s.replace('AcadianVersionTag = "AcadianV12.3.8"',
              'AcadianVersionTag = "AcadianV12.3.9"', 1)

open(dst, "w").write(s)
print("wrote", dst, "len", len(s),
      "tag_ok", 'AcadianV12.3.9' in s,
      "knob_ok", "12.3.9: ops$CSI_SCALE knob" in s,
      "scale_ok", "CSI <- CSI * ops$CSI_SCALE" in s)
