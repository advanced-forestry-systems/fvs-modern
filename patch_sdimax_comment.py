#!/usr/bin/env python3
# Correct the SDIMAX disable comment: the root cause was NOT field order (it is
# species, value, which is correct) but column misalignment -- the keyword was
# written as "SDIMAX"+10 spaces (16-char prefix), pushing the fields out of their
# fixed 10-col slots. _format_sdimax_keywords is now fixed (10-col keyword).
p = "/users/PUOM0008/crsfaaron/fvs-modern/config/config_loader.py"
s = open(p).read()
old = (
'        # SDIMAX keyword: DISABLED 2026-06-16 (WO-1, A. Weiskittel sign-off).\n'
'        # The per-species emission used the wrong FVS keyword field order\n'
'        # (species index in field 1, which FVS reads as the max-SDI VALUE),\n'
'        # so FVS set MAX SDI ~= 1 for all species and over-thinned TPH by\n'
'        # 25-35 percent across variants while only restating native defaults.\n'
'        # Re-enable ONLY with corrected field order (field 1 = value, field 2\n'
'        # = species) AND the revised localized max-SDI values, verified to bind\n'
'        # (per-species stand max is BA-weighted with retention flags).\n'
)
new = (
'        # SDIMAX keyword: FORMAT FIXED 2026-06-29 (was DISABLED 2026-06-16, WO-1).\n'
'        # Root cause (initre.f90 option 89): the field ORDER is correct (species,\n'
'        # value); the real bug was the keyword written as "SDIMAX"+10 spaces, a\n'
'        # 16-char prefix that pushed species/value out of their fixed 10-col fields\n'
'        # so FVS misread them (garbage MAX SDI -> over-thinning). _format_sdimax_\n'
'        # keywords now left-justifies the keyword to 10 cols (matches the tested\n'
'        # %-10s%10d%10.1f in sdimax_binding_test.py). Localized max-SDI surfaces\n'
'        # exist (brms_SDImax_site_specific.csv, alphaearth maxsdi_*_lookup.csv);\n'
'        # NA dropouts handled by make_sdifix. Note sdimax is often non-binding\n'
'        # (growth-engine-dominated, per the 2026-06-10 audit). Re-enable per variant\n'
'        # via config "_emit_sdimax": true once the NE binding test confirms binding.\n'
)
if old in s:
    s = s.replace(old, new); open(p, "w").write(s); print("comment patched")
else:
    print("comment block not found (already patched or drifted)")
