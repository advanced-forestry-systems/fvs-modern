# Max-SDI surface: found, root-caused, fixed, and validated (2026-06-29)

Triggered by "search Cardinal home and scratch for a max SDI surface." The surfaces exist; the SDIMAX
keyword bug is diagnosed, fixed, and the fix validated to bind. This resolves the SDIMAX half of task #30.

## 1. Localized max-SDI surfaces found (they exist)
- fvs-modern/diagnostics_2026-06-16/brms_SDImax_site_specific.csv  (site-specific brms SDImax)
- Disturbance/brms.SDImax.1-27-24.csv  and  conus_hcs/data/raw/sdi_max/brms.SDImax.1-27-24.csv (brms quantile surface)
- alphaearth CONUS lookups (scratch/Disturbance_run/alphaearth/conus/): maxsdi_ecoregion_fortype_lookup.csv,
  maxsdi_ecoregion_sppgrp_lookup.csv, maxsdi_fortype_lookup.csv, maxsdi_typegroup_lookup.csv, maxsdi_sppgrp_lookup.csv
- fia_cem_projections/config/: sdimax_brms_plot.csv, sdimax_brms_county_fortyp.csv, sdimax_by_l3_ecoregion.csv, etc.
- TREEMAP/ME.SDImax.tif, SiteIndex/TREEMAP_SDImax.tif (raster surfaces)
- The per-variant calibrated SDImax already lives in config/calibrated/*.json (the categories SDImax keys).

## 2. The real root cause (NOT field order)
The 2026-06-16 disable comment blamed "wrong field order (species in field 1 read as value)." That diagnosis
was inaccurate. Per sdimax_binding_test.py (citing initre.f90 option 89), the field ORDER is correct
(species, then value). The actual bug: config_loader emitted the keyword as "SDIMAX" + 10 literal spaces (a
16-char prefix), which pushed the species and value out of their fixed 10-column fields, so FVS parsed
garbage (MAX SDI ~ garbage -> over-thinning). FVS keyword fields are fixed 10-col: cols 1-10 keyword, 11-20
species, 21-30 value.

## 3. The fix (applied)
config_loader._format_sdimax_keywords line:
  was: f"SDIMAX          {i + 1:10d}{val:10.1f}"          # "SDIMAX" + 10 spaces = 16-char prefix (broken)
  now: f"{'SDIMAX':<10}{i + 1:10d}{val:10.1f}"            # keyword left-justified to 10 cols (correct)
This matches the tested format in sdimax_binding_test.py: "%-10s%10d%10.1f" % ("SDIMAX", species, value).
The disable comment was corrected to the true root cause.

## 4. Validation: the corrected format BINDS (NE binding test, 40 FIA plots, 50 yr)
| arm | mean final TPH | median |
|---|---|---|
| default (no SDIMAX) | 862.6 | 838.3 |
| SDImax x0.6 | 492.3 | 400.2 |
| SDImax x1.0 (native calibrated) | 755.7 | 733.7 |
| SDImax x1.6 | 1221.6 | 1215.1 |
Lower max SDI -> lower final density, monotonically; native (x1.0) constrains density below the unconstrained
default (756 vs 863). Verdict: BINDING. The corrected SDIMAX works correctly on the eastern/CONUS variant
(the old over-thinning was the misaligned format, not the values).

Caveat (from the 2026-06-10 audit): in the WEST (WS/CA carbon to 2100) the ceiling was non-binding -- those
stands are growth-engine-dominated and never reach the limit. So SDIMAX matters where stands approach the
self-thinning limit (denser eastern stands here) and is secondary where growth dynamics dominate.

## 5. Re-enable (ready, one flag)
Re-enable is gated by the per-config flag `_emit_sdimax` (default False). With the format now fixed and binding
confirmed, set `"_emit_sdimax": true` in the CONUS configs to emit corrected SDIMAX keywords. Left as Aaron's
explicit go since it changes behavior across ~26 configs and was a deliberate WO-1 sign-off disable; the
blocker (broken format) is now removed and the fix is validated. NA-species dropouts are handled by
make_sdifix (calibrated_sdifix configs).

## 6. Relation to options 2/3 (the Garcia BA carrying-capacity, the other half of #30)
SDIMAX is the native per-species self-thinning ceiling. Options 2/3 use the fvs-conus Garcia stand layer
(self-thinning power law + monomolecular BA carrying capacity toward Gmax = f(SDImax)) via the projector, which
is the richer density regulator. The BA carrying-capacity hook for the native engine remains the second half of
#30; the localized SDImax surfaces above feed Gmax there too.
