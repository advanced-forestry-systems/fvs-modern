# ADK (Adirondack) variant: build status
2026-06-18. Scaffolded from Acadian via add_variant.sh, then built and verified.

- src-converted/adk/ created from acd/ (full variant source + calibration species-map template).
- Build manifest fix: add_variant.sh did not create src-converted/bin/FVSadk_sourceList.txt, so the
  first build skipped ADK. Created it from FVSacd_sourceList.txt with ../acd/ to ../adk/ path
  substitution. add_variant.sh should generate this manifest for every new variant; flagged.
- FVSadk.so builds clean (8.7 MB) and loads via ctypes (ADK_SO_LOADS).

Calibration next step: the per-species DG calibration harness (calib_ne.py) invokes the standalone
executable FVSiable, but the variant build produces the .so (the intended deployment path), not a
standalone executable. ADK calibration should run through the .so via fvs2py/rFVS, or compute the
per-species multipliers from NY Adirondack FIA observed growth against the acd-seeded default. NY FIA
provides the Adirondack remeasurement data; no external dependency.

## Calibration attempt (2026-06-18) and the real dependency

Tried the shortcut of calibrating the Acadian engine directly on NY (Adirondack) FIA (calib_ne.py
VAR=acd STATES=NY). It failed building the per-plot input table (sqlite syntax error on a degenerate
standinit): NY trees do not load cleanly through the Acadian variant species machinery, so the Acadian
engine cannot be used as-is on NY plots. A proper ADK calibration therefore needs (1) the ADK
species_map.csv edited for the NY/Adirondack species set, and (2) a way to RUN FVSadk: either build the
FVSadk standalone executable (the build script currently produces only the .so) or drive the .so in
process (the same fvs2py path as Route A). ADK is built and deployable as a .so; its calibration shares
the Route A / executable-build infrastructure dependency and is a focused variant-development task, not a
one-shot run.
