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
