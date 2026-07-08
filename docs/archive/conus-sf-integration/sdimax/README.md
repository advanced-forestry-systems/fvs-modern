# Localized maximum SDI

Maximum SDI is not observable; it is estimated. The FVS species-weighted maximum is biased
about 28 percent high and has near-zero plot-level skill (bias-corrected R2 0.02). A localized,
FIA-derived maximum predicts observed FIA self-thinning about 85 percent better (deviance
explained 0.11 vs 0.06), in every region. See `../CONUS_MAXSDI_TECHNICAL_REPORT.md`.

`localized_sdimax.py` returns a per-stand maximum SDI (trees/ha) from, in order of fidelity:
the TreeMap 2022 30 m SDImax raster (Zenodo 10.5281/zenodo.19509367), the brms FIA plot table,
or a forest-type + geography fallback. It also emits the per-stand FVS SDIMAX keyword block.
Model-agnostic: any growth-and-yield engine consumes the per-stand value as its density limit.
