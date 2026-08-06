#!/usr/bin/env python3
"""Diagnose why the fvs_regional calibration block emits no BAIMULT/HTGMULT/BAMAX."""
import json, sys, os
import numpy as np

PROJECT_ROOT = "/users/PUOM0008/crsfaaron/fvs-modern"
CONFIG_DIR = os.environ.get("CONFIG_DIR", f"{PROJECT_ROOT}/config")
sys.path.insert(0, PROJECT_ROOT)
from config.config_loader import FvsConfigLoader

out = {}
loader = FvsConfigLoader("ne", version="calibrated", config_dir=CONFIG_DIR)
cfg = loader.config
cats = cfg.get("categories", {})
pre = cfg.get("calibration_multipliers") or {}

out["config_path"] = str(getattr(loader, "config_path", "unknown"))
out["top_level_keys"] = sorted(cfg.keys())
out["category_keys"] = sorted(cats.keys())
out["calibration_multipliers_keys"] = sorted(pre.keys())
out["_emit_sdimax"] = cfg.get("_emit_sdimax", True)


def describe(name, arr):
    if arr is None:
        return {"present": False}
    a = np.asarray(arr, dtype=float)
    finite = a[np.isfinite(a)]
    return {
        "present": True,
        "n": int(a.size),
        "n_nan": int(np.sum(~np.isfinite(a))),
        "min": float(np.min(finite)) if finite.size else None,
        "max": float(np.max(finite)) if finite.size else None,
        "n_gt_0p01_from_1": int(np.sum(np.abs(a - 1.0) > 0.01)) if finite.size else 0,
        "first10": [None if not np.isfinite(x) else round(float(x), 6) for x in a[:10]],
    }


# raw precomputed
for k in ("mort_multiplier", "dds_multiplier", "htg_multiplier"):
    out[f"raw_{k}"] = describe(k, pre.get(k))

mort = loader._array_or_none(pre.get("mort_multiplier"))
if mort is None:
    mort = loader._compute_mortality_multipliers(cats)
    out["mort_source"] = "legacy_compute"
else:
    out["mort_source"] = "precomputed"
out["mort_multiplier"] = describe("mort", mort)

dds = loader._array_or_none(pre.get("dds_multiplier"))
if dds is not None:
    growth = np.sqrt(np.clip(dds, 0.01, 100.0))
    out["dg_source"] = "precomputed_dds_sqrt"
else:
    growth = loader._compute_growth_multipliers(cats)
    out["dg_source"] = "legacy_compute"
out["growth_multiplier"] = describe("growth", growth)

hg = loader._array_or_none(pre.get("htg_multiplier"))
if hg is None:
    hg = loader._compute_height_multipliers(cats)
    out["hg_source"] = "legacy_compute"
else:
    out["hg_source"] = "precomputed"
out["htg_multiplier"] = describe("htg", hg)

bam = loader._find_bamax_param(cats)
out["bamax_param"] = None if bam is None else (str(type(bam)), str(bam)[:200])

sdi = loader._find_sdi_param(cats)
out["sdi_param_n"] = None if sdi is None else len(sdi)

# What actually gets emitted
blocks = {}
blocks["sdimax"] = loader._format_sdimax_keywords(sdi, False) if sdi is not None else ""
blocks["bamax"] = loader._format_bamax_keywords(bam, False) if bam is not None else ""
blocks["mortmult"] = loader._format_mortmult_keywords(mort, False) if mort is not None else ""
blocks["baimult"] = loader._format_baimult_keywords(growth, False) if growth is not None else ""
blocks["htgmult"] = loader._format_htgmult_keywords(hg, False) if hg is not None else ""
out["emitted_record_counts"] = {
    k: (0 if not v.strip() else len([ln for ln in v.splitlines() if ln.strip() and not ln.startswith("!!")]))
    for k, v in blocks.items()
}
out["sample_records"] = {k: (v.splitlines()[:3] if v.strip() else []) for k, v in blocks.items()}

# Column dissection of one record of each family
def cols(rec):
    r = rec.ljust(80)
    return {
        "keyword_1_8": repr(r[0:8]),
        "f1_11_20": repr(r[10:20]),
        "f2_21_30": repr(r[20:30]),
        "f3_31_40": repr(r[30:40]),
        "f4_41_50": repr(r[40:50]),
        "f5_51_60": repr(r[50:60]),
    }

out["column_dissection"] = {}
for k, v in blocks.items():
    recs = [ln for ln in v.splitlines() if ln.strip() and not ln.startswith("!!")]
    if recs:
        out["column_dissection"][k] = cols(recs[0])

print(json.dumps(out, indent=2, default=str))
