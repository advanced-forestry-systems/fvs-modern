#!/usr/bin/env python3
"""
Validate calibrated FVS configs against the adopted calibration.

This is the guardrail that prevents a recurrence of issue #54, where calibrated
configs silently failed to carry the components their metadata claimed, so only
SDIMAX reached FVS at runtime. It checks every config/calibrated/<variant>.json
for a consistent calibration_multipliers block and flags mismatches against the
authoritative equation_availability_full.csv.

Checks per variant:
  - maxsp and categories.species_definitions.FIAJSP are present.
  - calibration_multipliers carries the five component arrays
    (htdbh, mort, cr, dds, htg), each of length maxsp.
  - no NaN / inf values; all finite and positive (multipliers).
  - availability consistency: where a component is adopted (per the availability
    table) AND is SPCD-keyed (HD/MORT/CR), the array must carry real factors
    (not all 1.0). DG/HI may legitimately be all 1.0 where the fit was sparse,
    so those are warnings, not errors.

Exit code 0 if all hard checks pass, 1 otherwise. Intended for CI.

Usage:
  python config/validate_calibrated.py
  python config/validate_calibrated.py --config-dir config --quiet
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

COMPONENTS = {
    "htdbh_multiplier": "HD",
    "mort_multiplier": "MORT",
    "cr_multiplier": "CR",
    "dds_multiplier": "DG",
    "htg_multiplier": "HI",
}
SPCD_KEYED = {"htdbh_multiplier", "mort_multiplier", "cr_multiplier"}


def load_availability(calibration_dir: Path) -> dict[str, dict[str, bool]]:
    f = calibration_dir / "data" / "equation_availability_full.csv"
    table: dict[str, dict[str, bool]] = {}
    if not f.exists():
        return table
    with open(f) as fh:
        for row in csv.DictReader(fh):
            v = row["variant"].strip().lower()
            table[v] = {
                k: str(row.get(k, "")).strip().upper() == "TRUE"
                for k in ("HD", "MORT", "CR", "DG", "SDI", "HI")
            }
    return table


def validate_variant(path: Path, avail: dict[str, bool]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    d = json.load(open(path))
    maxsp = d.get("maxsp")
    if not isinstance(maxsp, int) or maxsp < 1:
        errors.append("missing or invalid maxsp")
        return errors, warnings

    cats = d.get("categories", {})
    sd = cats.get("species_definitions", {})
    if not isinstance(sd, dict) or "FIAJSP" not in sd:
        errors.append("categories.species_definitions.FIAJSP missing")

    cm = d.get("calibration_multipliers")
    if not isinstance(cm, dict):
        errors.append("calibration_multipliers block missing")
        return errors, warnings

    for key, comp in COMPONENTS.items():
        arr = cm.get(key)
        if not isinstance(arr, list):
            errors.append(f"{key} missing or not a list")
            continue
        if len(arr) != maxsp:
            errors.append(f"{key} length {len(arr)} != maxsp {maxsp}")
        bad = [x for x in arr if not isinstance(x, (int, float))
               or math.isnan(x) or math.isinf(x) or x <= 0]
        if bad:
            errors.append(f"{key} has {len(bad)} non-finite/non-positive values")
        off = sum(1 for x in arr if isinstance(x, (int, float)) and abs(x - 1.0) > 1e-6)
        adopted = avail.get(comp, False) if avail else None
        if adopted and off == 0:
            msg = f"{comp} adopted but {key} is all 1.0 (no factors propagated)"
            (errors if key in SPCD_KEYED else warnings).append(msg)
        if adopted is False and off > 0:
            warnings.append(f"{comp} not adopted but {key} carries {off} factors")
    return errors, warnings


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config-dir", default=str(Path(__file__).parent))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    config_dir = Path(args.config_dir)
    calibration_dir = config_dir.parent / "calibration"
    avail_table = load_availability(calibration_dir)
    if not avail_table and not args.quiet:
        print("WARN: equation_availability_full.csv not found; availability checks skipped")

    calibrated = sorted(p for p in (config_dir / "calibrated").glob("*.json")
                        if not p.stem.endswith("_draws"))
    if not calibrated:
        print("ERROR: no calibrated configs found")
        return 1

    total_err = 0
    for p in calibrated:
        v = p.stem
        errs, warns = validate_variant(p, avail_table.get(v, {}))
        total_err += len(errs)
        if errs:
            print(f"[FAIL] {v}")
            for e in errs:
                print(f"    ERROR: {e}")
        elif not args.quiet:
            print(f"[OK]   {v}" + (f"  ({len(warns)} warnings)" if warns else ""))
        if warns and not args.quiet:
            for w in warns:
                print(f"    warn: {w}")

    print(f"\n{len(calibrated)} configs checked, {total_err} hard errors")
    return 1 if total_err else 0


if __name__ == "__main__":
    sys.exit(main())
