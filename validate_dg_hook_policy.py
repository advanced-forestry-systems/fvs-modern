#!/usr/bin/env python3
"""Validate the per-variant DG native-hook (DGDRIVER) policy.

Confirms:
  1. config_loader imports and the policy table round-trips (JSON load).
  2. All 23 variants present with the expected classification counts.
  3. A HOOK variant (ne) emits 'DGDRIVER 2'; a NATIVE variant (bm) emits none.
  4. generate_keywords is backward compatible: without apply_dg_hook_policy no
     DGDRIVER line appears; with it, ne gets one and bm does not.
"""
import sys, json, collections
from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parent / "config"
sys.path.insert(0, str(CONFIG_DIR))

import config_loader  # noqa: E402  (import test)
from config_loader import FvsConfigLoader

print("[1] config_loader imported OK")

# Round-trip the policy table.
table = FvsConfigLoader.load_dg_hook_policy_table(CONFIG_DIR)
variants = table["variants"]
counts = collections.Counter(v["classification"] for v in variants.values())
print(f"[2] policy table round-trips: {len(variants)} variants, {dict(counts)}")
assert len(variants) == 23, "expected 23 variants"
assert counts["HOOK"] == 7 and counts["NATIVE"] == 12 and counts["DEFER"] == 4

def kw_for(v):
    ld = FvsConfigLoader(v, version="calibrated", config_dir=CONFIG_DIR)
    return ld.dg_hook_policy(v), ld.dgdriver_keyword(v)

pol_ne, kw_ne = kw_for("ne")
pol_bm, kw_bm = kw_for("bm")
print(f"[3] ne  -> class={pol_ne['classification']:6s} dgdriver_keyword={kw_ne!r}")
print(f"    bm  -> class={pol_bm['classification']:6s} dgdriver_keyword={kw_bm!r}")
assert kw_ne is not None and kw_ne.split() == ["DGDRIVER", "2"], f"ne keyword wrong: {kw_ne!r}"
assert kw_bm is None, f"bm should emit no DGDRIVER, got {kw_bm!r}"

# Backward-compat check on generate_keywords.
ld_ne = FvsConfigLoader("ne", version="calibrated", config_dir=CONFIG_DIR)
kw_default = ld_ne.generate_keywords(include_comments=False)
kw_policy  = ld_ne.generate_keywords(include_comments=False, apply_dg_hook_policy=True)
assert "DGDRIVER" not in kw_default, "default generate_keywords must NOT emit DGDRIVER"
assert "DGDRIVER" in kw_policy, "opt-in generate_keywords must emit DGDRIVER for ne"
print("[4] generate_keywords: default has no DGDRIVER; opt-in emits it for ne")

ld_bm = FvsConfigLoader("bm", version="calibrated", config_dir=CONFIG_DIR)
kw_bm_policy = ld_bm.generate_keywords(include_comments=False, apply_dg_hook_policy=True)
assert "DGDRIVER" not in kw_bm_policy, "bm opt-in must still emit no DGDRIVER"
print("    generate_keywords: opt-in for bm still emits no DGDRIVER (native)")

print("\nALL CHECKS PASSED")
print("  ne (HOOK)   -> emitted keyword line: '%s'" % kw_ne)
print("  bm (NATIVE) -> no DGDRIVER keyword (native FVS DG)")
