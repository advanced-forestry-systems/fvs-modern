"""
Inventory-mode keyfile generator for FVS-PN/SN/IE/other non-eastern variants.

The DATABASE/STANDSQL keyword path used by the original Bakuzis runner
works for FVS-NE and FVS-ACD because those variants' DBS readers
expect Maine-style state/county fields. Western and southern variants
(PN, SN, IE, etc.) read a different set of fields from the same
fvs_standinit table, and the synthetic generator's output is rejected
with a Fortran EOF error in errgro.f90 line 55.

This module bypasses the DATABASE path entirely. It generates a
self-contained keyfile that uses INVENTORY mode (STDIDENT + DESIGN +
STDINFO + SITECODE + INVYEAR + NUMCYCLE + TREEDATA) with embedded
tree records in the FVS DEFAULT tree record layout (the TREFMT
initialized in each variant's blkdat.f90 and read in
base/intree.f90), so no custom TREEFMT keyword is emitted:
  (I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,...)
  cols 1-4:   plot ID
  cols 5-7:   tree number (tree ID is cols 1-7 combined)
  cols 8-13:  PROB (trees per acre this record represents, given the
              per-acre DESIGN emitted below)
  col 14:     history code (1 = live)
  cols 15-17: species code (mapped from FIA SPCD)
  cols 18-21: DBH in inches (F4.1)
  cols 25-27: total height in feet (F3.0)
  col 35:     crown ratio code 1-9 (10 percent classes)

Patched 2026-08-27: TREEFMT/TREEDATA/STDINFO/SITECODE/DESIGN fixes
after PN/SN smoke test. Specifically: (1) dropped the custom TREEFMT
keyword, whose single-line form consumed downstream records as format
text, in favor of the default layout; (2) TREEDATA now carries 15 in
field 1 so tree records are read from the keyword file rather than an
external .tre dataset; (3) tree records follow the default variable
order (plot, tree, PROB, history, species, DBH, DG, HT, THT, HTG, ICR
code); (4) STDINFO field 3 now carries stand age (field 2, habitat/PV
code, is left blank) and site index moved to a SITECODE keyword; (5)
DESIGN now specifies a single 1-acre fixed plot (expansion 1.0) so
PROB values are trees per acre, replacing the 11-plot BAF-40 prism
design copied from test decks that inflated TPA roughly ninefold.

Author: A. Weiskittel
Date: 2026-04-25
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd


# Mapping from FIA SPCD to FVS 2-letter species code per variant.
# These come from the calibrated JSON files' species_definitions.JSP
# array. We only need the species used by the bakuzis SPECIES_GROUPS_*
# definitions in bakuzis_uncertainty_comparison.py.
_CONFIG_DIR = Path(__file__).resolve().parents[2] / "config" / "calibrated"
_XW_CACHE: dict = {}

FIA_TO_FVS_SP = {
    "pn": {
        15:  "WF", 17:  "GF", 19:  "AF", 98:  "SS",
        108: "LP", 116: "JP", 119: "WP", 122: "PP",
        202: "DF", 263: "WH",
    },
    "sn": {
        110: "SP", 111: "SA", 121: "LL", 129: "WP",
        131: "LP", 132: "VP", 221: "BY", 311: "FM",
        313: "BE", 316: "RM",
    },
    "ie": {
        15:  "GF", 17:  "GF", 19:  "AF", 73:  "WL",
        93:  "ES", 108: "LP", 119: "WP", 122: "PP",
        202: "DF",
    },
    # Eastern fallback (kept for API symmetry; eastern variants use
    # DATABASE path and do not call this module)
    "ne": {
        12:  "BF", 91:  "NS", 95:  "BS", 97:  "RS",
        129: "WP", 241: "WA", 261: "EH", 316: "RM",
        318: "SM", 371: "YB", 531: "AB", 802: "WO",
        833: "RO",
    },
}


# Standard variant default for STDINFO field 1 (forest code). FVS uses
# this internally to pick parameter sets even without state/county.
# These are well-known representative national forests for each variant.
STDINFO_DEFAULTS = {
    "pn": {"forest": 612, "fortyp_default": 201},   # Willamette NF
    "sn": {"forest": 803, "fortyp_default": 161},   # Talladega NF
    "ie": {"forest": 110, "fortyp_default": 201},   # Idaho Panhandle NF
    "ne": {"forest": 919, "fortyp_default": 504},   # White Mountain NF (919 per ne/blkdat.f90; 902 unrecognised)
    "ls": {"forest": 904, "fortyp_default": 805},   # Huron-Manistee NF; northern hardwoods
    "acd": {"forest": 919, "fortyp_default": 504},
}


def _variant_crosswalk(variant: str) -> dict:
    """FIA SPCD -> FVS 2-letter code for a variant, read from the variant's
    calibrated config (categories/species_definitions FIAJSP -> JSP), which is
    the full species list the compiled variant recognises. Falls back to the
    hand table above only when the config is absent. Cached per variant."""
    v = variant.lower()
    if v in _XW_CACHE:
        return _XW_CACHE[v]
    table = None
    cfg = _CONFIG_DIR / f"{v}.json"
    if cfg.exists():
        try:
            sd = json.loads(cfg.read_text())["categories"]["species_definitions"]
            table = {}
            for raw, code in zip(sd["FIAJSP"], sd["JSP"]):
                raw = str(raw).strip(); code = str(code).strip()
                if not raw.isdigit() or not code:        # blank slots exist in several variants
                    continue
                spcd = int(raw)
                if spcd > 0 and spcd not in table:      # first mapping wins (e.g. ls 125 RN/RP)
                    table[spcd] = code
        except Exception as exc:  # pragma: no cover
            warnings.warn(f"inventory_keyfile: could not read crosswalk from {cfg}: {exc}")
            table = None
    if table is None:
        table = FIA_TO_FVS_SP.get(v, FIA_TO_FVS_SP["pn"])
    _XW_CACHE[v] = table
    return table


def fia_to_fvs_code(spcd: int, variant: str, strict: bool = True):
    """Return the FVS 2-letter species code for an FIA species code.

    Reads the variant's full crosswalk from config/calibrated/<variant>.json.
    With strict=True (default) an unmapped SPCD raises KeyError so callers
    pre-filter stems rather than silently relabelling them; strict=False
    returns None. The pre-September-2026 behaviour (fall back to the first
    entry in a 10 to 13 species hand table, which on real Northeast data
    relabelled most stems as balsam fir) is gone.
    """
    table = _variant_crosswalk(variant)
    if spcd in table:
        return table[spcd]
    if strict:
        raise KeyError(f"FIA SPCD {spcd} has no FVS code in variant {variant!r} crosswalk ({len(table)} entries)")
    return None


def filter_mapped_species(spcds, variant: str):
    """Return (kept_index_list, dropped_spcd_counter) for a sequence of SPCDs."""
    from collections import Counter
    table = _variant_crosswalk(variant)
    keep, dropped = [], Counter()
    for i, s in enumerate(spcds):
        if s in table:
            keep.append(i)
        else:
            dropped[s] += 1
    return keep, dropped


def format_tree_record(
    tree_num: int,
    sp_code: str,
    dbh: float,
    ht: float,
    cr: float,
    prob: float = 1.0,
    plot: int = 1,
) -> str:
    """Return one tree record in the FVS default TREEFMT layout.

    Default format per variant blkdat.f90 / base/intree.f90:
      (I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,...)

    Layout written here:
      cols 1-4:   plot ID
      cols 5-7:   tree number within plot
      cols 8-13:  PROB with explicit decimal (TPA under the per-acre
                  DESIGN emitted by make_inventory_keyfile)
      col 14:     history code, 1 = live
      cols 15-17: species code, left-justified
      cols 18-21: DBH like " 6.2" or "11.5"
      cols 22-24: DG, blank
      cols 25-27: height like " 30" or "120"
      cols 28-34: THT and HTG, blank
      col 35:     crown ratio code 1-9 (code n = ((n-1)*10, n*10] pct)
    """
    sp = (sp_code.strip() + "   ")[:3]
    icr = max(1, min(9, (int(cr) + 9) // 10))
    line = (
        f"{plot:4d}"                    # cols 1-4: plot ID
        f"{tree_num:3d}"                # cols 5-7: tree number
        f"{prob:6.1f}"                  # cols 8-13: PROB (F6.0)
        f"1"                            # col 14: history = live
        f"{sp:3s}"                      # cols 15-17: species (A3)
        f"{dbh:4.1f}"                   # cols 18-21: DBH (F4.1)
        f"   "                          # cols 22-24: DG blank
        f"{ht:3.0f}"                    # cols 25-27: HT (F3.0)
        f"       "                      # cols 28-34: THT, HTG blank
        f"{icr:1d}"                     # col 35: ICR code (I1)
    )
    return line


def make_inventory_keyfile(
    stand_df: pd.DataFrame,
    tree_df: pd.DataFrame,
    stand_id: str,
    variant: str,
    inv_year: int = 2000,
    num_cycles: int = 20,
    calibration_keywords: str = "",
) -> str:
    """Generate an FVS keyfile in INVENTORY mode for the given stand.

    Returns a string containing the complete keyfile suitable for
    writing to disk and loading via fvs2py.

    The stand_df is expected to have at least: site_index, age,
    aspect, slope, elevft, forest_type. The tree_df is expected to
    have: tree_count (TPA), species (FIA SPCD), diameter, ht, crratio.
    """
    v = variant.lower()
    s = stand_df.iloc[0]
    forest = STDINFO_DEFAULTS.get(v, STDINFO_DEFAULTS["pn"])["forest"]
    fortyp_default = STDINFO_DEFAULTS.get(v, STDINFO_DEFAULTS["pn"])["fortyp_default"]
    fortyp = int(s.get("forest_type", fortyp_default) or fortyp_default)

    age = float(s.get("age", 40))
    si = float(s.get("site_index", 60))
    aspect = float(s.get("aspect", 0))
    slope = float(s.get("slope", 15))
    elev_h = float(s.get("elevft", 1000)) / 100.0  # STDINFO wants elev/100

    # STDINFO fields are (forest, habitat/PV code, age, aspect, slope,
    # elev in 100s of feet). Field 2 is left blank; age goes in field 3.
    # Site index is carried by a separate SITECODE keyword below.
    stdinfo = (
        "STDINFO   "
        f"{forest:>10.0f}"
        + " " * 10
        + f"{age:>10.0f}"
        f"{aspect:>10.1f}"
        f"{slope:>10.1f}"
        f"{elev_h:>10.1f}"
    )

    # SITECODE: field 1 = site species (first tree record's species),
    # field 2 = site index, field 3 = 1 makes it the site species.
    site_sp = fia_to_fvs_code(int(tree_df.iloc[0]["species"]), v)
    sitecode = (
        "SITECODE  "
        f"{site_sp:>10s}"
        f"{si:>10.1f}"
        f"{1.0:>10.1f}"
    )

    lines = []
    lines.append("STDIDENT")
    lines.append(stand_id)
    lines.append("")
    # DESIGN: BAF = -1 means large trees on a 1-acre fixed plot,
    # small-tree plot 1 acre, 1 inventory point, 0 nonstocked, so each
    # record's PROB is trees per acre with expansion 1.0 (base/notre.f90).
    lines.append(
        "DESIGN    "
        f"{-1.0:>10.1f}"
        f"{1.0:>10.1f}"
        f"{999.0:>10.1f}"
        f"{1.0:>10.1f}"
        f"{0.0:>10.1f}"
    )
    lines.append(stdinfo)
    lines.append(sitecode)
    lines.append(f"INVYEAR     {float(inv_year):>10.1f}")
    lines.append(f"NUMCYCLE    {float(num_cycles):>10.1f}")
    lines.append("TIMEINT            0         5")
    lines.append("ECHOSUM")
    lines.append("CALBSTAT")

    # Calibration keywords (GROWMULT/MORTMULT/SDIMAX/HTGMULT) come
    # before TREEDATA so they are parsed during keyword phase.
    if calibration_keywords:
        lines.append(calibration_keywords)

    # Tree records use the FVS default TREEFMT (no TREEFMT keyword).
    # TREEDATA field 1 = 15 reads the records from the keyword file
    # itself; blank would default to dataset unit 2 (external .tre).
    lines.append("TREEDATA  " + f"{15.0:>10.0f}")

    # Tree records, one per row of tree_df. tree_count is the TPA the
    # record represents and passes through PROB unchanged under the
    # per-acre DESIGN above.
    for tree_num, row in enumerate(tree_df.itertuples(), start=1):
        spcd = int(row.species)
        sp_code = fia_to_fvs_code(spcd, v)
        tpa = getattr(row, "tree_count", None)
        tpa = 1.0 if tpa is None or pd.isna(tpa) else float(tpa)
        lines.append(
            format_tree_record(
                tree_num=tree_num,
                sp_code=sp_code,
                dbh=float(row.diameter),
                ht=float(row.ht),
                cr=float(row.crratio),
                prob=tpa,
            )
        )

    lines.append("-999")  # End of TREEDATA marker
    lines.append("")
    lines.append("PROCESS")
    lines.append("STOP")

    return "\n".join(lines) + "\n"


__all__ = [
    "fia_to_fvs_code",
    "format_tree_record",
    "make_inventory_keyfile",
    "FIA_TO_FVS_SP",
    "STDINFO_DEFAULTS",
]
