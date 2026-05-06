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
STDINFO + INVYEAR + NUMCYCLE + TREEFMT + TREEDATA) with embedded tree
records in a simple known-good fixed-width format.

The TREEFMT used here is intentionally simple:
  (I4,T6,A2,T9,F4.1,T14,F3.0,T18,F2.0)
which encodes per record:
  cols 1-4:  tree number (I4)
  cols 6-7:  2-letter species code (A2, mapped from FIA SPCD)
  cols 9-12: DBH in inches with one decimal (F4.1)
  cols 14-16: total height in feet (F3.0)
  cols 18-19: live crown ratio percent (F2.0)

Author: A. Weiskittel
Date: 2026-04-25
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

import pandas as pd


# Mapping from FIA SPCD to FVS 2-letter species code per variant.
# These come from the calibrated JSON files' species_definitions.JSP
# array. We only need the species used by the bakuzis SPECIES_GROUPS_*
# definitions in bakuzis_uncertainty_comparison.py.
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
    "ne": {"forest": 902, "fortyp_default": 504},   # White Mountain NF
    "acd": {"forest": 902, "fortyp_default": 504},
}


def fia_to_fvs_code(spcd: int, variant: str) -> str:
    """Return the FVS 2-letter species code for an FIA species code.

    Falls back to the first entry in the variant's table if not found,
    so that synthetic stands always have valid records even for
    species codes that slipped past the SPECIES_GROUPS_* curation.
    """
    table = FIA_TO_FVS_SP.get(variant.lower(), FIA_TO_FVS_SP["pn"])
    if spcd in table:
        return table[spcd]
    # Fallback to the first species in the table
    return list(table.values())[0]


def format_tree_record(
    tree_num: int,
    sp_code: str,
    dbh: float,
    ht: float,
    cr: float,
) -> str:
    """Return one tree record matching TREEFMT (I4,T6,A2,T9,F4.1,T14,F3.0,T18,F2.0).

    Layout:
      cols 1-4:  tree number, right-justified
      col 5:     blank
      cols 6-7:  species code (left-justified within 2 chars)
      col 8:     blank
      cols 9-12: DBH like " 6.2" or "11.5"
      col 13:    blank
      cols 14-16: height like " 30" or "120"
      col 17:    blank
      cols 18-19: CR percent like "60"
    """
    sp = (sp_code + "  ")[:2]
    cr_int = max(0, min(99, int(round(cr))))
    line = (
        f"{tree_num:4d}"               # cols 1-4: I4
        f" "                            # col 5
        f"{sp:2s}"                      # cols 6-7
        f" "                            # col 8
        f"{dbh:4.1f}"                   # cols 9-12: F4.1
        f" "                            # col 13
        f"{ht:3.0f}"                    # cols 14-16: F3.0
        f" "                            # col 17
        f"{cr_int:2d}"                  # cols 18-19: F2.0 as int
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

    # STDINFO line: 6 fields, 10 chars each, fixed width
    stdinfo = (
        "STDINFO   "
        f"{forest:>10.0f}"
        f"{age:>10.0f}"
        f"{si:>10.1f}"
        f"{aspect:>10.1f}"
        f"{slope:>10.1f}"
        f"{elev_h:>10.1f}"
    )

    lines = []
    lines.append("STDIDENT")
    lines.append(stand_id)
    lines.append("")
    lines.append("DESIGN                                        11.0       1.0")
    lines.append(stdinfo)
    lines.append(f"INVYEAR     {float(inv_year):>10.1f}")
    lines.append(f"NUMCYCLE    {float(num_cycles):>10.1f}")
    lines.append("TIMEINT            0         5")
    lines.append("ECHOSUM")
    lines.append("CALBSTAT")

    # Calibration keywords (GROWMULT/MORTMULT/SDIMAX/HTGMULT) come
    # before TREEDATA so they are parsed during keyword phase.
    if calibration_keywords:
        lines.append(calibration_keywords)

    # TREEFMT specifies the column layout of the embedded TREEDATA records.
    lines.append("TREEFMT          (I4,T6,A2,T9,F4.1,T14,F3.0,T18,F2.0)")
    lines.append("TREEDATA")

    # Tree records, one per row of tree_df.
    for tree_num, row in enumerate(tree_df.itertuples(), start=1):
        spcd = int(row.species)
        sp_code = fia_to_fvs_code(spcd, v)
        # tree_count is TPA per record; each record is one expanded tree
        lines.append(
            format_tree_record(
                tree_num=tree_num,
                sp_code=sp_code,
                dbh=float(row.diameter),
                ht=float(row.ht),
                cr=float(row.crratio),
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
