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

Patched 2026-08-28 (session 38, Green Diamond MTT, FVS Method 3 OPEN 4
follow-up): STDINFO field 2 (habitat type / plant association code)
was left blank above with no note on what that does. For the IE
(Inland Empire) variant specifically, per the official FVS-IE Variant
Overview section 3.3, a blank or unrecognized habitat-type code is
NOT an error, FVS silently assigns the documented default, code 260
(PSME/PHMA). Section 3.3 also states there are 95 valid IE habitat
type codes, cross-walked to 30 original North Idaho codes in the
variant guide's Appendix A (table 11.1.1), and that plant-association
codes cross-walk to habitat codes via table 11.1.2.

This session inspected the one candidate field that might have
supplied a real per-stand habitat code, VEG_LBL on the Green Diamond
Montana property's MttStands2024.shp (69 unique values across all
10,851 stands, none blank, values like DF13/LP13/WL13/XX00/CX00).
That is a Green Diamond internal dominant-species-plus-size/stocking
cover-type label (2-letter species code + 2-digit size/stocking
suffix), not an IE habitat-type or plant-association code, and it
does not cross-walk to one by any documented FVS-IE method. No other
field in the shapefile schema is habitat/plant-association-shaped
either. Forcing VEG_LBL into field 2 would silently fabricate a
habitat code, so it was not done; field 2 stays blank here, which is
what makes IE fall through to its own documented 260 (PSME/PHMA)
default rather than this module guessing. `log_ie_habitat_default()`
below just makes that fallback explicit and auditable per stand
instead of silent, since the house rule against unexamined defaults
still applies even to an officially sanctioned one. If a real
per-stand habitat-type or plant-association source is ever found
(a Green Diamond soils/ecoclass layer, for instance), populate field
2 directly rather than extending VEG_LBL matching, since VEG_LBL
itself carries no productivity information to crosswalk from.

Patched 2026-09-06 (NCASI Model Evaluation Phase II, sessions 10 and
11): four defects found by the FVS-NE, FVS-LS and FVS-SN state-scope
runs on real FIA data and worked around at run time in each runner
are now fixed here so the workarounds become unnecessary.
  1. FIA_TO_FVS_SP["ne"] held 13 species and fia_to_fvs_code() fell
     back silently to the FIRST entry of the table on any unmapped
     SPCD, which on real Northeastern FIA data relabelled most stems
     as balsam fir. The NE table is now the variant's full crosswalk
     (107 mapped codes from a 108-slot table, slot 71 being blank in
     blkdat.f90), and the silent fallback is gone. fia_to_fvs_code()
     raises KeyError on an unmapped code, and make_inventory_keyfile()
     pre-filters tree_df to mapped codes, counts and drops the rest,
     and never relabels. Dropped stems are reported through the
     unmapped-stem counters below.
  2. STDINFO_DEFAULTS["ne"] used location code 902, which FVS-NE does
     not recognise (ne/blkdat.f90 IFORCD is 911, 919, 920, 921, 922).
     Now 919, White Mountain NF, the code the NE run used.
  3. There was no "ls" key in either table, so a Lake States request
     silently resolved the Pacific Northwest species table and
     location code 612. Added the 66-entry LS crosswalk (FIA SPCD 125,
     red pine, appears twice in blkdat.f90 as RN and RP; the later
     entry RP wins, both being the variant's red pine) and location
     code 904, Huron-Manistee NF. 904 maps to IFOR = 3 in ls/forkod.f90,
     which lands on CASE DEFAULT in every merchantability branch of
     ls/sitset.f90, so no single national forest's merchantability
     overrides are imposed, whereas 903 (Chippewa) and 909 (Superior)
     each carry such overrides.
  4. FIA_TO_FVS_SP["sn"] held 10 of 90 species and STDINFO_DEFAULTS
     ["sn"] used 803, which is not a valid five-digit Region 8
     region/forest/district code. sn/forkod.f90 computes IFORDI =
     KODFOR/100 = 8, matches no JFOR entry, calls ERRGRO and then reads
     IFOR before assignment. Now the full 90-entry SN crosswalk and
     80301 (Chattahoochee-Oconee NF, Armuchee District), which resolves
     to IFOR = 3 and lands on CASE DEFAULT in every sn/sitset.f90
     merchantability branch.
  Also, both lookups used .get(variant, <pn table>), so any unknown
  variant silently resolved the Pacific Northwest tables. Unknown
  variants now raise KeyError.
The crosswalk tables below are literal copies of each variant's
config/<variant>.json categories.species_definitions (FIAJSP[i] ->
JSP[i]), which is the blkdat.f90 table as extracted, and were checked
identical between config/ and config/calibrated/ on 2026-09-06.

Patched 2026-09-07 (NCASI Model Evaluation Phase II, session 13):
SPCTRN second-stage species resolution. FVS itself resolves an input
species code in two stages (base/intree.f90). Stage one is the
variant's FIAJSP table, which is what FIA_TO_FVS_SP reproduces. Stage
two, called whenever stage one misses, is SPCTRN in
src-converted/vls/spctrn.f90, a 562-row translation table (ASPT) whose
columns are FVS alpha code, FIA code, USDA PLANTS symbol, and the
target alpha code for the CS, LS, NE and SN variants respectively; a
code with no row at all lands on the variant's "other" slot (CS 85,
LS 49 OH, NE 98 OH, SN 90 OT). The state-scope runners pre-translated
to alpha codes through this module and so bypassed stage two, which
is why FVS-SN dropped 4.35 to 6.73 percent of Southeastern stems
(hickories recorded at species level, eastern redcedar 068, willow
oak 831) that the engine would have resolved: the SN FIAJSP table
keys several aggregate slots by FIA genus code (HI 400, JU 057) that
field crews never record. See active-projects/ncasi-modeleval/fvs-se/
se-unmapped-stems-diagnosis_DRAFT_2026-09-06.md.
  The table is parsed from the Fortran source at import time
  (SPCTRN_SOURCE, SPCTRN_ROW_COUNT), never retyped. Resolution order
  per stem is FIAJSP, then the matching variant's SPCTRN column, then
  drop with a count (or, with spctrn_other_fallback=True, the variant's
  own "other" slot, which is what the engine does). Only the matching
  variant's column is applied; the CS column is parsed but unused
  since this module has no CS table. Per-variant counts of stems
  resolved by each stage and dropped are kept in the same counters as
  before and are returned by get_species_resolution_counts().
  fia_to_fvs_code() is unchanged (FIAJSP only, raises on a miss);
  resolve_fvs_code() is the two-stage lookup. make_inventory_keyfile()
  takes use_spctrn (default False, so existing callers and tests see
  identical output) and the runners should pass use_spctrn=True.

Author: A. Weiskittel
Date: 2026-04-25
"""

from __future__ import annotations

import json
import os
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd


# Mapping from FIA SPCD to FVS 2-letter species code per variant.
# The pn and ie tables are still the short curated subsets that the
# bakuzis SPECIES_GROUPS_* definitions in
# bakuzis_uncertainty_comparison.py need. The ne, ls and sn tables are
# the full variant crosswalks, generated 2026-09-06 from
# config/<variant>.json categories.species_definitions (FIAJSP[i] ->
# JSP[i], entries kept in variant index order, blank slots skipped).
# Since 2026-09-06 there is no silent fallback for a code missing from
# any of these tables (see fia_to_fvs_code), so a short table now
# means dropped stems rather than relabelled ones.
FIA_TO_FVS_SP = {
    "pn": {
        15:  "WF", 17:  "GF", 19:  "AF", 98:  "SS",
        108: "LP", 116: "JP", 119: "WP", 122: "PP",
        202: "DF", 263: "WH",
    },
    # FVS-SN, all 90 species, from config/sn.json (blkdat.f90 table).
    "sn": {
        10: "FR",    57: "JU",    90: "PI",    107: "PU",
        110: "SP",   111: "SA",   115: "SR",   121: "LL",
        123: "TM",   126: "PP",   128: "PD",   129: "WP",
        131: "LP",   132: "VP",   221: "BY",   222: "PC",
        260: "HM",   311: "FM",   313: "BE",   316: "RM",
        317: "SV",   318: "SM",   330: "BU",   370: "BB",
        372: "SB",   391: "AH",   400: "HI",   450: "CA",
        460: "HB",   471: "RD",   491: "DW",   521: "PS",
        531: "AB",   540: "AS",   541: "WA",   543: "BA",
        544: "GA",   552: "HL",   555: "LB",   580: "HA",
        591: "HY",   601: "BN",   602: "WN",   611: "SU",
        621: "YP",   650: "MG",   651: "CT",   652: "MS",
        653: "MV",   654: "ML",   660: "AP",   680: "MB",
        691: "WT",   693: "BG",   694: "TS",   701: "HH",
        711: "SD",   721: "RA",   731: "SY",   740: "CW",
        743: "BT",   762: "BC",   802: "WO",   806: "SO",
        812: "SK",   813: "CB",   819: "TO",   820: "LK",
        822: "OV",   824: "BJ",   825: "SN",   826: "CK",
        827: "WK",   832: "CO",   833: "RO",   834: "QS",
        835: "PO",   837: "BO",   838: "LO",   901: "BK",
        920: "WI",   931: "SS",   950: "BD",   970: "EL",
        971: "WE",   972: "AE",   975: "RL",   299: "OS",
        998: "OH",   999: "OT",
    },
    "ie": {
        15:  "GF", 17:  "GF", 19:  "AF", 73:  "WL",
        93:  "ES", 108: "LP", 119: "WP", 122: "PP",
        202: "DF",
    },
    # FVS-NE, 107 mapped codes from the 108-slot blkdat.f90 table
    # (slot 71 is blank in the variant), from config/ne.json. Note
    # FIA 125 (red pine) is RN in this variant.
    "ne": {
        12: "BF",    71: "TA",    94: "WS",    97: "RS",
        91: "NS",    95: "BS",    90: "PI",    125: "RN",
        129: "WP",   131: "LP",   132: "VP",   241: "WC",
        43: "AW",    68: "RC",    57: "JU",    261: "EH",
        260: "HM",   100: "OP",   105: "JP",   110: "SP",
        123: "TM",   126: "PP",   128: "PD",   130: "SC",
        299: "OS",   316: "RM",   318: "SM",   314: "BM",
        317: "SV",   371: "YB",   372: "SB",   373: "RB",
        375: "PB",   379: "GB",   400: "HI",   403: "PH",
        405: "SL",   407: "SH",   409: "MH",   531: "AB",
        540: "AS",   541: "WA",   543: "BA",   544: "GA",
        545: "PA",   621: "YP",   611: "SU",   651: "CT",
        746: "QA",   741: "BP",   742: "EC",   743: "BT",
        744: "PY",   762: "BC",   802: "WO",   823: "BR",
        826: "CK",   835: "PO",   800: "OK",   806: "SO",
        817: "QI",   827: "WK",   830: "PN",   832: "CO",
        804: "SW",   825: "SN",   833: "RO",   812: "SK",
        837: "BO",   813: "CB",   330: "BU",   332: "YY",
        374: "WR",   462: "HK",   521: "PS",   591: "HY",
        601: "BN",   602: "WN",   641: "OO",   650: "MG",
        653: "MV",   660: "AP",   691: "WT",   693: "BG",
        711: "SD",   712: "PW",   731: "SY",   831: "WL",
        901: "BK",   922: "BL",   931: "SS",   951: "BW",
        952: "WB",   970: "EL",   972: "AE",   975: "RL",
        998: "OH",   313: "BE",   315: "ST",   341: "AI",
        356: "SE",   391: "AH",   491: "DW",   500: "HT",
        701: "HH",   760: "PL",   761: "PR",
    },
    # FVS-LS, 66 mapped codes from the 68-slot blkdat.f90 table, from
    # config/ls.json. One slot is blank and FIA 125 (red pine, Pinus
    # resinosa Aiton) appears twice, as RN then RP; the later entry RP
    # wins, and both are the variant's own red pine, so this is not a
    # substitution. Added 2026-09-06; before that there was no "ls" key
    # and a Lake States run silently resolved the pn table.
    "ls": {
        105: "JP",   130: "SC",   125: "RP",   129: "WP",
        94: "WS",    91: "NS",    12: "BF",    95: "BS",
        71: "TA",    241: "WC",   261: "EH",   299: "OS",
        68: "RC",    543: "BA",   544: "GA",   742: "EC",
        317: "SV",   316: "RM",   762: "BC",   972: "AE",
        975: "RL",   977: "RE",   371: "YB",   951: "BW",
        318: "SM",   314: "BM",   531: "AB",   541: "WA",
        802: "WO",   804: "SW",   823: "BR",   826: "CK",
        833: "RO",   837: "BO",   809: "NP",   402: "BH",
        403: "PH",   407: "SH",   743: "BT",   746: "QA",
        741: "BP",   375: "PB",   601: "BN",   602: "WN",
        701: "HH",   901: "BK",   998: "OH",   313: "BE",
        315: "ST",   319: "MM",   391: "AH",   421: "AC",
        462: "HK",   491: "DW",   500: "HT",   660: "AP",
        693: "BG",   731: "SY",   761: "PR",   763: "CC",
        760: "PL",   920: "WI",   922: "BL",   923: "DM",
        931: "SS",   935: "MA",
    },
}


# Standard variant default for STDINFO field 1 (forest location code).
# FVS uses this internally to pick parameter sets even without
# state/county. Each code must be one the variant's forkod.f90 actually
# recognises; an unrecognised code is not a benign default (see the
# 2026-09-06 docstring entry for what 902 and 803 did).
STDINFO_DEFAULTS = {
    "pn": {"forest": 612, "fortyp_default": 201},   # Willamette NF
    # 80301 = Chattahoochee-Oconee NF, Armuchee District, a valid
    # five-digit Region 8 region/forest/district code (IFORDI 803 is
    # JFOR(3) in sn/forkod.f90, district 1 passes the 803 filter,
    # IFOR = 3 hits CASE DEFAULT in sn/sitset.f90). Was 803, which
    # forkod.f90 rejects and then reads IFOR unassigned.
    "sn": {"forest": 80301, "fortyp_default": 161},
    "ie": {"forest": 110, "fortyp_default": 201},   # Idaho Panhandle NF
    # 919 = White Mountain NF, in ne/blkdat.f90 IFORCD (911, 919, 920,
    # 921, 922). Was 902, which FVS-NE does not recognise.
    "ne": {"forest": 919, "fortyp_default": 504},
    # 904 = Huron-Manistee NF, in ls/blkdat.f90 IFORCD (902, 903, 904,
    # 906, 907, 909, 910, 913); IFOR = 3 hits CASE DEFAULT in every
    # ls/sitset.f90 merchantability branch. Added 2026-09-06, with the
    # same fortyp default (901, aspen) the LS state-scope runner used.
    "ls": {"forest": 904, "fortyp_default": 901},
    "acd": {"forest": 902, "fortyp_default": 504},
}

# FVS-IE (Inland Empire) documented habitat-type default. Per section
# 3.3 of the official FVS-IE Variant Overview: "If the habitat type
# code is blank or not recognized, the default 260 (PSME/PHMA) will be
# assigned." This is a real, officially sanctioned FVS fallback, not a
# guess invented here -- but per house convention it still needs to be
# tracked, not silently relied on. Added 2026-08-28 (session 38, Green
# Diamond MTT), after confirming this project's own MttStands2024.shp
# has no field that cross-walks to an IE habitat-type or plant-
# association code (VEG_LBL is a Green Diamond internal species/size/
# stocking cover-type label, checked and ruled out this session).
IE_DEFAULT_HABITAT_CODE = 260  # PSME/PHMA

# Module-level counters so a batch driver (e.g. one call per stand)
# can report what fraction of stands actually got a real habitat-type
# code in STDINFO field 2 versus fell through to the IE default. Reset
# with reset_ie_habitat_default_counts() between separate batch runs.
_IE_HABITAT_DEFAULT_COUNTS = {"total_ie_stands": 0, "defaulted_stands": 0}


def reset_ie_habitat_default_counts() -> None:
    """Reset the IE habitat-default tracking counters before a batch run."""
    _IE_HABITAT_DEFAULT_COUNTS["total_ie_stands"] = 0
    _IE_HABITAT_DEFAULT_COUNTS["defaulted_stands"] = 0


def log_ie_habitat_default(stand_id: str, habitat_code: Optional[int]) -> None:
    """Record whether an IE stand got a real habitat code or the default.

    Call once per stand for variant == "ie". habitat_code is whatever
    would be written to STDINFO field 2; pass None (or leave the field
    blank, as this module currently does for every variant) to record
    a fall-through to IE_DEFAULT_HABITAT_CODE. Does not change any FVS
    behavior -- the 260 default already applies internally whether or
    not this is called -- it only makes the fallback auditable instead
    of silent, per the house rule against unexamined defaults.
    """
    _IE_HABITAT_DEFAULT_COUNTS["total_ie_stands"] += 1
    if habitat_code is None:
        _IE_HABITAT_DEFAULT_COUNTS["defaulted_stands"] += 1


def get_ie_habitat_default_summary() -> str:
    """Return a human-readable count/percentage of IE stands defaulted
    to habitat code 260 (PSME/PHMA) versus given a real code."""
    n = _IE_HABITAT_DEFAULT_COUNTS["total_ie_stands"]
    d = _IE_HABITAT_DEFAULT_COUNTS["defaulted_stands"]
    if n == 0:
        return "No IE stands logged yet (call log_ie_habitat_default per stand)."
    pct = 100.0 * d / n
    return (
        f"IE habitat type: {d} of {n} stands ({pct:.1f}%) defaulted to "
        f"code {IE_DEFAULT_HABITAT_CODE} (PSME/PHMA) per FVS-IE section "
        f"3.3; {n - d} stand(s) carried a real habitat/plant-association "
        f"code in STDINFO field 2."
    )


def species_table(variant: str) -> dict:
    """Return the FIA SPCD -> FVS alpha crosswalk for a variant.

    Raises KeyError for a variant with no table. Before 2026-09-06 an
    unknown variant silently resolved the Pacific Northwest table,
    which is how a Lake States request came to be one keystroke away
    from writing PN species codes into an LS keyfile.
    """
    v = variant.lower()
    if v not in FIA_TO_FVS_SP:
        raise KeyError(
            f"No FIA-to-FVS species crosswalk for variant {variant!r}; "
            f"known variants are {sorted(FIA_TO_FVS_SP)}. Add the variant's "
            "table from config/<variant>.json species_definitions rather "
            "than reusing another variant's."
        )
    return FIA_TO_FVS_SP[v]


def stdinfo_defaults(variant: str) -> dict:
    """Return the STDINFO location-code defaults for a variant.

    Raises KeyError for a variant with no entry, for the same reason as
    species_table(): a silent fall-through to the PN entry (612) hands
    the variant a location code its forkod.f90 does not recognise.
    """
    v = variant.lower()
    if v not in STDINFO_DEFAULTS:
        raise KeyError(
            f"No STDINFO_DEFAULTS entry for variant {variant!r}; known "
            f"variants are {sorted(STDINFO_DEFAULTS)}. Add a location code "
            "that the variant's forkod.f90 recognises."
        )
    return STDINFO_DEFAULTS[v]


def fia_to_fvs_code(spcd: int, variant: str) -> str:
    """Return the FVS 2-letter species code for an FIA species code.

    Raises KeyError when the code is not in the variant's table.
    Before 2026-09-06 this fell back silently to the first entry in
    the table, which on real FIA data relabelled every unmapped stem
    as that species (balsam fir in NE, shortleaf pine in SN). Callers
    that want drop-not-raise semantics should pre-filter with
    species_table(variant), as make_inventory_keyfile() now does.
    """
    table = species_table(variant)
    if spcd in table:
        return table[spcd]
    raise KeyError(
        f"FIA SPCD {spcd} is not in the {variant.lower()} variant species "
        f"crosswalk ({len(table)} entries); it must be dropped or mapped "
        "explicitly, never relabelled."
    )


# ---------------------------------------------------------------------
# SPCTRN: FVS's second-stage species translation table.
#
# src-converted/vls/spctrn.f90 (linked into FVSne, FVSls and FVSsn by
# src-converted/bin/FVS*_sourceList.txt) declares CHARACTER*8
# ASPT(562,7) and fills it from DATA statements, one row per line:
#   'FR ','010','ABIES   ','OS ','BF ','BF ','FR ',   !Abies
# column 1 FVS alpha code, 2 FIA code (three characters, zero padded,
# blank when the row is a PLANTS-symbol-only entry), 3 USDA PLANTS
# symbol, 4 CS target, 5 LS target, 6 NE target, 7 SN target. The
# routine walks rows in order and takes the first row whose column 1,
# 2 or 3 equals the input code, so first occurrence wins (there are no
# duplicate FIA codes in the table as of 2026-09-07). It then looks the
# target alpha up in the variant's NSP list; a target absent from the
# variant leaves ISPC1 = 0 and the 2026-05-16 fork guard drops the tree
# silently, so this module treats such a row as unresolved too. A code
# with no row at all becomes the variant's "other" index (CS 85, LS 49,
# NE 98, SN 90), which is the SPCTRN_OTHER_FALLBACK alpha below.
# ---------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parents[2]
SPCTRN_SOURCE = _REPO_ROOT / "src-converted" / "vls" / "spctrn.f90"
SPCTRN_ROW_COUNT = 562          # MAXASPT in spctrn.f90; parse asserts it
SPCTRN_FIA_ROW_COUNT = 291      # rows carrying a non-blank FIA code
SPCTRN_COLUMNS = ("alpha", "fia", "plants", "cs", "ls", "ne", "sn")
SPCTRN_VARIANT_COLUMNS = {"cs": "cs", "ls": "ls", "ne": "ne", "sn": "sn"}
# The alpha code of the "other" slot SPCTRN returns for a code with no
# table row, per variant (JSP index 85 CS, 49 LS = OH, 98 NE = OH,
# 90 SN = OT). CS is listed for completeness; there is no CS table here.
SPCTRN_OTHER_FALLBACK = {"ls": "OH", "ne": "OH", "sn": "OT"}
# Alpha codes the variant's JSP list carries that FIA_TO_FVS_SP cannot,
# because the same FIA code appears twice in FIAJSP and the later slot
# wins. LS slot 3 RN and slot 7 RP are both FIA 125 (red pine); the
# crosswalk keeps RP, but SPCTRN sends the southern pines (110, 111,
# 121, 125, 131, 136, 144) to RN, which the engine models at slot 3.
# Checked against config/<variant>.json on 2026-09-07: NE and SN have
# no such slot.
VARIANT_ALPHAS_NOT_IN_CROSSWALK = {"ls": {"RN"}}

_SPCTRN_ROW_RE = re.compile(
    r"^\s*'([^']*)','([^']*)','([^']*)','([^']*)','([^']*)','([^']*)','([^']*)'",
    re.M,
)


@lru_cache(maxsize=None)
def _parse_spctrn(path: str = str(SPCTRN_SOURCE)) -> tuple:
    """Parse the ASPT DATA rows out of spctrn.f90.

    Returns a tuple of 7-tuples in table order (see SPCTRN_COLUMNS),
    every field stripped. Asserts the row count equals SPCTRN_ROW_COUNT
    so a future edit to the Fortran that changes the table size fails
    loudly here rather than silently shifting the mapping.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    rows = tuple(tuple(f.strip() for f in m) for m in _SPCTRN_ROW_RE.findall(text))
    if len(rows) != SPCTRN_ROW_COUNT:
        raise RuntimeError(
            f"SPCTRN parse of {path} found {len(rows)} ASPT rows, expected "
            f"{SPCTRN_ROW_COUNT} (MAXASPT); the Fortran table changed or the "
            "regex missed rows. Update SPCTRN_ROW_COUNT deliberately."
        )
    return rows


@lru_cache(maxsize=None)
def spctrn_table(variant: str) -> dict:
    """Return the SPCTRN FIA SPCD -> target alpha map for one variant.

    Only rows with a numeric FIA code are kept (first occurrence wins,
    matching the Fortran loop), only the matching variant's target
    column is read, and a target alpha that the variant's own FIAJSP
    table does not carry is omitted because the engine would drop that
    tree at the ISPC1 guard rather than model it. Raises KeyError for a
    variant SPCTRN does not translate (anything but cs, ls, ne, sn).
    """
    v = variant.lower()
    if v not in SPCTRN_VARIANT_COLUMNS:
        raise KeyError(
            f"SPCTRN translates only {sorted(SPCTRN_VARIANT_COLUMNS)}; variant "
            f"{variant!r} has no target column in spctrn.f90."
        )
    col = SPCTRN_COLUMNS.index(SPCTRN_VARIANT_COLUMNS[v])
    known = None
    if v in FIA_TO_FVS_SP:
        known = set(species_table(v).values()) | VARIANT_ALPHAS_NOT_IN_CROSSWALK.get(v, set())
    out: dict = {}
    for row in _parse_spctrn():
        fia = row[1]
        if not fia.isdigit():
            continue
        code = int(fia)
        if code in out:
            continue                      # first row wins, as in the loop
        target = row[col]
        if not target or target == "XX":
            continue
        if known is not None and target not in known:
            continue                      # engine would hit the ISPC1 guard
        out[code] = target
    return out


def resolve_fvs_code(
    spcd: int,
    variant: str,
    use_spctrn: bool = True,
    other_fallback: bool = False,
) -> tuple:
    """Two-stage species resolution, the way base/intree.f90 does it.

    Returns (alpha, stage) where stage is "fiajsp" when the variant's
    own crosswalk carries the code, "spctrn" when the translation table
    supplied it, "other" when nothing matched and other_fallback is set
    (the variant's OH/OT slot, which is what the engine assigns), or
    None with alpha None when the stem must be dropped. FIAJSP always
    wins, so any code fia_to_fvs_code() already maps resolves exactly
    as before.
    """
    v = variant.lower()
    table = species_table(v)              # raises for an unknown variant
    spcd = int(spcd)
    if spcd in table:
        return table[spcd], "fiajsp"
    if use_spctrn and v in SPCTRN_VARIANT_COLUMNS:
        hit = spctrn_table(v).get(spcd)
        if hit is not None:
            return hit, "spctrn"
        if other_fallback and v in SPCTRN_OTHER_FALLBACK:
            return SPCTRN_OTHER_FALLBACK[v], "other"
    return None, None


# Module-level counters so a batch driver can report how many stems
# each resolution stage handled and how many were dropped, per variant,
# instead of that happening silently. Reset with
# reset_unmapped_stem_counts() between batch runs. Mirrors the per-state
# counts the state-scope runners log.
_UNMAPPED_STEM_COUNTS: dict = {}
_STAGE_KEYS = ("fiajsp", "spctrn", "other")


def reset_unmapped_stem_counts() -> None:
    """Reset the per-variant species resolution counters."""
    _UNMAPPED_STEM_COUNTS.clear()


def _new_counter() -> dict:
    return {
        "stands": 0,
        "total_stems": 0,
        "fiajsp_stems": 0,
        "spctrn_stems": 0,
        "other_stems": 0,
        "dropped_stems": 0,
        "spctrn_by_spcd": Counter(),   # SPCD -> stems the table rescued
        "dropped_by_spcd": Counter(),  # SPCD -> stems dropped
    }


def get_species_resolution_counts(variant: Optional[str] = None) -> dict:
    """Return the per-variant species resolution counts.

    With variant=None returns {variant: counts}; with a variant returns
    that variant's counts (an empty counter if nothing was logged).
    Each counts dict carries stands, total_stems, fiajsp_stems,
    spctrn_stems, other_stems, dropped_stems, and the Counters
    spctrn_by_spcd and dropped_by_spcd keyed by FIA SPCD. Copies are
    returned, so a caller cannot corrupt the running totals.
    """
    def _copy(c: dict) -> dict:
        return {k: (Counter(v) if isinstance(v, Counter) else v) for k, v in c.items()}
    if variant is None:
        return {v: _copy(c) for v, c in _UNMAPPED_STEM_COUNTS.items()}
    return _copy(_UNMAPPED_STEM_COUNTS.get(variant.lower(), _new_counter()))


def get_unmapped_stem_summary() -> str:
    """Return a human-readable per-variant count of stems resolved by
    each stage and dropped from keyfiles."""
    if not _UNMAPPED_STEM_COUNTS:
        return "No stands logged yet (make_inventory_keyfile records per call)."
    parts = []
    for v, c in sorted(_UNMAPPED_STEM_COUNTS.items()):
        n, d = c["total_stems"], c["dropped_stems"]
        pct = 100.0 * d / n if n else 0.0
        extra = ""
        if c["spctrn_stems"] or c["other_stems"]:
            extra = (
                f", {c['fiajsp_stems']} via FIAJSP, {c['spctrn_stems']} via "
                f"SPCTRN, {c['other_stems']} to the variant other slot"
            )
        parts.append(
            f"{v}: {d} of {n} stems ({pct:.2f}%) dropped as SPCD not in the "
            f"variant species list (never relabelled){extra}, across "
            f"{c['stands']} stand(s)"
        )
    return "; ".join(parts) + "."


def _record_resolution(variant: str, stages: pd.Series, spcds: pd.Series) -> None:
    c = _UNMAPPED_STEM_COUNTS.setdefault(variant, _new_counter())
    c["stands"] += 1
    c["total_stems"] += int(len(stages))
    for key in _STAGE_KEYS:
        c[f"{key}_stems"] += int((stages == key).sum())
    dropped = stages.isna()
    c["dropped_stems"] += int(dropped.sum())
    c["dropped_by_spcd"].update(spcds[dropped].astype(int).tolist())
    c["spctrn_by_spcd"].update(spcds[stages == "spctrn"].astype(int).tolist())


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
    use_spctrn: bool = False,
    spctrn_other_fallback: bool = False,
) -> str:
    """Generate an FVS keyfile in INVENTORY mode for the given stand.

    Returns a string containing the complete keyfile suitable for
    writing to disk and loading via fvs2py.

    The stand_df is expected to have at least: site_index, age,
    aspect, slope, elevft, forest_type. The tree_df is expected to
    have: tree_count (TPA), species (FIA SPCD), diameter, ht, crratio.

    Species resolution is FIAJSP first (FIA_TO_FVS_SP), then, when
    use_spctrn is True, the variant's column of FVS's own SPCTRN table
    (see resolve_fvs_code), then drop with a count. With
    spctrn_other_fallback=True a code with no SPCTRN row goes to the
    variant's OH/OT slot instead of being dropped, which is exactly what
    the engine does when handed the FIA code directly. The default,
    use_spctrn=False, reproduces the 2026-09-06 behaviour byte for byte.
    The state-scope runners should pass use_spctrn=True and read the
    per-variant counts back with get_species_resolution_counts().
    """
    v = variant.lower()
    s = stand_df.iloc[0]
    defaults = stdinfo_defaults(v)          # raises for an unknown variant
    forest = defaults["forest"]
    fortyp_default = defaults["fortyp_default"]
    fortyp = int(s.get("forest_type", fortyp_default) or fortyp_default)

    # Resolve every tree record's SPCD up front, count the stage that
    # resolved it, and drop the rest. Unresolved stems are counted and
    # dropped, never relabelled. A stand with no resolvable stems at all
    # is an error, not a keyfile.
    species_table(v)                        # raises for an unknown variant
    n_total = len(tree_df)
    spcds = tree_df["species"].astype(int)
    resolved = [
        resolve_fvs_code(c, v, use_spctrn=use_spctrn,
                         other_fallback=spctrn_other_fallback)
        for c in spcds
    ]
    alphas = pd.Series([r[0] for r in resolved], index=tree_df.index, dtype=object)
    stages = pd.Series([r[1] for r in resolved], index=tree_df.index, dtype=object)
    _record_resolution(v, stages, spcds)
    mapped_mask = alphas.notna()
    n_dropped = int(n_total - mapped_mask.sum())
    if n_dropped:
        tree_df = tree_df.loc[mapped_mask]
        alphas = alphas.loc[mapped_mask]
    if len(tree_df) == 0:
        raise ValueError(
            f"Stand {stand_id}: none of its {n_total} tree records carry an "
            f"FIA SPCD in the {v} variant crosswalk; refusing to write a "
            "keyfile with no trees."
        )

    age = float(s.get("age", 40))
    si = float(s.get("site_index", 60))
    aspect = float(s.get("aspect", 0))
    slope = float(s.get("slope", 15))
    elev_h = float(s.get("elevft", 1000)) / 100.0  # STDINFO wants elev/100

    # STDINFO fields are (forest, habitat/PV code, age, aspect, slope,
    # elev in 100s of feet). Field 2 (habitat/PV code) is populated from
    # stand_df["habitat_code"] when the caller supplies a real IE
    # habitat-type or plant-association code (int); otherwise it is left
    # blank, same as before 2026-08-28. Left blank, FVS-IE section 3.3's
    # own documented default (code 260, PSME/PHMA) applies internally.
    # No VEG_LBL-derived value is ever written here -- VEG_LBL was
    # checked this session and is a Green Diamond internal cover-type
    # label, not a habitat-type or plant-association code, so it is not
    # a legitimate source for this field. Age goes in field 3. Site
    # index is carried by a separate SITECODE keyword below.
    habitat_code = s.get("habitat_code", None)
    habitat_code = None if habitat_code is None or (isinstance(habitat_code, float) and pd.isna(habitat_code)) else int(habitat_code)
    habitat_field = f"{habitat_code:>10d}" if habitat_code is not None else " " * 10
    if v == "ie":
        log_ie_habitat_default(stand_id, habitat_code)
    stdinfo = (
        "STDINFO   "
        f"{forest:>10.0f}"
        + habitat_field
        + f"{age:>10.0f}"
        f"{aspect:>10.1f}"
        f"{slope:>10.1f}"
        f"{elev_h:>10.1f}"
    )

    # SITECODE: field 1 = site species (first tree record's species),
    # field 2 = site index, field 3 = 1 makes it the site species.
    site_sp = str(alphas.iloc[0])
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
    for tree_num, (row, sp_code) in enumerate(zip(tree_df.itertuples(), alphas), start=1):
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
    "resolve_fvs_code",
    "spctrn_table",
    "get_species_resolution_counts",
    "SPCTRN_SOURCE",
    "SPCTRN_ROW_COUNT",
    "SPCTRN_FIA_ROW_COUNT",
    "SPCTRN_OTHER_FALLBACK",
    "species_table",
    "stdinfo_defaults",
    "format_tree_record",
    "make_inventory_keyfile",
    "reset_unmapped_stem_counts",
    "get_unmapped_stem_summary",
    "log_ie_habitat_default",
    "get_ie_habitat_default_summary",
    "reset_ie_habitat_default_counts",
    "IE_DEFAULT_HABITAT_CODE",
    "FIA_TO_FVS_SP",
    "STDINFO_DEFAULTS",
]
