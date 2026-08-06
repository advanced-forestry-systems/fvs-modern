#!/usr/bin/env python3
"""
fvs_input_guards.py -- DIAGNOSTIC-ONLY input-side guards for the FVS Northeast driver.

Two defects motivate this module.  Neither changes FVS model output on valid
input; both are driver-side and diagnostic.

DEFECT A -- invalid species codes are swallowed.
    FIA SPCD 11, 20, 263 and 264 are real FIA species that the Northeast
    variant does not carry, and 9999 is not an FIA code at all.  All five are
    silently remapped by base/intree.f90 to the internal code OH (other
    hardwood), announced only as a NOTE in the listing:
        NOTE: INPUT SPECIES CODE (263     )WAS SET TO (OH )
    (format statement at src-converted/base/intree.f90:421).  FVS then returns
    success and all five produce bit-identical output.  Genuinely wrong input
    passes without complaint.

DEFECT B -- StandID aliasing.
    The legacy driver built  sid = f"{variant.upper()}_{cn % 10_000_000:07d}",
    which collapses 3,321 FIA control numbers onto 50 distinct StandIDs.  It is
    latent today (every run_stand call gets its own tempfile.mkdtemp database
    and nothing downstream keys on StandID) but is a trap for any later join.
    The replacement is a dense sequential index plus a CN lookup table.  The
    full FIA CN cannot be used as an identifier inside FVS because
    src-converted/dbsqlite/dbstreesin.f90 declares
        INTEGER ITREE,IPLOT,HISTORY,CRWNR,DMG1,DMG2,DMG3,SVR1,SVR2   (line 20)
    and binds TREE_ID with fsql3_colint at line 97, i.e. a default 4-byte
    signed integer with maximum 2147483647; FIA CNs are 14 to 15 digits.

Dependency-light: standard library only.  pandas is optional and used only if
already importable.

--------------------------------------------------------------------------
SPECIES-LIST PROVENANCE  (derived from source, not from memory)
--------------------------------------------------------------------------
NE_SPCD / NE_SPCD_TO_ALPHA
    Source file : /users/PUOM0008/crsfaaron/fvs-modern/src-converted/ne/blkdat.f90
    DATA JSP    : lines 126-142  (108 internal alpha codes, MAXSP = 108)
    DATA FIAJSP : lines 144-160  (the parallel FIA SPCD crosswalk)
    MAXSP = 108 declared at src-converted/ne/common/PRGPRM.f90:12
    Slot 71 is blank in both arrays (an unused reserved slot), so the NE
    variant carries 107 distinct FIA species codes.
    Extracted 2026-08-06.

KNOWN_FIA_SPCD
    Union of the DATA FIAJSP arrays of every FVS variant blkdat.f90 in
    /users/PUOM0008/crsfaaron/fvs-modern/src-converted/*/blkdat.f90
    (24 variants: acd adk ak bm ca ci cr cs ec em ie kt ls nc ne oc op pn sn
    so tt ut wc ws).  248 distinct codes.  This is the "is it a real FIA
    species code" universe.  It is a lower bound on the true FIA SPCD list,
    which is why codes outside it but inside the structural FIA range
    (1..999) are reported as OUT_OF_REGION rather than INVALID.
    Extracted 2026-08-06.
"""
from __future__ import annotations

import csv
import os
import warnings
from collections import Counter, OrderedDict

__all__ = [
    "NE_SPCD", "NE_SPCD_TO_ALPHA", "KNOWN_FIA_SPCD",
    "VALID", "OUT_OF_REGION", "INVALID",
    "SpeciesCodeError", "OutOfRegionSpeciesWarning",
    "classify_spcd", "validate_spcd", "validate_spcd_series",
    "StandIDMap", "build_standid_map", "check_standid_collisions",
    "write_standid_lookup", "allocate_tree_id",
    "TREE_ID_MAX_SIGNED32", "MAX_TREES_PER_STAND",
]

SPCD_SOURCE_FILE = ("/users/PUOM0008/crsfaaron/fvs-modern/src-converted/ne/"
                    "blkdat.f90")
SPCD_SOURCE_LINES = "JSP 126-142; FIAJSP 144-160"

# --------------------------------------------------------------------------
# NE species table, in JSP index order.  ('' marks the blank slot 71.)
# Transcribed verbatim from ne/blkdat.f90 DATA JSP / DATA FIAJSP.
# --------------------------------------------------------------------------
_NE_TABLE = (
    (12, 'BF'),  (71, 'TA'),  (94, 'WS'),  (97, 'RS'),  (91, 'NS'),  (95, 'BS'),  (90, 'PI'),
    (125, 'RN'), (129, 'WP'), (131, 'LP'), (132, 'VP'), (241, 'WC'), (43, 'AW'),  (68, 'RC'),
    (57, 'JU'),  (261, 'EH'), (260, 'HM'), (100, 'OP'), (105, 'JP'), (110, 'SP'), (123, 'TM'),
    (126, 'PP'), (128, 'PD'), (130, 'SC'), (299, 'OS'), (316, 'RM'), (318, 'SM'), (314, 'BM'),
    (317, 'SV'), (371, 'YB'), (372, 'SB'), (373, 'RB'), (375, 'PB'), (379, 'GB'), (400, 'HI'),
    (403, 'PH'), (405, 'SL'), (407, 'SH'), (409, 'MH'), (531, 'AB'), (540, 'AS'), (541, 'WA'),
    (543, 'BA'), (544, 'GA'), (545, 'PA'), (621, 'YP'), (611, 'SU'), (651, 'CT'), (746, 'QA'),
    (741, 'BP'), (742, 'EC'), (743, 'BT'), (744, 'PY'), (762, 'BC'), (802, 'WO'), (823, 'BR'),
    (826, 'CK'), (835, 'PO'), (800, 'OK'), (806, 'SO'), (817, 'QI'), (827, 'WK'), (830, 'PN'),
    (832, 'CO'), (804, 'SW'), (825, 'SN'), (833, 'RO'), (812, 'SK'), (837, 'BO'), (813, 'CB'),
    (None, ''),  (330, 'BU'), (332, 'YY'), (374, 'WR'), (462, 'HK'), (521, 'PS'), (591, 'HY'),
    (601, 'BN'), (602, 'WN'), (641, 'OO'), (650, 'MG'), (653, 'MV'), (660, 'AP'), (691, 'WT'),
    (693, 'BG'), (711, 'SD'), (712, 'PW'), (731, 'SY'), (831, 'WL'), (901, 'BK'), (922, 'BL'),
    (931, 'SS'), (951, 'BW'), (952, 'WB'), (970, 'EL'), (972, 'AE'), (975, 'RL'), (998, 'OH'),
    (313, 'BE'), (315, 'ST'), (341, 'AI'), (356, 'SE'), (391, 'AH'), (491, 'DW'), (500, 'HT'),
    (701, 'HH'), (760, 'PL'), (761, 'PR'),
)

NE_SPCD_TO_ALPHA = OrderedDict((c, a) for (c, a) in _NE_TABLE if c is not None)
NE_SPCD = frozenset(NE_SPCD_TO_ALPHA)

# Internal code every unmatched input species is silently remapped to.
FALLBACK_ALPHA = "OH"

# --------------------------------------------------------------------------
# Union of DATA FIAJSP over all 24 variant blkdat.f90 files (see provenance).
# --------------------------------------------------------------------------
KNOWN_FIA_SPCD = frozenset((
    10, 11, 12, 15, 17, 18, 19, 20, 21, 22, 41, 42, 43, 57, 62, 63, 64, 65, 66, 68,
    69, 71, 72, 73, 81, 90, 91, 92, 93, 94, 95, 96, 97, 98, 100, 101, 102, 103, 104, 105,
    106, 107, 108, 109, 110, 111, 113, 114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 137, 142, 143, 201, 202, 211, 212, 221, 222, 231, 241, 242, 251,
    260, 261, 263, 264, 298, 299, 311, 312, 313, 314, 315, 316, 317, 318, 319, 321, 322, 324, 330, 331,
    332, 333, 341, 350, 351, 352, 356, 361, 370, 371, 372, 373, 374, 375, 376, 379, 391, 400, 401, 402,
    403, 404, 405, 407, 408, 409, 421, 431, 450, 460, 461, 462, 471, 475, 478, 491, 492, 500, 521, 531,
    540, 541, 542, 543, 544, 545, 546, 552, 555, 571, 580, 591, 600, 601, 602, 611, 621, 631, 641, 650,
    651, 652, 653, 654, 660, 680, 690, 691, 693, 694, 701, 711, 712, 721, 730, 731, 740, 741, 742, 743,
    744, 745, 746, 747, 748, 749, 760, 761, 762, 763, 768, 800, 801, 802, 803, 804, 805, 806, 807, 809,
    810, 811, 812, 813, 814, 815, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 830, 831,
    832, 833, 834, 835, 836, 837, 838, 839, 843, 901, 920, 922, 923, 928, 931, 935, 950, 951, 952, 970,
    971, 972, 974, 975, 977, 981, 998, 999,
))

# Structural range of an FIA SPCD.  FIA species codes are 1..999 (genus-level
# and unknown codes included); anything outside cannot be an FIA SPCD.
FIA_SPCD_MIN, FIA_SPCD_MAX = 1, 999

VALID = "VALID"
OUT_OF_REGION = "OUT_OF_REGION"
INVALID = "INVALID"


class SpeciesCodeError(ValueError):
    """Raised when an input species code cannot be an FIA SPCD."""


class OutOfRegionSpeciesWarning(UserWarning):
    """Real FIA species the NE variant does not carry; FVS will remap it to OH."""


# ==========================================================================
# Species validation
# ==========================================================================
def _coerce(spcd):
    """Return int(spcd) or None if it cannot be an integer species code."""
    if spcd is None:
        return None
    if isinstance(spcd, bool):
        return None
    if isinstance(spcd, int):
        return spcd
    if isinstance(spcd, float):
        return int(spcd) if float(spcd).is_integer() else None
    try:
        s = str(spcd).strip()
        return int(s) if s else None
    except (TypeError, ValueError):
        return None


def classify_spcd(spcd):
    """Non-raising primitive.

    Returns dict(code, cls, alpha, reason).
      VALID          in the NE variant FIAJSP table
      OUT_OF_REGION  a real FIA species code the NE variant does not carry
                     (or a code inside the structural FIA range that no FVS
                     variant table lists); FVS will silently remap it to OH
      INVALID        cannot be an FIA SPCD at all
    """
    code = _coerce(spcd)
    if code is None:
        return {"code": spcd, "cls": INVALID, "alpha": None,
                "reason": "not an integer species code"}
    if code in NE_SPCD:
        return {"code": code, "cls": VALID, "alpha": NE_SPCD_TO_ALPHA[code],
                "reason": "in NE variant FIAJSP table"}
    if not (FIA_SPCD_MIN <= code <= FIA_SPCD_MAX):
        return {"code": code, "cls": INVALID, "alpha": None,
                "reason": "outside the structural FIA SPCD range %d-%d"
                          % (FIA_SPCD_MIN, FIA_SPCD_MAX)}
    if code in KNOWN_FIA_SPCD:
        return {"code": code, "cls": OUT_OF_REGION, "alpha": FALLBACK_ALPHA,
                "reason": "real FIA species carried by another FVS variant but "
                          "not by NE; FVS will silently set it to OH"}
    return {"code": code, "cls": OUT_OF_REGION, "alpha": FALLBACK_ALPHA,
            "reason": "inside the FIA SPCD range but listed by no FVS variant "
                      "table; FVS will silently set it to OH"}


def validate_spcd(spcd, strict=False, warn=True):
    """Validate one species code.

    Default: raise SpeciesCodeError on INVALID, warn on OUT_OF_REGION.
    strict=True: raise on OUT_OF_REGION as well.
    Returns the classify_spcd() dict on success.
    """
    r = classify_spcd(spcd)
    if r["cls"] == INVALID:
        raise SpeciesCodeError(
            "FVS input guard: species code %r is not a valid FIA SPCD (%s). "
            "The FVS NE variant would silently remap it to %s and return "
            "success, so this input must be rejected before FVS is launched."
            % (spcd, r["reason"], FALLBACK_ALPHA))
    if r["cls"] == OUT_OF_REGION:
        msg = ("FVS input guard: species code %d is out of region for the NE "
               "variant (%s)." % (r["code"], r["reason"]))
        if strict:
            raise SpeciesCodeError(msg + "  strict=True.")
        if warn:
            warnings.warn(msg, OutOfRegionSpeciesWarning, stacklevel=2)
    return r


def validate_spcd_series(codes, strict=False, warn=True, label=""):
    """Validate an iterable of species codes.

    Returns a structured summary a driver can log:
        {n, n_valid, n_out_of_region, n_invalid, counts_by_class,
         valid_codes, out_of_region_codes, invalid_codes, details, label}
    Raising behaviour matches validate_spcd, but the WHOLE series is
    classified first so the error names every offending code at once.
    """
    codes = list(codes)
    details = [classify_spcd(c) for c in codes]
    by_class = Counter(d["cls"] for d in details)
    oor = Counter(d["code"] for d in details if d["cls"] == OUT_OF_REGION)
    bad = Counter(str(d["code"]) for d in details if d["cls"] == INVALID)
    good = Counter(d["code"] for d in details if d["cls"] == VALID)
    summary = {
        "label": label,
        "n": len(codes),
        "n_valid": by_class.get(VALID, 0),
        "n_out_of_region": by_class.get(OUT_OF_REGION, 0),
        "n_invalid": by_class.get(INVALID, 0),
        "counts_by_class": dict(by_class),
        "valid_codes": dict(good),
        "out_of_region_codes": dict(oor),
        "invalid_codes": dict(bad),
        "details": details,
    }
    if summary["n_invalid"]:
        raise SpeciesCodeError(
            "FVS input guard%s: %d record(s) carry species codes that are not "
            "valid FIA SPCDs: %s.  FVS NE would silently remap every one of "
            "them to %s and return success."
            % (" [%s]" % label if label else "", summary["n_invalid"],
               ", ".join("%s x%d" % (k, v) for k, v in sorted(bad.items())),
               FALLBACK_ALPHA))
    if summary["n_out_of_region"]:
        msg = ("FVS input guard%s: %d record(s) carry FIA species codes the NE "
               "variant does not carry and will silently remap to %s: %s."
               % (" [%s]" % label if label else "", summary["n_out_of_region"],
                  FALLBACK_ALPHA,
                  ", ".join("%d x%d" % (k, v) for k, v in sorted(oor.items()))))
        if strict:
            raise SpeciesCodeError(msg + "  strict=True.")
        if warn:
            warnings.warn(msg, OutOfRegionSpeciesWarning, stacklevel=2)
    return summary


# ==========================================================================
# StandID: dense, collision-free by construction, plus a CN lookup table
# ==========================================================================
# src-converted/dbsqlite/dbstreesin.f90:20 declares ITREE as a default INTEGER
# and line 97 binds TREE_ID with fsql3_colint, i.e. a 4-byte signed integer.
TREE_ID_MAX_SIGNED32 = 2147483647
MAX_TREES_PER_STAND = 9999          # tree_id = stand_index * 10000 + seq
_TREE_ID_STRIDE = MAX_TREES_PER_STAND + 1


def allocate_tree_id(stand_index, tree_seq):
    """Collision-free TREE_ID that provably fits in signed 32-bit.

    stand_index is 1-based (as produced by build_standid_map), tree_seq is
    1-based within the stand.  Never put a raw FIA CN here: a 14- to 15-digit
    CN overflows the INTEGER that dbstreesin.f90 reads TREE_ID into.
    """
    if stand_index < 1:
        raise ValueError("stand_index must be >= 1, got %r" % (stand_index,))
    if not (1 <= tree_seq <= MAX_TREES_PER_STAND):
        raise ValueError("tree_seq must be 1..%d, got %r"
                         % (MAX_TREES_PER_STAND, tree_seq))
    tid = stand_index * _TREE_ID_STRIDE + tree_seq
    if tid > TREE_ID_MAX_SIGNED32:
        raise OverflowError(
            "TREE_ID %d exceeds the signed 32-bit maximum %d that "
            "dbsqlite/dbstreesin.f90 reads TREE_ID into (INTEGER ITREE, "
            "line 20; fsql3_colint, line 97)." % (tid, TREE_ID_MAX_SIGNED32))
    return tid


class StandIDMap(object):
    """Dense sequential StandID scheme with a CN lookup table."""

    def __init__(self, variant, pairs, width):
        self.variant = variant
        self.width = width
        self.pairs = list(pairs)                       # [(StandID, CN), ...]
        self.by_cn = OrderedDict((cn, sid) for (sid, cn) in self.pairs)
        self.index_by_cn = OrderedDict(
            (cn, i) for i, (sid, cn) in enumerate(self.pairs, start=1))
        self.cn_by_standid = OrderedDict((sid, cn) for (sid, cn) in self.pairs)

    def __len__(self):
        return len(self.pairs)

    def standid(self, cn):
        return self.by_cn[cn]

    def index(self, cn):
        return self.index_by_cn[cn]

    def tree_id(self, cn, tree_seq):
        return allocate_tree_id(self.index_by_cn[cn], tree_seq)

    def n_distinct(self):
        return len(set(sid for (sid, _cn) in self.pairs))

    def max_collision_group(self):
        c = Counter(sid for (sid, _cn) in self.pairs)
        return max(c.values()) if c else 0

    def rows(self):
        return [{"StandID": sid, "CN": cn, "StandIndex": i}
                for i, (sid, cn) in enumerate(self.pairs, start=1)]

    def to_dataframe(self):
        try:
            import pandas as pd
        except ImportError:
            return self.rows()
        return pd.DataFrame(self.rows())


def build_standid_map(cns, variant, width=7):
    """Dense sequential StandIDs: NE_0000001 ... NE_000NNNN.

    Collision-free by construction: the identifier is the position in the
    de-duplicated CN list, not a lossy function of the CN.  Insertion order of
    the first appearance of each CN is preserved, so reruns on the same plot
    list are reproducible.  Width auto-widens if the list outgrows it.
    """
    seen = OrderedDict()
    for cn in cns:
        key = str(cn).strip()
        if key and key not in seen:
            seen[key] = None
    n = len(seen)
    need = len(str(n)) if n else 1
    width = max(width, need)
    pfx = str(variant).upper()
    pairs = [("%s_%0*d" % (pfx, width, i), cn)
             for i, cn in enumerate(seen, start=1)]
    return StandIDMap(variant=pfx, pairs=pairs, width=width)


def check_standid_collisions(cns, sid_fn):
    """Report collisions produced by an arbitrary StandID function.

    sid_fn takes one CN (as given) and returns the StandID string.  Use it to
    audit the legacy lambda
        lambda cn: f"NE_{int(cn) % 10_000_000:07d}"
    Returns dict(n_input, n_unique_cn, n_distinct_standids, max_group_size,
                 n_cns_in_collision, largest_group_standid, group_sizes).
    """
    uniq = list(OrderedDict((str(c).strip(), None) for c in cns
                            if str(c).strip()))
    sids = [sid_fn(c) for c in uniq]
    groups = Counter(sids)
    if groups:
        largest_sid, max_group = groups.most_common(1)[0]
    else:
        largest_sid, max_group = None, 0
    return {
        "n_input": len(list(cns)),
        "n_unique_cn": len(uniq),
        "n_distinct_standids": len(groups),
        "max_group_size": max_group,
        "n_cns_in_collision": sum(v for v in groups.values() if v > 1),
        "largest_group_standid": largest_sid,
        "group_sizes": dict(Counter(groups.values())),
    }


def write_standid_lookup(path, mapping):
    """Write the StandID <-> CN lookup table alongside the run outputs.

    mapping is a StandIDMap or any iterable of (StandID, CN).  Contains no
    spatial attribute of any kind; it is a bare control-number crosswalk.
    """
    if isinstance(mapping, StandIDMap):
        rows = mapping.rows()
    else:
        rows = [{"StandID": s, "CN": c, "StandIndex": i}
                for i, (s, c) in enumerate(mapping, start=1)]
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["StandID", "CN", "StandIndex"])
        w.writeheader()
        w.writerows(rows)
    return path



def legacy_standid(cn, variant="ne"):
    """The legacy, aliasing scheme, kept only so the defect stays testable."""
    return "%s_%07d" % (str(variant).upper(), int(cn) % 10_000_000)


# ==========================================================================
# Self-test
# ==========================================================================
PLOT_LIST_3321 = ("/fs/scratch/PUOM0008/crsfaaron/balfix_run_20260804/"
                  "plot_list_3321.csv")
PLOT_LIST_3321_MD5 = "5f43e67952ab71135b4b159e730dfd77"


def _load_plot_list(path=PLOT_LIST_3321):
    with open(path) as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        col = 0
        if "PLOT" in header:
            col = header.index("PLOT")
        return [r[col].strip() for r in rdr if r and r[col].strip()]


def _selftest():
    results = []

    def chk(name, got, want):
        ok = (got == want)
        results.append((ok, name, got, want))
        print("%-4s %-52s got=%-22r want=%r"
              % ("PASS" if ok else "FAIL", name, got, want))
        return ok

    print("=" * 78)
    print("fvs_input_guards self-test")
    print("  species source : %s" % SPCD_SOURCE_FILE)
    print("  species lines  : %s" % SPCD_SOURCE_LINES)
    print("=" * 78)

    # --- Defect A: species classification -------------------------------
    chk("NE species count", len(NE_SPCD), 107)
    chk("NE table slots (incl. blank 71)", len(_NE_TABLE), 108)
    chk("998 maps to OH", NE_SPCD_TO_ALPHA.get(998), "OH")
    for c in (12, 97, 316, 802, 998):
        chk("classify %d" % c, classify_spcd(c)["cls"], VALID)
    for c in (11, 20, 263, 264):
        chk("classify %d" % c, classify_spcd(c)["cls"], OUT_OF_REGION)
    for c in (9999, 0, -5, 100000, "abc", None, ""):
        chk("classify %r" % (c,), classify_spcd(c)["cls"], INVALID)

    raised = False
    try:
        validate_spcd(9999)
    except SpeciesCodeError:
        raised = True
    chk("validate_spcd(9999) raises", raised, True)

    raised = False
    try:
        validate_spcd(263)
    except SpeciesCodeError:
        raised = True
    chk("validate_spcd(263) does NOT raise by default", raised, False)

    raised = False
    try:
        validate_spcd(263, strict=True)
    except SpeciesCodeError:
        raised = True
    chk("validate_spcd(263, strict=True) raises", raised, True)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        validate_spcd(263)
        chk("validate_spcd(263) warns", len(w), 1)

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        s = validate_spcd_series([97, 97, 12, 11, 20, 263, 264])
    chk("series n_valid", s["n_valid"], 3)
    chk("series n_out_of_region", s["n_out_of_region"], 4)
    chk("series n_invalid", s["n_invalid"], 0)

    raised = False
    try:
        validate_spcd_series([97, 12, 9999])
    except SpeciesCodeError:
        raised = True
    chk("series with 9999 raises", raised, True)

    # --- Defect B: TREE_ID width ----------------------------------------
    chk("tree_id(1,1)", allocate_tree_id(1, 1), 10001)
    chk("tree_id(3321,9999) fits int32",
        allocate_tree_id(3321, 9999) <= TREE_ID_MAX_SIGNED32, True)
    chk("a 15-digit FIA CN overflows int32",
        int("107545100010661") > TREE_ID_MAX_SIGNED32, True)

    # --- Defect B: collision reproduction on the real 3,321 CN list ------
    if os.path.exists(PLOT_LIST_3321):
        cns = _load_plot_list()
        chk("plot list length", len(cns), 3321)
        legacy = check_standid_collisions(cns, lambda c: legacy_standid(c, "ne"))
        chk("legacy distinct StandIDs", legacy["n_distinct_standids"], 50)
        chk("legacy max collision group", legacy["max_group_size"], 147)
        chk("legacy CNs in a collision group", legacy["n_cns_in_collision"], 3320)
        new = build_standid_map(cns, "ne")
        chk("new distinct StandIDs", new.n_distinct(), 3321)
        chk("new max collision group", new.max_collision_group(), 1)
        chk("new first StandID", new.pairs[0][0], "NE_0000001")
        chk("new last StandID", new.pairs[-1][0], "NE_0003321")
    else:
        print("SKIP plot-list collision reproduction (%s not readable)"
              % PLOT_LIST_3321)

    n_fail = sum(1 for ok, *_ in results if not ok)
    print("-" * 78)
    print("%d checks, %d failures -- %s"
          % (len(results), n_fail, "ALL PASS" if n_fail == 0 else "FAILURES"))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(_selftest())
