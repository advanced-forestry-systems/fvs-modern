"""PROPOSED, NOT APPLIED (2026-10-29). Patch conus_catalog_runs.py so live FIA tree records with a null
TPA_UNADJ are dropped before the FVS_TreeInit rows are built.

Diagnosis (catalog/diag_y0/diag_y0_summary.csv, diag_y0_trees.csv, diag_nulltpa_stands.csv): the engine
reproduces its input exactly at year 0 (tree count, sum of TPA and basal area match the FVS_TreeInit rows
for the four test stands). The excess over BALIVE is already in the FVS_TreeInit rows: FIA TREE records with
STATUSCD = 1 and TPA_UNADJ = NULL (P2A_GRM_FLG = Y, legacy subplot numbers 8, 1xx, 2xx, on 2002 to 2008
SN, CS, LS and NE plots) contribute nothing to BALIVE, but fia_stand_generator.build_tree_init fills a null
tree_count with 6.0 TPA, so each such record enters the engine as a 6 TPA tree. 250 of the 256 eastern
default/native seed-gate failures in catalog_seedgate_v4.csv carry at least one such record, and dropping
those records puts all 250 inside the 20 percent gate (fail_after_drop = 0 for ACD, CS, LS, NE, SN). The
other 6 failures have BALIVE below 3 ft2/ac, where the engine's integer summary BA cannot meet a 20 percent
relative tolerance: a gate artefact, not over-seeding.

Usage, only on a fresh re-run of the catalog (never on the archived catalog files):
    python3 patch_overseed_proposal.py conus_catalog_runs.py
"""
import sys
p = sys.argv[1]; t = open(p).read()
old = '''        t = tree[(tree["PLT_CN"] == cn) & (tree["CONDID"] == condid)]
        if t.empty:
'''
new = '''        t = tree[(tree["PLT_CN"] == cn) & (tree["CONDID"] == condid)]
        # Live records with a null TPA_UNADJ (P2A_GRM_FLG = Y: trees carried from the prior periodic
        # inventory for growth, removal and mortality only, on legacy subplot numbers) contribute nothing
        # to BALIVE, but build_tree_init fills a null tree_count with 6.0 TPA. Drop them (2026-10-29:
        # 250 of 256 eastern seed-gate failures were this).
        n_null_tpa = int(t["TPA_UNADJ"].isna().sum()); t = t[t["TPA_UNADJ"].notna()]
        if t.empty:
'''
assert t.count(old) == 1, "anchor not found"
t = t.replace(old, new)
old2 = "n_relabel=n_relabel, n_drop=n_drop,"
assert t.count(old2) == 1, "gate anchor not found"
t = t.replace(old2, "n_relabel=n_relabel, n_drop=n_drop, n_null_tpa=n_null_tpa,")
open(p, "w").write(t); print("patched", p)
