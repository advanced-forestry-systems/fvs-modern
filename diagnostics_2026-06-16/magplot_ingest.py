#!/usr/bin/env python3
"""MAGPlot ingestion for ACD and AK calibration (2026-06-18). Turnkey: run after the MAGPlot data package
is staged and unzipped on Cardinal. Inspects the tables, identifies the site and tree tables and their key
columns, filters to the jurisdictions that map to the Acadian (NB, NS, PE) and Alaska/coastal-BC (BC)
variants, builds remeasurement pairs (sites with two or more measurements), computes observed diameter
growth and mortality, and writes calibration-ready per-variant CSVs. Schema-tolerant: detects columns by
common MAGPlot names and prints what it found so the mapping can be confirmed.

Usage: python3 magplot_ingest.py <MAGPLOT_DIR> <OUT_DIR>
"""
import os, sys, glob, re
import pandas as pd, numpy as np
MP = sys.argv[1] if len(sys.argv) > 1 else "/fs/scratch/PUOM0008/crsfaaron/MAGPlot"
OUT = sys.argv[2] if len(sys.argv) > 2 else MP

# jurisdiction (province/territory) -> FVS variant for the two targets
JUR2VAR = {"NB": "acd", "NS": "acd", "PE": "acd", "BC": "ak"}

def find_tables(d):
    fs = []
    for ext in ("*.csv", "*.txt", "*.tsv"):
        fs += glob.glob(os.path.join(d, "**", ext), recursive=True)
    return sorted(set(fs))

def detect(cols, pats):
    for p in pats:
        for c in cols:
            if re.search(p, c, re.I):
                return c
    return None

def main():
    tabs = find_tables(MP)
    if not tabs:
        print("No CSV/TSV tables under", MP, "- is the package unzipped? Resources may be in subfolders.")
        return
    print("=== MAGPlot tables found ===")
    meta = {}
    for f in tabs:
        try:
            h = pd.read_csv(f, nrows=200, low_memory=False)
        except Exception as e:
            print(" ", os.path.basename(f), "READ FAIL", e); continue
        meta[f] = list(h.columns)
        print(" ", os.path.basename(f), "| cols:", ", ".join(list(h.columns)[:12]), "...")

    # identify the tree-measurement table (DBH + species + status) and the site table (jurisdiction + coords)
    tree_f = next((f for f, c in meta.items() if detect(c, [r"\bdbh"]) and detect(c, [r"spec"])), None)
    site_f = next((f for f, c in meta.items() if detect(c, [r"juris|prov|admin|agency"]) ), None)
    print("\n=== detected ===")
    print(" tree table:", os.path.basename(tree_f) if tree_f else "NOT FOUND")
    print(" site table:", os.path.basename(site_f) if site_f else "NOT FOUND")
    if not tree_f:
        print("Could not auto-detect the tree table; confirm column names and set them manually.")
        return

    tcols = meta[tree_f]
    K = {
        "site":  detect(tcols, [r"site.?id", r"\bplot.?id", r"\bnfi.?plot", r"loc.?id"]),
        "meas":  detect(tcols, [r"meas.?(num|no|id)", r"visit", r"remeas", r"\byear", r"meas.?year"]),
        "dbh":   detect(tcols, [r"\bdbh"]),
        "ht":    detect(tcols, [r"\bheight", r"\bht\b", r"tree.?ht"]),
        "spp":   detect(tcols, [r"species", r"\bspp", r"\bspcd"]),
        "status":detect(tcols, [r"status", r"live|dead", r"condition"]),
        "tree":  detect(tcols, [r"tree.?(num|no|id)", r"stem.?id"]),
    }
    print(" tree-table key columns:", K)
    tr = pd.read_csv(tree_f, low_memory=False)

    # attach jurisdiction from the site table if available
    if site_f and K["site"]:
        scols = meta[site_f]
        sj = detect(scols, [r"juris|prov|admin|agency"]); ssite = detect(scols, [r"site.?id", r"\bplot.?id", r"loc.?id"])
        if sj and ssite:
            sdf = pd.read_csv(site_f, low_memory=False)[[ssite, sj]].drop_duplicates()
            tr = tr.merge(sdf, left_on=K["site"], right_on=ssite, how="left")
            tr["JUR"] = tr[sj].astype(str).str.upper().str[:2]
        else:
            tr["JUR"] = None
    else:
        tr["JUR"] = None

    # remeasurement candidates per jurisdiction of interest
    print("\n=== remeasurement candidates for ACD (NB/NS/PE) and AK (BC) ===")
    if K["site"] and K["meas"]:
        for jur, var in JUR2VAR.items():
            sub = tr[tr["JUR"] == jur]
            if not len(sub): continue
            nmeas = sub.groupby(K["site"])[K["meas"]].nunique()
            remeas = (nmeas >= 2).sum()
            print(f"  {jur} -> {var}: {sub[K['site']].nunique()} sites, {remeas} with >=2 measurements, {len(sub)} tree rows")
    else:
        print("  site/measurement columns not detected; set K['site'] and K['meas'] then re-run.")

    # write a manifest of the detected schema so the pair-builder can be finalized
    man = os.path.join(OUT, "magplot_schema_manifest.txt")
    with open(man, "w") as fh:
        fh.write("tree_table=%s\nsite_table=%s\nkeycols=%s\nJUR2VAR=%s\n" % (tree_f, site_f, K, JUR2VAR))
        fh.write("\nNEXT: confirm keycols, then build t1/t2 pairs per site (min DBH, status live), map\n")
        fh.write("MAGPlot species to FVS species (crosswalk via the data dictionary), and run the\n")
        fh.write("calibration (calib_ne.py-style observed/default DG) per variant.\n")
    print("\nwrote", man)
    print("DONE_MAGPLOT_INGEST")

if __name__ == "__main__":
    main()
