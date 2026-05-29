#!/usr/bin/env python3
"""Pull CFI plot tree lists at YEAR=2000 for backfill plots and write
in AGM-input format. Targets Cedar A+B (PLOT 1106) and Mixedwood A+B
(PLOTS 1101, 1102, 1109)."""
import pandas as pd
od = "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"
out_dir = f"{od}/cfi_backfill_input"
import os; os.makedirs(out_dir, exist_ok=True)

tr = pd.read_csv(f"{od}/TREE.csv")
sm = pd.read_csv(f"{od}/silc_cfi_plot_strata_map.csv")
EXPF = 5.0  # SILC CFI 1/5-ac fixed plots

# Backfill targets
backfill = [
    ("Cedar", "A+B (high)", [1106]),
    ("Mixedwood", "A+B (high)", [1101, 1102, 1109]),
]

manifest_rows = []
for ft, dc, plots in backfill:
    for p in plots:
        t = tr[(tr.PLOT==p) & (tr.MEASYEAR==2000) & (tr.STATUSCD==1) &
               tr.DIA_IN.notna() & (tr.DIA_IN >= 4.5)]
        if len(t) == 0:
            print(f"  {ft} / {dc} PLOT {p}: 0 trees!")
            continue
        out = pd.DataFrame({
            "STAND": f"CFI_{int(p):04d}",
            "YEAR": 2000,
            "PLOT": p,
            "TREE": t["TREE"].values,
            "SP": t["COMMON_NAME"].values,
            "SPCD": t["SPCD"].values,
            "DBH": t["DIA_IN"].values,
            "HT": t["HT_FT"].values,
            "HCB": None, "EXPF": EXPF, "Form": 1, "Risk": 1,
        })
        fn = f"backfill_{int(p):04d}_2000_tree.csv"
        out.to_csv(f"{out_dir}/{fn}", index=False)
        manifest_rows.append({
            "forest_type": ft, "density_class": dc, "PLOT": int(p),
            "tree_list_file": f"cfi_backfill_input/{fn}", "n_trees": len(out),
        })
        print(f"  {ft} / {dc} PLOT {p}: {len(out)} trees -> {fn}")

manifest = pd.DataFrame(manifest_rows)
manifest.to_csv(f"{od}/cfi_backfill_manifest.csv", index=False)
print(f"\nWrote {len(manifest)} pair tree files to {out_dir}/")
