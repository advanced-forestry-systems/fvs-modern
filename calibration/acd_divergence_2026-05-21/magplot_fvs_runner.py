#!/usr/bin/env python3
"""
magplot_fvs_runner.py
=====================
Ingest Canadian MAGPlot (New Brunswick) tree lists into the Fortran FVS engine
(FVS-NE / FVS-ACD, default + calibrated) by reusing the proven DATABASE/STANDSQL
path in silc_fvs_runner_v2.run_one (which avoids the Sequential-READ-after-EOF
inventory error).

Root cause of the prior block: fvs2py needs Python >= 3.11 (enum.StrEnum). The
cluster default is 3.9, so the import failed before any inventory could load.
Run this under:  module load python/3.12

Smoke-test mode (default): convert a handful of NB plots and run FVS-NE default
to confirm the engine ingests MAGPlot trees. Scale up by raising --nstands and
adding variants/configs.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd

# fvs2py uses enum.StrEnum (3.11+) and typing.ParamSpec (3.10+). The cluster
# default python3 is 3.9, which fails to import fvs2py. Fail fast with guidance.
if sys.version_info < (3, 11):
    raise SystemExit(
        "fvs2py requires Python >= 3.11 (enum.StrEnum, typing.ParamSpec). "
        "On Cardinal run:  module load python/3.12")

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", "/users/PUOM0008/crsfaaron/fvs-modern")
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "fvs2py"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "microfvs"))
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-conus/python")

# reuse the proven inventory build + run path
from silc_fvs_runner_v2 import run_one, build_standinit_df  # noqa: E402

ACRES_PER_HA = 2.4710538147
INCH_PER_CM = 1 / 2.54
FT_PER_M = 1 / 0.3048

# MAGPlot genus.species (species_gs) -> FIA SPCD. Generics/unknowns -> other.
GS_TO_SPCD = {
    "ABIE.BAL": 12, "PICE.MAR": 95, "ACER.RUB": 316, "BETU.PAP": 375,
    "PICE.RUB": 97, "POPU.TRE": 746, "THUJ.OCC": 241, "BETU.ALL": 371,
    "ACER.SAH": 318, "ACER.SAC": 318, "PICE.GLA": 94, "BETU.POP": 379,
    "FAGU.GRA": 531, "PINU.BAN": 105, "PRUN.PEN": 761, "PINU.STR": 129,
    "ACER.PEN": 315, "ACER.SPI": 319, "LARI.LAR": 71, "FRAX.AME": 541,
    "POPU.GRA": 743, "POPU.BAL": 741, "FRAX.NIG": 543, "SALI.SPP": 920,
    "TSUG.CAN": 261, "PINU.RES": 125, "PRUN.VIR": 763, "SORB.AME": 935,
    "OSTR.VIR": 701, "PICE.ABI": 91, "QUER.SPP": 833, "PRUN.SER": 762,
    "PINU.SYL": 130, "ULMU.SPP": 972, "CRAT.SPP": 500, "TILI.AME": 951,
    "FRAX.PEN": 544, "JUGL.CIN": 601,
    # generics / shrubs / unknown -> FVS "other" buckets
    "ALNU.SPP": 998, "ILEX.MUC": 998, "VIBU.CAS": 998, "AMEL.SPP": 998,
    "MALU.SPP": 998, "SAMB.RAC": 998, "SAMB.NIG": 998, "CORN.ALT": 998,
    "VIBU.LAN": 998, "UNKN.SPP": 998, "GENH.SPP": 998, "GENC.SPP": 298,
}

def build_treeinit_magplot(plot_trees: pd.DataFrame, stand_id: str) -> pd.DataFrame:
    rows = []
    tid = 1
    for _, r in plot_trees.iterrows():
        spcd = GS_TO_SPCD.get(str(r["species_gs"]).strip().upper())
        if spcd is None:
            continue
        try:
            dbh_cm = float(r["dbh"])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(dbh_cm) or dbh_cm <= 0:
            continue
        try:
            stem_ha = float(r["stem_ha"])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(stem_ha) or stem_ha <= 0:
            continue
        ht_ft = 0.0
        try:
            h = float(r["height"])
            if np.isfinite(h) and h > 0:
                ht_ft = h * FT_PER_M
        except (TypeError, ValueError):
            pass
        rows.append({
            "stand_id":   stand_id,
            "plot_id":    1,
            "tree_id":    tid,
            "tree_count": round(stem_ha / ACRES_PER_HA, 5),  # TPH -> TPA
            "species":    spcd,
            "diameter":   round(dbh_cm * INCH_PER_CM, 3),
            "ht":         round(ht_ft, 1),
            "crratio":    40,
        })
        tid += 1
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--magdir", default="/users/PUOM0008/crsfaaron/magplot")
    ap.add_argument("--nstands", type=int, default=6)
    ap.add_argument("--min-trees", type=int, default=10)
    ap.add_argument("--variant", default="ne")
    ap.add_argument("--config", default="default")
    ap.add_argument("--num-cycles", type=int, default=2)
    args = ap.parse_args()

    trees = pd.read_csv(os.path.join(args.magdir, "magp_trees_nb.csv"),
                        usecols=["magp_site_id","plot_id","meas_num","species_gs",
                                 "dbh","height","stem_ha","tree_status"],
                        low_memory=False)
    hdr = pd.read_csv(os.path.join(args.magdir, "magp_tree_header_nb.csv"))

    # initial measurement, live trees only
    t0 = trees[(trees["meas_num"] == 0) & (trees["tree_status"] == "L")].copy()
    t0 = t0[t0["dbh"].notna() & (t0["dbh"] > 0)]
    cov = t0["species_gs"].astype(str).str.upper().isin(GS_TO_SPCD).mean()
    print(f"[magplot] NB initial live trees: {len(t0):,}  species crosswalk coverage: {cov*100:.1f}%")

    # year per plot from header (median meas_year for meas_num==0)
    hdr0 = hdr[hdr["meas_num"] == 0]
    yr_lookup = dict(zip(
        hdr0["magp_site_id"].astype(str) + "_" + hdr0["plot_id"].astype(str),
        hdr0["meas_year"]))

    grp = t0.groupby(["magp_site_id", "plot_id"])
    picked, results = 0, []
    for (site, plot), g in grp:
        if len(g) < args.min_trees:
            continue
        stand_id = f"{site}"[:25]
        meta = {"elev_m": 200.0, "csi_m": 12.0,
                "inv_year": int(yr_lookup.get(f"{site}_{plot}", 2021))}
        stand_df = build_standinit_df(stand_id, meta, variant=args.variant)
        tree_df = build_treeinit_magplot(g, stand_id)
        if tree_df.empty:
            continue
        res = run_one(stand_df, tree_df, stand_id, args.variant, args.config,
                      num_cycles=args.num_cycles, cycle_length=5)
        if "error" in res:
            print(f"  {stand_id}: ERROR {res['error']}")
            continue
        s = res["summary"]
        if s is None or len(s) == 0:
            print(f"  {stand_id}: ran but EMPTY summary ({len(tree_df)} trees in)")
            continue
        ba_col = next((c for c in s.columns if c.lower() in ("ba","baa","tcuft_ba","ba_ft2")), None)
        ba0 = s.iloc[0][ba_col] if ba_col else float("nan")
        ban = s.iloc[-1][ba_col] if ba_col else float("nan")
        print(f"  {stand_id}: OK  {len(tree_df)} trees in  rows={len(s)}  "
              f"BA[{ba_col}] {ba0} -> {ban}")
        results.append((stand_id, len(tree_df), ba0, ban))
        picked += 1
        if picked >= args.nstands:
            break

    print(f"\n[magplot] SMOKE TEST: {picked} stands ingested + projected through FVS-{args.variant.upper()} ({args.config})")
    if picked > 0:
        print("[magplot] RESULT: fvs2py ingests MAGPlot tree lists; the inventory blocker is the Python<3.11 StrEnum import, resolved under python/3.12.")
    else:
        print("[magplot] RESULT: no stands ran; inspect errors above.")

if __name__ == "__main__":
    main()
