#!/usr/bin/env python3
"""
5-model cohort scorecard: size-aware GOMPIT vs size-blind baseline.
Bias = (PRED - OBS)/OBS * 100, per row.  median% and aggregate% by model,
overall and by forest type (SW/MW/HW), for BA, MCuFt, BdFt.  Then the DELTA
(gompsize - baseline).  Also checks that the GOMPIT-driven model predictions
actually changed vs baseline (size term live).
"""
import numpy as np, pandas as pd

BASE = "/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_cohort/scorecard_ne_5model_cohort.csv"
SIZE = "/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_cohort_gompsize/scorecard_ne_5model_cohort.csv"

MODELS = ["fvs_base", "fvs_regional", "organon", "conus_spdep", "conus_climate"]
METRICS = [("BA", "BA_PRED", "BA_OBS"),
           ("MCuFt", "MCuFt_PRED", "MCuFt_OBS"),
           ("BdFt", "BdFt_PRED", "BdFt_OBS")]

def load(p):
    d = pd.read_csv(p)
    return d

def bias_table(d):
    """returns dict[(model, stratum, metric)] = (median%, agg%, n)"""
    out = {}
    for m in MODELS:
        dm = d[d["model"] == m]
        for stratum in ["OVERALL", "SW", "MW", "HW"]:
            ds = dm if stratum == "OVERALL" else dm[dm["FORTYPE"] == stratum]
            for name, pc, oc in METRICS:
                pred = ds[pc].to_numpy(float); obs = ds[oc].to_numpy(float)
                ok = np.isfinite(pred) & np.isfinite(obs) & (obs > 0)
                pred, obs = pred[ok], obs[ok]
                if len(pred) == 0:
                    out[(m, stratum, name)] = (np.nan, np.nan, 0); continue
                med = float(np.median((pred - obs) / obs * 100.0))
                agg = float((pred.sum() - obs.sum()) / obs.sum() * 100.0)
                out[(m, stratum, name)] = (med, agg, len(pred))
    return out

def print_block(title, tab):
    print("\n" + "=" * 96)
    print(title)
    print("=" * 96)
    hdr = f"{'model':14} {'stratum':8} " + " ".join(f"{n+'_med':>10} {n+'_agg':>10}" for n, _, _ in METRICS) + f"  {'n':>6}"
    print(hdr)
    for m in MODELS:
        for stratum in ["OVERALL", "SW", "MW", "HW"]:
            cells = []
            n = 0
            for name, _, _ in METRICS:
                med, agg, nn = tab[(m, stratum, name)]; n = nn
                cells.append(f"{med:10.2f} {agg:10.2f}")
            print(f"{m:14} {stratum:8} " + " ".join(cells) + f"  {n:6d}")

def print_delta(base, size):
    print("\n" + "=" * 96)
    print("DELTA (size-aware minus baseline; positive agg-BA delta = deficit closing)")
    print("=" * 96)
    hdr = f"{'model':14} {'stratum':8} " + " ".join(f"{n+'_dmed':>10} {n+'_dagg':>10}" for n, _, _ in METRICS)
    print(hdr)
    for m in MODELS:
        for stratum in ["OVERALL", "SW", "MW", "HW"]:
            cells = []
            for name, _, _ in METRICS:
                bmed, bagg, _ = base[(m, stratum, name)]
                smed, sagg, _ = size[(m, stratum, name)]
                cells.append(f"{smed-bmed:10.2f} {sagg-bagg:10.2f}")
            print(f"{m:14} {stratum:8} " + " ".join(cells))

def changed_check(db, ds):
    print("\n" + "=" * 96)
    print("PREDICTIONS CHANGED vs BASELINE?  (mean |dBA_PRED| per plot; 0 => identical)")
    print("=" * 96)
    key = ["PLOT", "STATE", "MEASYEAR1", "MEASYEAR2", "model"]
    mb = db[key + ["BA_PRED", "MCuFt_PRED", "BdFt_PRED"]].merge(
        ds[key + ["BA_PRED", "MCuFt_PRED", "BdFt_PRED"]], on=key, suffixes=("_b", "_s"))
    for m in MODELS:
        g = mb[mb["model"] == m]
        dba = (g["BA_PRED_s"] - g["BA_PRED_b"]).abs()
        dmc = (g["MCuFt_PRED_s"] - g["MCuFt_PRED_b"]).abs()
        dbf = (g["BdFt_PRED_s"] - g["BdFt_PRED_b"]).abs()
        nchg = int((dba > 1e-6).sum())
        print(f"  {m:14} n={len(g):5d}  meanABS dBA={dba.mean():9.4f}  dMCuFt={dmc.mean():9.3f}  "
              f"dBdFt={dbf.mean():9.2f}  rows_changed={nchg}")

def main():
    db = load(BASE); ds = load(SIZE)
    print(f"baseline rows={len(db)}  gompsize rows={len(ds)}")
    tb = bias_table(db); tsz = bias_table(ds)
    print_block("BASELINE (size-blind GOMPIT)  scorecard_conus_cohort", tb)
    print_block("SIZE-AWARE GOMPIT  scorecard_conus_cohort_gompsize", tsz)
    print_delta(tb, tsz)
    changed_check(db, ds)

if __name__ == "__main__":
    main()
