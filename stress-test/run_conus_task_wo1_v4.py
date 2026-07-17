#!/usr/bin/env python3
"""CONUS-wide FVS projection task (v4): full re-run with the v3 expansion fix.

DIAGNOSIS (2026-06-21, supersedes the v3 "join-broken" framing):
  The "42/5000 matched" reported for CR batch b0 is NOT a join bug and NOT a
  stand-universe mismatch between the engine and the fvs-conus equation arm.
  It is a BATCH-ORDERING ARTIFACT. standinit_CR.csv batch b0 = rows [0:5000] are
  ALL California (STATE=6). California is the one CR state where the standinit
  STAND_CN/STAND_ID universe (INV_YEAR 2021 FIA pull) does NOT match the treeinit
  STAND_CN universe (older FIA pull): only ~1.1% (86/8170) of CA standinit stands
  exist in CA_FVS_TREEINIT_PLOT.csv. So b0 in isolation matches ~42.

  Aggregated over ALL 64 CR batches, engine v2 matched 52,423 stands-with-trees
  and PROJECTED 51,256 -- byte-identical to the fvs-conus eq arm's 51,256. The
  two arms were ALWAYS on the same stand universe. Both arms inner-join
  standinit<->treeinit on STAND_CN and drop the same non-overlapping stands.
  Overall standinit_CR overlap with treeinit = 16.6% (52,423/315,529), and that
  16.6% is exactly what BOTH arms project. The non-overlap is a genuine FIA-pull
  vintage difference (different plots/cycles), NOT fixable by CN normalization.

  v3 was an INCOMPLETE run: its SLURM array was cancelled (SIGTERM) after only
  task 0 (SN b0) plus a CR b0 test, so out_conus_engine_v3 has just two batch
  files. The v3 expansion fix itself is CORRECT: on the 42 matched CR b0 stands,
  engine v3 year-0 BA == eq year-0 BA to within 0.024% mean (0.39% max).

v4 = v3 driver (per-acre expansion fix retained: inv_plot_size=1, num_plots=1,
brk_dbh=999, basal_area_factor=0) run over the FULL 381-task manifest, writing
to out_conus_engine_v4. Ledgers now also record per-STATE match counts so a
single-batch ledger is interpretable. No data files modified; v2/v3 untouched.
"""
from __future__ import annotations
import argparse, json, logging, os, sys, time
import numpy as np, pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-conus/python")
import perseus_100yr_projection as P  # noqa: E402

# ---- SEEDING FIX: wrap build_fvs_standinit to force per-acre plot design ----
_orig_build_standinit = P.build_fvs_standinit
def build_fvs_standinit_peracre(plot_data, stand_id, variant):
    sdf = _orig_build_standinit(plot_data, stand_id, variant)
    # supplied TREE_COUNT is already per-acre TPA (TPA_UNADJ); take it verbatim.
    sdf["inv_plot_size"] = 1.0
    sdf["num_plots"] = 1
    sdf["brk_dbh"] = 999.0
    sdf["basal_area_factor"] = 0.0
    return sdf
P.build_fvs_standinit = build_fvs_standinit_peracre  # used by all calls below

ARMS = (((None, "gompit"),) if os.environ.get("GOMPIT_ARM")
        else ((None, "default"), ("calibrated", "calibrated")))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("conus")

FIPS = {1:"AL",2:"AK",4:"AZ",5:"AR",6:"CA",8:"CO",9:"CT",10:"DE",12:"FL",13:"GA",
        16:"ID",17:"IL",18:"IN",19:"IA",20:"KS",21:"KY",22:"LA",23:"ME",24:"MD",
        25:"MA",26:"MI",27:"MN",28:"MS",29:"MO",30:"MT",31:"NE",32:"NV",33:"NH",
        34:"NJ",35:"NM",36:"NY",37:"NC",38:"ND",39:"OH",40:"OK",41:"OR",42:"PA",
        44:"RI",45:"SC",46:"SD",47:"TN",48:"TX",49:"UT",50:"VT",51:"VA",53:"WA",
        54:"WV",55:"WI",56:"WY"}

TPA_PER_HA = 2.4710538
BA_CONST = np.pi / (4.0 * 144.0)

def stand_metrics(tl: pd.DataFrame) -> dict:
    if tl is None or tl.empty:
        return {"BA_FT2AC": 0.0, "QMD_IN": 0.0, "TPH": 0.0}
    sum_ba_in2 = 0.0; sum_tpa = 0.0
    for _, tree in tl.iterrows():
        dbh = tree.get("DBH", tree.get("Dbh", 0))
        tpa = tree.get("TPA", tree.get("Tpa", 1.0))
        if pd.isna(dbh) or dbh < 1.0 or pd.isna(tpa):
            continue
        d2 = float(dbh) * float(dbh); t = float(tpa)
        sum_ba_in2 += d2 * t; sum_tpa += t
    if sum_tpa <= 0.0:
        return {"BA_FT2AC": 0.0, "QMD_IN": 0.0, "TPH": 0.0}
    ba = BA_CONST * sum_ba_in2
    qmd = np.sqrt(sum_ba_in2 / sum_tpa)
    tph = sum_tpa * TPA_PER_HA
    return {"BA_FT2AC": round(float(ba),4), "QMD_IN": round(float(qmd),4), "TPH": round(float(tph),4)}

def treeinit_for_stand(fvs_rows: pd.DataFrame, stand_id: str) -> pd.DataFrame:
    recs = []
    for i, t in enumerate(fvs_rows.itertuples(index=False)):
        d = getattr(t, "DIAMETER", np.nan)
        if pd.isna(d) or float(d) < 1.0:
            continue
        ht = getattr(t, "HT", 0); cr = getattr(t, "CRRATIO", 0)
        recs.append({
            "stand_id": stand_id,
            "plot_id": int(float(getattr(t, "PLOT_ID", 1) or 1)),
            "tree_id": i + 1,
            "tree_count": float(getattr(t, "TREE_COUNT", 1.0) or 1.0),
            "species": int(float(getattr(t, "SPECIES", 0) or 0)),
            "diameter": round(float(d), 1),
            "ht": round(float(ht), 0) if pd.notna(ht) and ht not in ("", None) and float(ht) > 0 else 0,
            "crratio": int(float(cr)) if pd.notna(cr) and cr not in ("", None) and float(cr) > 0 else 0,
        })
    return pd.DataFrame(recs)

def lookup_task(manifest, task_id):
    with open(manifest) as fh:
        for line in fh:
            idx, variant, batch_id, batch_size = line.rstrip("\n").split("\t")
            if int(idx) == task_id:
                return variant, int(batch_id), int(batch_size)
    raise SystemExit(f"task {task_id} not in {manifest}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--standinit-dir", required=True)
    ap.add_argument("--treeinit-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--task-id", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", "-1")))
    a = ap.parse_args()
    os.makedirs(a.output_dir, exist_ok=True)
    variant, batch_id, batch_size = lookup_task(a.manifest, a.task_id)
    t0 = time.time()
    done_csv = os.path.join(a.output_dir, f"conus_{variant.lower()}_b{batch_id}.csv")
    done_ledger = os.path.join(a.output_dir, f"ledger_{variant.lower()}_b{batch_id}.json")
    if os.path.exists(done_csv) and os.path.exists(done_ledger) and os.path.getsize(done_csv) > 0:
        log.info(f"task {a.task_id}: {variant} b{batch_id} already done; skipping"); return
    log.info(f"task {a.task_id}: {variant} batch {batch_id} [SEEDING FIX v3]")

    si_all = pd.read_csv(os.path.join(a.standinit_dir, f"standinit_{variant}.csv"), low_memory=False)
    si = si_all.iloc[batch_id*batch_size:(batch_id+1)*batch_size].reset_index(drop=True)
    si["STAND_CN"] = si["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")

    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    rows, failures = [], []
    n_with_trees = n_proj = n_ba0 = 0
    state_comp = {}
    cache = {}
    for state_fips, grp in si.groupby("STATE"):
        try:
            state = FIPS[int(float(state_fips))]
        except (KeyError, ValueError, TypeError):
            failures.append({"stand_cn": None, "stage": "state_map", "detail": str(state_fips)}); continue
        tfile = os.path.join(a.treeinit_dir, f"{state}_FVS_TREEINIT_PLOT.csv")
        if not os.path.exists(tfile):
            failures.append({"stand_cn": None, "stage": "no_treeinit", "detail": state, "n": int(len(grp))})
            state_comp.setdefault(state, [0,0]); state_comp[state][0] += int(len(grp)); continue
        if state not in cache:
            tt = pd.read_csv(tfile, low_memory=False)
            tt["STAND_CN"] = tt["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
            cache[state] = {k: v for k, v in tt.groupby("STAND_CN")}
        by_cn = cache[state]
        state_comp.setdefault(state, [0, 0])
        state_comp[state][0] += int(len(grp))
        for _, stand in grp.iterrows():
            cn = stand["STAND_CN"]
            fvs_rows = by_cn.get(cn)
            if fvs_rows is None or fvs_rows.empty:
                continue
            n_with_trees += 1
            state_comp[state][1] += 1
            sid = f"S{cn}"
            inv_year = int(float(stand.get("INV_YEAR") or 2010))
            plot_data = {"INVYR": inv_year, "LAT": stand.get("LATITUDE"), "LON": stand.get("LONGITUDE"),
                         "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
                         "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50}
            try:
                sdf = P.build_fvs_standinit(plot_data, sid, variant.lower())  # patched -> per-acre
                tdf = treeinit_for_stand(fvs_rows, sid)
            except Exception as e:
                failures.append({"stand_cn": cn, "stage": "build", "detail": str(e)}); continue
            if tdf.empty:
                continue
            n_proj += 1
            for cfg, label in ARMS:
                try:
                    fr = P.run_fvs_projection(sdf, tdf, sid, variant.lower(),
                                              config_version=cfg, num_cycles=20, cycle_length=5)
                    for cy, tl in sorted(fr["treelists"].items()):
                        py = cy - inv_year
                        if py < 0: continue
                        agb = P.compute_plot_agb(tl, nsbe)
                        sm = stand_metrics(tl)
                        if py == 0 and label == "default" and sm["BA_FT2AC"] == 0.0:
                            n_ba0 += 1
                            failures.append({"stand_cn": cn, "stage": "ba0_year0", "detail": "BA=0 after fix"})
                        rows.append({"STAND_CN": cn, "STATE": state, "YEAR": cy, "PROJ_YEAR": py,
                                     "VARIANT": variant.upper(), "CONFIG": label,
                                     "AGB_TONS_AC": round(float(agb), 4),
                                     "BA_FT2AC": sm["BA_FT2AC"], "QMD_IN": sm["QMD_IN"], "TPH": sm["TPH"]})
                except Exception as e:
                    failures.append({"stand_cn": cn, "config": label, "stage": "project", "detail": str(e)})

    tag = f"{variant.lower()}_b{batch_id}"
    pd.DataFrame(rows).to_csv(os.path.join(a.output_dir, f"conus_{tag}.csv"), index=False)
    ledger = {"task_id": a.task_id, "variant": variant.upper(), "batch_id": batch_id,
              "seeding_fix": "inv_plot_size=1,num_plots=1,brk_dbh=999 (per-acre passthrough)",
              "n_stands_in_batch": int(len(si)), "n_stands_with_trees": n_with_trees,
              "n_stands_projected": n_proj, "n_year0_ba0_default": n_ba0,
              "n_output_rows": len(rows), "n_failures": len(failures),
              "state_composition": {k:{"n_in_batch":v[0],"n_matched":v[1]} for k,v in sorted(state_comp.items())},
              "elapsed_sec": round(time.time()-t0, 1), "failures": failures[:500]}
    with open(os.path.join(a.output_dir, f"ledger_{tag}.json"), "w") as f:
        json.dump(ledger, f, indent=2)
    log.info(f"task {a.task_id} done: {n_with_trees} matched, {n_proj} projected, "
             f"{len(rows)} rows, {n_ba0} year0-BA0, {len(failures)} failures, {ledger['elapsed_sec']}s")

if __name__ == "__main__":
    main()
