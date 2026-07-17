#!/usr/bin/env python3
"""Arm 3 (fvs-conus species-DEPENDENT, parametric + residual uncertainty).

CONUS adaptation of perseus_uncertainty_projection.py onto the proven WO-1
StandInit + FVS TreeInit subprocess path (run_conus_task_wo1.py). For a
subsample of K matched stands per variant, runs the calibrated point projection
plus N posterior draws (UncertaintyEngine). The across-draw spread of per-ha
carbon at each projection year is the PARAMETRIC band; a residual-variance term
is added downstream so intervals are predictive.

WO-1 decision (A. Weiskittel, 2026-06-20): the per-species SDIMAX block in the
draw keywords carries the same field-order bug WO-1 disabled in config_loader,
so it is STRIPPED here. Arm 3 varies the 5 non-density components (diameter
growth, height-diameter, height increment, mortality, crown ratio) and leaves
SDImax at the FVS default, matching the signed-off WO-1 base. Density spread
folds into the residual term.

Output: out_conus_arm3/arm3_<variant>.csv with columns
  STAND_CN, STATE, YEAR, PROJ_YEAR, VARIANT, KIND (point|draw), DRAW, AGB_TONS_AC
plus ledger_<variant>.json.

Usage (SLURM array task, one variant per task):
  python run_conus_arm3.py --manifest arm3_manifest.tsv \
     --standinit-dir standinit_by_variant --treeinit-dir .../treeinit_h \
     --output-dir .../out_conus_arm3 --task-id $SLURM_ARRAY_TASK_ID \
     --k-stands 120 --n-draws 30 --seed 42
"""
from __future__ import annotations
import argparse, json, logging, os, re, sys, time
import numpy as np, pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-conus/python")
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs_stress")
import perseus_100yr_projection as P  # noqa: E402
from perseus_uncertainty_projection import run_fvs_with_draw  # noqa: E402
from run_conus_task_wo1 import treeinit_for_stand, FIPS  # noqa: E402
from config.uncertainty import UncertaintyEngine  # noqa: E402
from config.config_loader import FvsConfigLoader  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("arm3")


def strip_sdimax(block: str) -> str:
    """Remove the SDIMAX keyword lines and their comment header from a draw
    keyword block (WO-1 decision: density draws disabled, FVS default SDImax)."""
    out = []
    for ln in block.splitlines():
        s = ln.strip()
        if s.startswith("SDIMAX"):
            continue
        if s.startswith("!!") and "SDI max" in s:
            continue
        out.append(ln)
    return "\n".join(out)


def lookup_task(manifest, task_id):
    with open(manifest) as fh:
        for line in fh:
            idx, variant = line.rstrip("\n").split("\t")[:2]
            if int(idx) == task_id:
                return variant
    raise SystemExit(f"task {task_id} not in {manifest}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--standinit-dir", required=True)
    ap.add_argument("--treeinit-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--task-id", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", "-1")))
    ap.add_argument("--k-stands", type=int, default=120)
    ap.add_argument("--n-draws", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    os.makedirs(a.output_dir, exist_ok=True)
    variant = lookup_task(a.manifest, a.task_id)
    vl, vu = variant.lower(), variant.upper()
    t0 = time.time()

    done_csv = os.path.join(a.output_dir, f"arm3_{vl}.csv")
    done_ledger = os.path.join(a.output_dir, f"ledger_{vl}.json")
    if os.path.exists(done_csv) and os.path.exists(done_ledger) and os.path.getsize(done_csv) > 0:
        log.info(f"{vu} already done; skipping"); return
    log.info(f"task {a.task_id}: variant {vu}, K={a.k_stands}, N={a.n_draws}")

    cfgdir = os.path.join(PROJECT_ROOT, "config")
    try:
        eng = UncertaintyEngine(vl, config_dir=cfgdir, seed=a.seed)
        if not eng.draws_available:
            raise RuntimeError("no draws")
        dft = FvsConfigLoader(vl, version="default", config_dir=cfgdir).config
    except Exception as e:
        log.error(f"{vu}: uncertainty engine init failed: {e}")
        with open(done_ledger, "w") as f:
            json.dump({"variant": vu, "error": str(e)}, f)
        pd.DataFrame().to_csv(done_csv, index=False)
        return

    si_all = pd.read_csv(os.path.join(a.standinit_dir, f"standinit_{variant}.csv"), low_memory=False)
    si_all["STAND_CN"] = si_all["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
    nsbe = P.NSBECalculator(P.NSBE_ROOT)

    # cache treeinit per state; collect matched stands until K reached
    cache = {}
    matched = []  # (cn, state, inv_year, stand_row, fvs_rows)
    rng = np.random.default_rng(a.seed)
    # shuffle stand order for a representative subsample across the variant
    order = rng.permutation(len(si_all))
    for ix in order:
        if len(matched) >= a.k_stands:
            break
        stand = si_all.iloc[ix]
        try:
            state = FIPS[int(float(stand["STATE"]))]
        except (KeyError, ValueError, TypeError):
            continue
        if state not in cache:
            tfile = os.path.join(a.treeinit_dir, f"{state}_FVS_TREEINIT_PLOT.csv")
            if not os.path.exists(tfile):
                cache[state] = None; continue
            tt = pd.read_csv(tfile, low_memory=False)
            tt["STAND_CN"] = tt["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
            cache[state] = {k: v for k, v in tt.groupby("STAND_CN")}
        by_cn = cache[state]
        if by_cn is None:
            continue
        fvs_rows = by_cn.get(stand["STAND_CN"])
        if fvs_rows is None or fvs_rows.empty:
            continue
        inv_year = int(float(stand.get("INV_YEAR") or 2010))
        matched.append((stand["STAND_CN"], state, inv_year, stand, fvs_rows))

    log.info(f"{vu}: {len(matched)} matched stands selected")
    # COMMON posterior draw indices applied to EVERY stand, so each draw is a
    # coherent parameter set for the whole population. The spread across draws of
    # the population per-ha mean is the parametric band (NOT mixed with stand
    # variation). Sampled once per variant.
    common_draws = rng.choice(eng.n_draws, size=min(a.n_draws, eng.n_draws), replace=False)
    rows, failures = [], []
    n_point = n_draw = 0
    for (cn, state, inv_year, stand, fvs_rows) in matched:
        sid = f"S{cn}"
        plot_data = {"INVYR": inv_year, "LAT": stand.get("LATITUDE"), "LON": stand.get("LONGITUDE"),
                     "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
                     "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50}
        try:
            sdf = P.build_fvs_standinit(plot_data, sid, vl)
            tdf = treeinit_for_stand(fvs_rows, sid)
        except Exception as e:
            failures.append({"stand_cn": cn, "stage": "build", "detail": str(e)}); continue
        if tdf.empty:
            continue
        # calibrated point
        try:
            fr = P.run_fvs_projection(sdf, tdf, sid, vl, config_version="calibrated",
                                      num_cycles=20, cycle_length=5)
            for cy, tl in sorted(fr["treelists"].items()):
                py = cy - inv_year
                if py < 0:
                    continue
                rows.append({"STAND_CN": cn, "STATE": state, "YEAR": cy, "PROJ_YEAR": py,
                             "VARIANT": vu, "KIND": "point", "DRAW": -1,
                             "AGB_TONS_AC": round(float(P.compute_plot_agb(tl, nsbe)), 4)})
            n_point += 1
        except Exception as e:
            failures.append({"stand_cn": cn, "stage": "point", "detail": str(e)})
        # posterior draws (SDIMAX stripped), common indices across all stands
        for di in common_draws:
            try:
                draw = eng.get_draw(int(di))
                kw = strip_sdimax(eng.generate_keywords_for_draw(draw, dft, draw_idx=int(di)))
                fr = run_fvs_with_draw(sdf, tdf, sid, vl, draw_keywords=kw,
                                       num_cycles=20, cycle_length=5)
                for cy, tl in sorted(fr["treelists"].items()):
                    py = cy - inv_year
                    if py < 0:
                        continue
                    rows.append({"STAND_CN": cn, "STATE": state, "YEAR": cy, "PROJ_YEAR": py,
                                 "VARIANT": vu, "KIND": "draw", "DRAW": int(di),
                                 "AGB_TONS_AC": round(float(P.compute_plot_agb(tl, nsbe)), 4)})
                n_draw += 1
            except Exception as e:
                failures.append({"stand_cn": cn, "draw": int(di), "stage": "draw", "detail": str(e)})

    pd.DataFrame(rows).to_csv(done_csv, index=False)
    ledger = {"variant": vu, "k_stands": len(matched), "n_point_runs": n_point,
              "n_draw_runs": n_draw, "n_draws_req": a.n_draws, "n_rows": len(rows),
              "n_failures": len(failures), "elapsed_sec": round(time.time() - t0, 1),
              "failures": failures[:200]}
    with open(done_ledger, "w") as f:
        json.dump(ledger, f, indent=2)
    log.info(f"{vu} done: {n_point} point, {n_draw} draw runs, {len(rows)} rows, "
             f"{len(failures)} fail, {ledger['elapsed_sec']}s")


if __name__ == "__main__":
    main()
