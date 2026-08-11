#!/usr/bin/env python3
"""
WO-1 over-thinning isolation diagnostic (recreated 2026-06-16).

Question: the calibrated fvs-modern config improves basal area but over-thins
density (TPH falls too far). Which emission causes it: the explicit SDIMAX
(density-limit) keyword block, or the MORTMULT (mortality) multipliers?

Method: four arms on the same FIA remeasurement pairs, projected one cycle of
length = observed interval, compared to the observed later measurement.
  1. default            : no calibration keywords (FVS internal behavior)
  2. calibrated         : full generate_keywords() block (production calibrated)
  3. cal_noSDI(9999)    : calibrated, but SDIMAX neutralized (set 9999 all spp)
  4. cal_noMORT         : calibrated, but MORTMULT lines removed

If arm 3 returns TPH bias toward the default arm while arm 2 over-thins, the
SDIMAX emission is the driver. If arm 4 does, mortality multipliers are.

Runs entirely through the subprocess FVS executable (works under py3.9). Does
NOT import the in-repo conflicted config_loader; uses a clean standalone copy
placed at ~/overthin_work/config/ by the launcher. No repo/git writes.

Env knobs: VAR, STATES, NSAMP, SEED, OUTCSV.
"""
import os, sys, math, re, tempfile, sqlite3

P = "/users/PUOM0008/crsfaaron/fvs-modern"
FIA = "/fs/scratch/PUOM0008/crsfaaron/FIA"
WORK = os.path.expanduser("~/overthin_work")
os.environ["FIA_DATA_DIR"] = FIA
os.environ["FVS_PROJECT_ROOT"] = P
os.environ["FVS_LIB_DIR"] = P + "/lib"
os.environ["FVS_CONFIG_DIR"] = P + "/config"

# WORK first so the clean config package shadows the repo's conflicted one
for p in [WORK, "/users/PUOM0008/crsfaaron/fvs-conus/python", P, P + "/calibration/python", P + "/calibration", P + "/deployment/fvs2py"]:
    sys.path.insert(0, p)

import pandas as pd
import numpy as np
from pathlib import Path
import fia_stand_generator as G
from config.config_loader import FvsConfigLoader  # clean copy from WORK
from perseus_100yr_projection import (
    build_fvs_treeinit, build_fvs_standinit,
    KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess,
)

VAR = os.environ.get("VAR", "ne")
STATES = os.environ.get("STATES", "CT,ME,MA,NH,NY,RI,VT").split(",")
NSAMP = int(os.environ.get("NSAMP", "150"))
SEED = int(os.environ.get("SEED", "11"))
OUTCSV = os.environ.get("OUTCSV", os.path.expanduser(f"~/overthin_work/overthin_{VAR}.csv"))
CFG_DIR = P + "/config"

M2HA = 0.2296      # ft2/ac -> m2/ha
TPHc = 2.4710538   # tpa -> tph
CMc = 2.54         # in -> cm

_REV = {"AL":1,"CA":6,"CO":8,"CT":9,"FL":12,"GA":13,"ID":16,"IL":17,"IN":18,"IA":19,
        "ME":23,"MA":25,"MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,
        "NM":35,"NY":36,"OR":41,"RI":44,"SC":45,"TN":47,"UT":49,"VT":50,"WA":53,
        "WI":55,"WY":56}
G.VARIANT_STATES[VAR] = tuple(_REV[x] for x in STATES if x in _REV)
_o = G._state_abbrev
_A = {v: k for k, v in _REV.items()}
G._state_abbrev = lambda c: _A.get(c) or _o(c)


def build_arms(var):
    """Return dict arm_label -> calibration keyword string.

    Arms:
      default     : no calibration keywords
      calibrated  : current production block (SDIMAX emitted species-first = BUG)
      calfix      : same block but SDIMAX re-emitted value-first (proposed fix)
    """
    loader = FvsConfigLoader(var.lower(), version="calibrated", config_dir=CFG_DIR)
    full = loader.generate_keywords(include_comments=False)
    maxsp = loader.config.get("maxsp", 0)
    sdi = loader._find_sdi_param(loader.config.get("categories", {})) or []

    blocks = {}
    for kw in ("SDIMAX", "BAMAX", "MORTMULT", "BAIMULT", "HTGMULT"):
        blocks[kw] = sum(1 for ln in full.splitlines() if ln.strip().startswith(kw))

    # calfix: strip the buggy SDIMAX lines, re-emit value-first (field1=value,
    # field2=species index) which the engine parses correctly. All other
    # calibrated keywords (MORTMULT/BAIMULT/HTGMULT) are preserved verbatim.
    non_sdi = [ln for ln in full.splitlines() if not ln.strip().startswith("SDIMAX")]
    sdifix = "\n".join(
        f"SDIMAX    {float(v):10.1f}{i + 1:10d}"
        for i, v in enumerate(sdi)
        if isinstance(v, (int, float)) and v > 0
    )
    calfix = "!! SDIMAX value-first (corrected field order)\n" + sdifix + "\n" + "\n".join(non_sdi)

    arms = {
        "default": "** DEFAULT PARAMETERS",
        "calibrated": full,
        "calfix": calfix,
    }
    return arms, blocks, maxsp


def run_with_kw(std, tdf, sid, var, kw, yrs):
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "FVS_Data.db")
        con = sqlite3.connect(db)
        std.to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tdf.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        kc = KEYFILE_TEMPLATE.format(
            stand_id=sid, db_path=db,
            calibration_keywords=kw if kw else "** DEFAULT PARAMETERS",
            num_cycles=1, cycle_length=yrs,
        )
        kp = os.path.join(tmp, f"{var}_{sid}.key")
        with open(kp, "w") as f:
            f.write(kc)
        exe = os.path.join(FVS_LIB_DIR, f"FVS{var.lower()}")
        return _run_via_subprocess(exe, kp, db, tmp)


def metr(df):
    d = df[(df.DIA > 0) & (df.TPA_UNADJ > 0)]
    if len(d) == 0:
        return (0, 0, 0)
    return ((d.TPA_UNADJ * 0.005454 * d.DIA ** 2).sum() * M2HA,
            d.TPA_UNADJ.sum() * TPHc,
            math.sqrt((d.DIA ** 2 * d.TPA_UNADJ).sum() / d.TPA_UNADJ.sum()) * CMc)


def g(r, k):
    try:
        return float(r.get(k, 0) or 0)
    except Exception:
        return 0.0


def stt(p, o):
    m = [(a, b) for a, b in zip(p, o) if a == a and b == b]
    k = len(m)
    if k == 0:
        return (float("nan"), float("nan"), 0)
    e = [a - b for a, b in m]
    mo = sum(b for _, b in m) / k
    return (100 * math.sqrt(sum(x * x for x in e) / k) / mo, 100 * sum(e) / k / mo, k)


def main():
    arms, blocks, maxsp = build_arms(VAR)
    print(f"VAR={VAR} maxsp={maxsp} NSAMP={NSAMP} SEED={SEED}", flush=True)
    print(f"keyword blocks in calibrated: {blocks}", flush=True)

    frames = []
    for ab in STATES:
        f = Path(FIA) / f"{ab}_PLOT.csv"
        if f.exists():
            frames.append(pd.read_csv(
                f, usecols=lambda c: c in ("CN", "PREV_PLT_CN", "MEASYEAR", "STATECD"),
                low_memory=False))
    plot = pd.concat(frames, ignore_index=True)
    yr = dict(zip(plot.CN.astype("int64"), plot.MEASYEAR))
    rem = plot.dropna(subset=["PREV_PLT_CN"]).copy()
    rem["PREV_PLT_CN"] = rem.PREV_PLT_CN.astype("int64")
    rem["CN"] = rem.CN.astype("int64")
    rem["interval"] = rem.apply(lambda r: r.MEASYEAR - yr.get(r.PREV_PLT_CN, np.nan), axis=1)
    rem = rem[(rem.interval >= 5) & (rem.interval <= 15)].sample(
        n=min(NSAMP, len(rem)), random_state=SEED)

    tr1 = G.load_fia_trees(VAR, rem.PREV_PLT_CN.tolist(), Path(FIA))
    tr2 = G.load_fia_trees(VAR, rem.CN.tolist(), Path(FIA))

    A = {lab: {x: [] for x in ["BA", "TPH", "QMD", "oBA", "oTPH", "oQMD"]} for lab in arms}
    rows = []
    n = 0
    for _, r in rem.iterrows():
        t1 = int(r.PREV_PLT_CN)
        tr = tr1[tr1.PLT_CN == t1]
        tro = tr2[tr2.PLT_CN == int(r.CN)]
        if len(tr) < 5 or len(tro) < 3:
            continue
        oBA, oTPH, oQMD = metr(tro)
        if oBA <= 0:
            continue
        sid = str(t1)
        yrs = int(r.interval)
        std = build_fvs_standinit({"INVYR": 2000, "STATECD": int(r.STATECD), "COUNTYCD": 0}, sid, VAR)
        std["inv_plot_size"] = 1.0
        std["brk_dbh"] = 99.0
        tdf = build_fvs_treeinit(tr, sid)
        n += 1
        rec = {"plt": t1, "interval": yrs, "oBA": oBA, "oTPH": oTPH, "oQMD": oQMD}
        for lab, kw in arms.items():
            try:
                res = run_with_kw(std, tdf, sid, VAR, kw, yrs)
                s = res.get("summary")
                if s is None or len(s) == 0:
                    continue
                l = s.iloc[-1]
                ba = g(l, "BA") * M2HA
                if ba <= 0:
                    continue
                tph = g(l, "Tpa") * TPHc
                qmd = g(l, "QMD") * CMc
                a = A[lab]
                a["BA"].append(ba); a["TPH"].append(tph); a["QMD"].append(qmd)
                a["oBA"].append(oBA); a["oTPH"].append(oTPH); a["oQMD"].append(oQMD)
                rec[f"{lab}_BA"] = ba; rec[f"{lab}_TPH"] = tph; rec[f"{lab}_QMD"] = qmd
            except Exception:
                pass
        rows.append(rec)
        if n % 25 == 0:
            print(f"  ...{n} stands done", flush=True)

    pd.DataFrame(rows).to_csv(OUTCSV, index=False)
    print(f"\nVAR={VAR} n_run={n}  (raw -> {OUTCSV})", flush=True)
    print(f"{'arm':<14}{'BA bias%':>10}{'BA RMSE%':>10}{'TPH bias%':>11}{'TPH RMSE%':>11}{'QMD bias%':>11}{'n':>6}")
    for lab in ["default", "calibrated", "calfix"]:
        b = stt(A[lab]["BA"], A[lab]["oBA"])
        t = stt(A[lab]["TPH"], A[lab]["oTPH"])
        q = stt(A[lab]["QMD"], A[lab]["oQMD"])
        print(f"{lab:<14}{b[1]:>10.1f}{b[0]:>10.1f}{t[1]:>11.1f}{t[0]:>11.1f}{q[1]:>11.1f}{t[2]:>6d}")
    print("DONE_OVERTHIN", flush=True)


if __name__ == "__main__":
    main()
