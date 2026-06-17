#!/usr/bin/env python3
"""
Combined density+growth benchmark (2026-06-16). Three arms per stand, one cycle at
the observed remeasurement interval, compared to the observed later measurement:
  default     : native FVS, no calibration
  calibrated  : per-species mortality/DG/crown/HD multipliers (SDIMAX suppressed)
  cal_locSDI  : calibrated multipliers PLUS the per-stand localized FIA max SDI
                (brms_SDImax) emitted with the corrected keyword format, at the
                co-calibrated per-variant level.
Reports BA, TPH, QMD bias and RMSE vs observed. Env: VAR, STATES, MAXSP, LEVEL, NSAMP, SEED, OUTCSV.
"""
import os, sys, math, tempfile, sqlite3

P = "/users/PUOM0008/crsfaaron/fvs-modern"
FIA = "/fs/scratch/PUOM0008/crsfaaron/FIA"
CONUS = "/users/PUOM0008/crsfaaron/fvs-conus"
WORK = os.path.expanduser("~/overthin_work")
os.environ["FIA_DATA_DIR"] = FIA
os.environ["FVS_PROJECT_ROOT"] = P
os.environ["FVS_LIB_DIR"] = P + "/lib"
os.environ["FVS_CONFIG_DIR"] = P + "/config"
for p in [WORK, CONUS + "/python", P, P + "/calibration/python"]:
    sys.path.insert(0, p)

import pandas as pd
import numpy as np
from pathlib import Path
import fia_stand_generator as G
from config.config_loader import FvsConfigLoader
from perseus_100yr_projection import (
    build_fvs_treeinit, build_fvs_standinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess,
)

VAR = os.environ.get("VAR", "ne")
STATES = os.environ.get("STATES", "CT,ME,MA,NH,NY,RI,VT").split(",")
MAXSP = int(os.environ.get("MAXSP", "108"))
LEVEL = float(os.environ.get("LEVEL", "1.8"))   # co-calibrated localized max-SDI scale
NSAMP = int(os.environ.get("NSAMP", "120"))
SEED = int(os.environ.get("SEED", "5"))
OUTCSV = os.environ.get("OUTCSV", os.path.expanduser(f"~/overthin_work/comb_{VAR}.csv"))
CFG_DIR = P + "/config"
M2HA = 0.2296; TPHc = 2.4710538; CMc = 2.54

_REV = {"AL":1,"CA":6,"CO":8,"CT":9,"FL":12,"GA":13,"ID":16,"IL":17,"IN":18,"IA":19,"ME":23,"MA":25,
        "MI":26,"MN":27,"MS":28,"MO":29,"MT":30,"NV":32,"NH":33,"NM":35,"NY":36,"OR":41,"RI":44,
        "SC":45,"TN":47,"UT":49,"VT":50,"WA":53,"WI":55,"WY":56}
G.VARIANT_STATES[VAR] = tuple(_REV[x] for x in STATES if x in _REV)
_o = G._state_abbrev; _A = {v: k for k, v in _REV.items()}
G._state_abbrev = lambda c: _A.get(c) or _o(c)

# calibrated keyword block (clean loader; SDIMAX already suppressed)
CAL = FvsConfigLoader(VAR, version="calibrated", config_dir=CFG_DIR).generate_keywords(include_comments=False)

# localized max SDI per stand: correct keyword format (keyword cols 1-10, species 11-20, value 21-30)
def locsdi_kw(sdi_value):
    v = sdi_value / TPHc  # brms metric -> FVS imperial SDI units (matches var_scale_diag)
    return "\n".join("%-10s%10d%10.1f" % ("SDIMAX", i, v) for i in range(1, MAXSP + 1))


def run_kw(std, tdf, sid, kw, ncyc=1, clen=10):
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "FVS_Data.db")
        con = sqlite3.connect(db)
        std.to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tdf.to_sql("fvs_treeinit", con, if_exists="replace", index=False)
        con.close()
        kc = KEYFILE_TEMPLATE.format(stand_id=sid, db_path=db,
                                     calibration_keywords=kw if kw else "** DEFAULT PARAMETERS",
                                     num_cycles=ncyc, cycle_length=clen)
        kp = os.path.join(tmp, "t.key"); open(kp, "w").write(kc)
        return _run_via_subprocess(os.path.join(FVS_LIB_DIR, "FVS" + VAR), kp, db, tmp)


def metr(df):
    d = df[(df.DIA > 0) & (df.TPA_UNADJ > 0)]
    if len(d) == 0:
        return (0, 0, 0)
    return ((d.TPA_UNADJ * 0.005454 * d.DIA ** 2).sum() * M2HA,
            d.TPA_UNADJ.sum() * TPHc,
            math.sqrt((d.DIA ** 2 * d.TPA_UNADJ).sum() / d.TPA_UNADJ.sum()) * CMc)


def g(r, k):
    try: return float(r.get(k, 0) or 0)
    except Exception: return 0.0


def stt(p, o):
    m = [(a, b) for a, b in zip(p, o) if a == a and b == b]
    k = len(m)
    if k == 0: return (float("nan"), float("nan"), 0)
    e = [a - b for a, b in m]; mo = sum(b for _, b in m) / k
    return (100 * math.sqrt(sum(x * x for x in e) / k) / mo, 100 * sum(e) / k / mo, k)


def main():
    print(f"VAR={VAR} MAXSP={MAXSP} LEVEL={LEVEL} NSAMP={NSAMP}")
    cols = ["CN", "PREV_PLT_CN", "MEASYEAR", "STATECD", "UNITCD", "COUNTYCD", "PLOT"]
    frames = [pd.read_csv(Path(FIA) / f"{ab}_PLOT.csv", usecols=lambda c: c in cols, low_memory=False)
              for ab in STATES if (Path(FIA) / f"{ab}_PLOT.csv").exists()]
    plot = pd.concat(frames, ignore_index=True)
    plot["plot_key"] = (plot.STATECD.astype(int).astype(str) + "-" + plot.UNITCD.astype(int).astype(str)
                        + "-" + plot.COUNTYCD.astype(int).astype(str) + "-" + plot.PLOT.astype(int).astype(str))
    cn2key = dict(zip(plot.CN.astype("int64"), plot.plot_key))
    yr = dict(zip(plot.CN.astype("int64"), plot.MEASYEAR))
    b = pd.read_csv(CONUS + "/data/brms_SDImax.csv")
    b.columns = [c.strip().strip('"') for c in b.columns]
    key2sdi = dict(zip(b.ID.astype(str), b["SDImax.mean"]))

    rem = plot.dropna(subset=["PREV_PLT_CN"]).copy()
    rem["PREV_PLT_CN"] = rem.PREV_PLT_CN.astype("int64"); rem["CN"] = rem.CN.astype("int64")
    rem["interval"] = rem.apply(lambda r: r.MEASYEAR - yr.get(r.PREV_PLT_CN, np.nan), axis=1)
    rem = rem[(rem.interval >= 5) & (rem.interval <= 15)]
    rem = rem[rem.plot_key.map(lambda k: k in key2sdi)]
    rem = rem.sample(n=min(NSAMP, len(rem)), random_state=SEED)

    tr1 = G.load_fia_trees(VAR, rem.PREV_PLT_CN.tolist(), Path(FIA))
    tr2 = G.load_fia_trees(VAR, rem.CN.tolist(), Path(FIA))

    arms = ["default", "calibrated", "cal_locSDI"]
    A = {a: {x: [] for x in ["BA", "TPH", "QMD", "oBA", "oTPH", "oQMD"]} for a in arms}
    rows = []; n = 0
    for _, r in rem.iterrows():
        t1 = int(r.PREV_PLT_CN); tr = tr1[tr1.PLT_CN == t1]; tro = tr2[tr2.PLT_CN == int(r.CN)]
        if len(tr) < 5 or len(tro) < 3: continue
        brms = key2sdi.get(cn2key.get(t1))
        if brms is None or not np.isfinite(brms): continue
        oBA, oTPH, oQMD = metr(tro)
        if oBA <= 0: continue
        sid = str(t1); yrs = int(r.interval)
        std = build_fvs_standinit({"INVYR": 2000, "STATECD": int(r.STATECD), "COUNTYCD": 0}, sid, VAR)
        std["inv_plot_size"] = 1.0; std["brk_dbh"] = 99.0
        tdf = build_fvs_treeinit(tr, sid)
        n += 1
        kwmap = {"default": "", "calibrated": CAL, "cal_locSDI": CAL + "\n" + locsdi_kw(brms * LEVEL)}
        rec = {"plt": t1, "interval": yrs, "brms": brms}
        for a in arms:
            try:
                res = run_kw(std, tdf, sid, kwmap[a], ncyc=1, clen=yrs)
                s = res.get("summary")
                if s is None or len(s) == 0: continue
                l = s.iloc[-1]; ba = g(l, "BA") * M2HA
                if ba <= 0: continue
                A[a]["BA"].append(ba); A[a]["TPH"].append(g(l, "Tpa") * TPHc); A[a]["QMD"].append(g(l, "QMD") * CMc)
                A[a]["oBA"].append(oBA); A[a]["oTPH"].append(oTPH); A[a]["oQMD"].append(oQMD)
                rec[f"{a}_BA"] = ba; rec[f"{a}_TPH"] = g(l, "Tpa") * TPHc
            except Exception:
                pass
        rows.append(rec)
    pd.DataFrame(rows).to_csv(OUTCSV, index=False)
    print(f"n_run={n}")
    print(f"{'arm':<12}{'BA bias%':>9}{'BA RMSE%':>9}{'TPH bias%':>10}{'TPH RMSE%':>10}{'QMD bias%':>10}{'n':>5}")
    for a in arms:
        ba = stt(A[a]["BA"], A[a]["oBA"]); t = stt(A[a]["TPH"], A[a]["oTPH"]); q = stt(A[a]["QMD"], A[a]["oQMD"])
        print(f"{a:<12}{ba[1]:>9.1f}{ba[0]:>9.1f}{t[1]:>10.1f}{t[0]:>10.1f}{q[1]:>10.1f}{t[2]:>5d}")
    print("DONE_COMBINED")


if __name__ == "__main__":
    main()
