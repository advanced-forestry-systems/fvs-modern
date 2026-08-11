#!/usr/bin/env python3
"""
Verify the corrected per-species SDIMAX keyword format binds (item 2, 2026-06-16).

Root cause (from vbase/initre.f90 option 89): SDIMAX field 1 = species code
(SPDECD; 0/blank = all), field 2 = max SDI value. FVS reads fixed 10-col fields:
keyword cols 1-10, field1 cols 11-20, field2 cols 21-30. The old emitter
("SDIMAX" + 10 spaces + species(10) + value(10)) pushed the species into cols
17-26, so field1 (11-20) read blank -> ALL species, and field2 (21-30) caught the
trailing digit of the species index -> "MAX SDI = 1". Hence the over-thinning.

Correct emission: "%-10s%10d%10.1f" % ("SDIMAX", species, value).

Test: project NE stands 50 yr (5 cycles x 10 yr). Arms: default (no SDIMAX), and
corrected per-species SDIMAX at native x {0.6, 1.0, 1.6}. If lower max SDI yields
lower final TPH monotonically, the format binds. (Note FVS auto-raises a max set
below initial density, sdichk.f90, so the effect shows over a long projection as
stands grow into the limit.)
"""
import os, sys, math, tempfile, sqlite3

P = "/users/PUOM0008/crsfaaron/fvs-modern"
FIA = "/fs/scratch/PUOM0008/crsfaaron/FIA"
WORK = os.path.expanduser("~/overthin_work")
os.environ["FIA_DATA_DIR"] = FIA
os.environ["FVS_PROJECT_ROOT"] = P
os.environ["FVS_LIB_DIR"] = P + "/lib"
os.environ["FVS_CONFIG_DIR"] = P + "/config"
for p in [WORK, "/users/PUOM0008/crsfaaron/fvs-conus/python", P, P + "/calibration/python", P + "/calibration", P + "/deployment/fvs2py"]:
    sys.path.insert(0, p)

import pandas as pd
import numpy as np
from pathlib import Path
import fia_stand_generator as G
from config.config_loader import FvsConfigLoader
from perseus_100yr_projection import (
    build_fvs_treeinit, build_fvs_standinit, KEYFILE_TEMPLATE, FVS_LIB_DIR, _run_via_subprocess,
)

VAR = "ne"
STATES = "CT,ME,MA,NH,NY,RI,VT".split(",")
NSAMP = int(os.environ.get("NSAMP", "60"))
SEED = 11
CFG_DIR = P + "/config"
TPHc = 2.4710538
M2HA = 0.2296

_REV = {"CT":9,"ME":23,"MA":25,"NH":33,"NY":36,"RI":44,"VT":50}
G.VARIANT_STATES[VAR] = tuple(_REV.values())
_o = G._state_abbrev; _A = {v: k for k, v in _REV.items()}
G._state_abbrev = lambda c: _A.get(c) or _o(c)

loader = FvsConfigLoader(VAR, version="calibrated", config_dir=CFG_DIR)
native_sdi = loader._find_sdi_param(loader.config.get("categories", {})) or []
maxsp = loader.config.get("maxsp", 0)


def sdimax_block(scale):
    # CORRECT format: keyword(10) + species(10) + value(10)
    lines = []
    for i, v in enumerate(native_sdi):
        if isinstance(v, (int, float)) and v > 0:
            lines.append("%-10s%10d%10.1f" % ("SDIMAX", i + 1, float(v) * scale))
    return "\n".join(lines)


def run(std, tdf, sid, kw, ncyc=5, clen=10):
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
        return _run_via_subprocess(os.path.join(FVS_LIB_DIR, "FVSne"), kp, db, tmp)


def g(r, k):
    try: return float(r.get(k, 0) or 0)
    except Exception: return 0.0


def main():
    arms = {"default": "", "sdi_x0.6": sdimax_block(0.6),
            "sdi_x1.0": sdimax_block(1.0), "sdi_x1.6": sdimax_block(1.6)}
    print(f"maxsp={maxsp} NSAMP={NSAMP}; native SDI head={[round(float(x),0) for x in native_sdi[:5]]}")

    frames = []
    for ab in STATES:
        f = Path(FIA) / f"{ab}_PLOT.csv"
        if f.exists():
            frames.append(pd.read_csv(f, usecols=lambda c: c in ("CN", "PREV_PLT_CN", "MEASYEAR", "STATECD"), low_memory=False))
    plot = pd.concat(frames, ignore_index=True)
    rem = plot.dropna(subset=["PREV_PLT_CN"]).copy()
    rem["CN"] = rem.CN.astype("int64")
    rem = rem.sample(n=min(NSAMP, len(rem)), random_state=SEED)
    tr = G.load_fia_trees(VAR, rem.CN.tolist(), Path(FIA))

    finals = {k: [] for k in arms}
    n = 0
    for _, r in rem.iterrows():
        cn = int(r.CN); t = tr[tr.PLT_CN == cn]
        if len(t) < 5:
            continue
        sid = str(cn)
        std = build_fvs_standinit({"INVYR": 2000, "STATECD": int(r.STATECD), "COUNTYCD": 0}, sid, VAR)
        std["inv_plot_size"] = 1.0; std["brk_dbh"] = 99.0
        tdf = build_fvs_treeinit(t, sid)
        n += 1
        for lab, kw in arms.items():
            try:
                res = run(std, tdf, sid, kw)
                s = res.get("summary")
                if s is None or len(s) == 0:
                    continue
                finals[lab].append(g(s.iloc[-1], "Tpa") * TPHc)
            except Exception:
                pass
    print(f"n_run={n}")
    print(f"{'arm':<12}{'mean final TPH':>16}{'median':>10}")
    for lab in ["default", "sdi_x0.6", "sdi_x1.0", "sdi_x1.6"]:
        v = finals[lab]
        if v:
            print(f"{lab:<12}{np.mean(v):>16.1f}{np.median(v):>10.1f}  (n={len(v)})")
    print("BINDING" if (finals['sdi_x0.6'] and finals['sdi_x1.6'] and
                        np.mean(finals['sdi_x0.6']) < np.mean(finals['sdi_x1.6'])) else "NO_CLEAR_BINDING")
    print("DONE_BIND")


if __name__ == "__main__":
    main()
