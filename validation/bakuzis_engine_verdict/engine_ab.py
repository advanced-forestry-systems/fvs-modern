#!/usr/bin/env python3
"""Engine-level A/B of gompmort survival forms.
Runs the SAME standard NE stand ~100 yr (20 cycles x 5 yr) under two libs:
  GOMPIT lib   (S_ann = 1 - exp(-exp(ETA)))
  EXP-HAZ lib  (HZ = exp(ETA); S = exp(-HZ*FINTL))
Both with FVS_GOMPIT=1 and the SAME Greg coefficient CSV, so the only
difference is the survival transform compiled into the .so.
Extracts TPA and BA trajectory from fvs.summary.
"""
import os, sys, tempfile, json
import numpy as np, pandas as pd
ROOT=os.path.expanduser("~/fvs-modern")
sys.path.insert(0, os.path.join(ROOT,"deployment","fvs2py"))
from fvs2py import FVS

LIB=os.environ["FVSNE_SO"]
CFG=os.path.join(ROOT,"config")
TAG=os.environ.get("AB_TAG","run")
OUT=os.environ.get("AB_OUT","/fs/scratch/PUOM0008/crsfaaron/wt-engine/ab_engine")
os.makedirs(OUT,exist_ok=True)

# Build a standard, fully-stocked even-aged NE stand: red spruce (SPCD 97) +
# balsam fir (SPCD 12) mix, ~600 TPA, small-medium diameters so density-driven
# mortality has room to act over 100 yr. FIA/NE tree-record format used by NE.key.
# TREEFMT columns (from NE.key): plot(I4) tree(I4) count(F6) hist(I1) spp(A3) dbh(F4.1)...
def stand_key():
    lines=[]
    lines+=["STDIDENT","ABSTD ENGINE-AB","STDINFO          922                   1",
            "DESIGN            -1         1","INVYEAR       1990.0","NUMCYCLE        20.0",
            "TIMEINT            0         5",
            "TREEFMT","(I4,I4,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,","6I2,2I1,I2,2I3,2I1,F3.0)",
            "TREEDATA        15.0"]
    # generate ~40 records spanning a diameter distribution, expansion via count
    rng=np.random.default_rng(7)
    recs=[]
    tno=1
    # inverse-J: many small, few large; total ~600 TPA over the records
    diam_bins=[(1.5,120),(2.5,110),(3.5,95),(4.5,80),(5.5,60),(6.5,45),(7.5,30),(8.5,20),(9.5,12),(11.0,8)]
    for (dbh,tpa) in diam_bins:
        for sp in ["097","012"]:   # red spruce, balsam fir (FIA codes as A3)
            cnt=tpa/2.0
            # crown ratio ~ 0.5 default; leave blank -> FVS estimates. Provide dbh only.
            plot=101
            line=f"{plot:4d}{tno:4d}{cnt:6.1f}0{sp:>3}{dbh:4.1f}"
            recs.append(line); tno+=1
    lines+=recs
    lines+=["-999","PROCESS","STOP"]
    return "\n".join(lines)+"\n"

def main():
    with tempfile.TemporaryDirectory() as tmp:
        kp=os.path.join(tmp,"ab.key"); open(kp,"w").write(stand_key())
        os.chdir(tmp)
        fvs=FVS(lib_path=LIB, config_version=None, config_dir=CFG)
        fvs.load_keyfile(kp)
        fvs.run()
        s=fvs.summary
        if s is None or len(s)==0:
            print(f"[{TAG}] SUMMARY EMPTY"); return
        print(f"[{TAG}] summary cols: {list(s.columns)}")
        s.to_csv(os.path.join(OUT,f"summary_{TAG}.csv"),index=False)
        # print key trajectory
        cols={c.lower():c for c in s.columns}
        yr = cols.get("year") or cols.get("inv_year") or list(s.columns)[0]
        tpa = cols.get("tpa") or cols.get("live_tpa")
        ba  = cols.get("ba") or cols.get("baa")
        print(f"[{TAG}] yr_col={yr} tpa_col={tpa} ba_col={ba}")
        show=[c for c in [yr,tpa,ba] if c]
        print(s[show].to_string(index=False))

if __name__=="__main__":
    main()
