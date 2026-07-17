#!/usr/bin/env python3
"""Empirical SDIMAX field-order test for arm 3. Runs FVS-NE on a few matched
stands under: (A) SDIMAX off (WO-1 base), (B) draw species-first (as-is, suspect
buggy), (C) draw value-first (proposed fix). Reports final-cycle AGB. If B is
materially below A and C, the species-first order over-thins (bug confirmed)."""
import os,sys,re,glob
import numpy as np, pandas as pd
sys.path.insert(0,"/users/PUOM0008/crsfaaron/fvs-modern")
sys.path.insert(0,"/users/PUOM0008/crsfaaron/fvs-conus/python")
import perseus_100yr_projection as P
from perseus_uncertainty_projection import run_fvs_with_draw
from config.uncertainty import UncertaintyEngine
from config.config_loader import FvsConfigLoader

VAR="ne"; STATE="NH"  # NE variant, NH state
SI="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant/standinit_NE.csv"
TI="/fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h/NH_FVS_TREEINIT_PLOT.csv"
FIPS_NH=33
si=pd.read_csv(SI,low_memory=False)
si["STAND_CN"]=si["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
si=si[si["STATE"].apply(lambda x:int(float(x))==FIPS_NH)]
tt=pd.read_csv(TI,low_memory=False)
tt["STAND_CN"]=tt["STAND_CN"].apply(lambda x:str(int(float(x))) if pd.notna(x) else "")
by=dict(tuple(tt.groupby("STAND_CN")))
nsbe=P.NSBECalculator(P.NSBE_ROOT)
eng=UncertaintyEngine(VAR,config_dir="/users/PUOM0008/crsfaaron/fvs-modern/config",seed=42)
dft=FvsConfigLoader(VAR,version="default",config_dir="/users/PUOM0008/crsfaaron/fvs-modern/config").config

def swap_sdimax(block):
    out=[]
    for ln in block.splitlines():
        m=re.match(r"^(SDIMAX\s+)(\S+)(\s+)(\S+)\s*$",ln)
        if m:
            # m groups: kw, f1(species_idx), spaces, f2(value) -> swap to value-first
            f1=m.group(2); f2=m.group(4)
            out.append(f"SDIMAX          {float(f2):10.1f}{int(float(f1)):10d}")
        else:
            out.append(ln)
    return "\n".join(out)

def finalAGB(sdf,tdf,sid,kw=None,cfg=None):
    if kw is None:
        fr=P.run_fvs_projection(sdf,tdf,sid,VAR,config_version=cfg,num_cycles=20,cycle_length=5)
    else:
        fr=run_fvs_with_draw(sdf,tdf,sid,VAR,draw_keywords=kw,num_cycles=20,cycle_length=5)
    yrs=sorted(fr["treelists"].keys())
    return P.compute_plot_agb(fr["treelists"][yrs[-1]],nsbe)

draw=eng.get_draw(0); kw_sp=eng.generate_keywords_for_draw(draw,dft,draw_idx=0); kw_val=swap_sdimax(kw_sp)
n=0
print("stand   A_off   B_spfirst  C_valfirst")
for cn,stand in si.iterrows():
    cnid=stand["STAND_CN"]; rows=by.get(cnid)
    if rows is None or rows.empty: continue
    sid=f"S{cnid}"
    iv=int(float(stand.get("INV_YEAR") or 2010))
    pd_={"INVYR":iv,"LAT":stand.get("LATITUDE"),"LON":stand.get("LONGITUDE"),"ELEV":stand.get("ELEVFT") or 500,"SLOPE":stand.get("SLOPE") or 10,"ASPECT":stand.get("ASPECT") or 180,"STDAGE":stand.get("AGE") or 50}
    try:
        sdf=P.build_fvs_standinit(pd_,sid,VAR)
        from run_conus_task_wo1 import treeinit_for_stand
        tdf=treeinit_for_stand(rows,sid)
        if tdf.empty: continue
        a=finalAGB(sdf,tdf,sid,cfg="calibrated")
        b=finalAGB(sdf,tdf,sid,kw=kw_sp)
        c=finalAGB(sdf,tdf,sid,kw=kw_val)
        print(f"{cnid[:8]} {a:8.2f} {b:9.2f} {c:10.2f}")
        n+=1
    except Exception as e:
        print("ERR",cnid[:8],str(e)[:80])
    if n>=3: break
print("done n=",n)
