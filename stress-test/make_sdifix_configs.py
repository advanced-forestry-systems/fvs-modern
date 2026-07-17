#!/usr/bin/env python3
# make_sdifix_configs.py - fix the NA max-SDI dropout in calibrated FVS configs.
#
# AUDIT finding (FVS_SDIMAX_AUDIT_2026-06-10): calibrated configs carry NA for many
# species' max SDI; the keyword writer skips NA species, so they silently revert to
# the FVS built-in variant DEFAULT (e.g. WS redwood -> 1052), over-allowing density
# (FVS plateaus near 85% of SDImax). Fix: fill each NA species with the variant's
# calibrated MEDIAN max SDI (a realized, FIA-consistent value) instead of letting it
# fall back to the high built-in default. Leaves calibrated (non-NA) values untouched.
# Writes corrected configs to config/calibrated_sdifix/<variant>.json.
import json, os, glob, statistics as st, sys
ROOT = sys.argv[1] if len(sys.argv)>1 else os.path.expanduser("~/fvs-modern")
SRC = os.path.join(ROOT,"config","calibrated")
DST = os.path.join(ROOT,"config","calibrated_sdifix"); os.makedirs(DST, exist_ok=True)

def find_holder(o, keys=("sdimax","SDICON")):
    # return (dict, key) that holds the sdimax list, searching recursively
    if isinstance(o, dict):
        for k,v in o.items():
            if k in keys and isinstance(v, list): return (o,k)
            r=find_holder(v,keys)
            if r: return r
    elif isinstance(o, list):
        for x in o:
            r=find_holder(x,keys)
            if r: return r
    return None

summary=[]
for f in sorted(glob.glob(os.path.join(SRC,"*.json"))):
    v=os.path.basename(f)[:-5]
    d=json.load(open(f)); h=find_holder(d)
    if not h: continue
    holder,key=h; arr=holder[key]
    nums=[float(x) for x in arr if str(x).upper() not in ("NA","NONE","")]
    if not nums: continue
    med=round(st.median(nums),1); na=0
    for i,x in enumerate(arr):
        if str(x).upper() in ("NA","NONE",""):
            arr[i]=med; na+=1
    if na:
        json.dump(d, open(os.path.join(DST,f"{v}.json"),"w"), indent=2)
        summary.append((v,na,len(arr),med))
print("variant  NA_filled  n_species  fill_value(median calib SDImax)")
for v,na,n,med in summary: print(f"{v:6}  {na:9}  {n:9}  {med}")
print(f"\nwrote {len(summary)} corrected configs -> {DST}")
