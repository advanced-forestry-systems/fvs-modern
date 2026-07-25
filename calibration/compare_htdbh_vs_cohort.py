#!/usr/bin/env python3
"""Compare Wykoff HT-DBH scorecard vs cohort baseline: 5-model bias + deltas."""
import numpy as np, pandas as pd, sys

BASE = "/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_cohort/scorecard_ne_5model_cohort.csv"
HTDB = "/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_htdbh/scorecard_ne_5model_htdbh.csv"
MODELS = ["fvs_base","fvs_regional","organon","conus_spdep","conus_climate"]
METRICS = [("BA","BA_PRED","BA_OBS"),("MCuFt","MCuFt_PRED","MCuFt_OBS"),("BdFt","BdFt_PRED","BdFt_OBS")]

def stats(g,p,o):
    ok = g[p].notna()&g[o].notna()&(g[o]>0)&(g[p]>0)
    n=int(ok.sum())
    if n==0: return n,float('nan'),float('nan')
    P=g.loc[ok,p].astype(float); O=g.loc[ok,o].astype(float)
    return n, float(np.median(100*(P/O-1))), float(100*(P.mean()/O.mean()-1))

def tbl(df,label):
    out={}
    for m in MODELS:
        g=df[df["model"]==m]
        row={}
        for name,p,o in METRICS:
            row[name]=stats(g,p,o)
        out[m]=(len(g),row)
    return out

b=pd.read_csv(BASE); h=pd.read_csv(HTDB)
print(f"rows: cohort={len(b)}  htdbh={len(h)}")
print(f"cohort FORTYPE counts: {b['FORTYPE'].value_counts().to_dict()}")
print(f"htdbh  FORTYPE counts: {h['FORTYPE'].value_counts().to_dict()}")

for strat in ["OVERALL","SW","MW","HW"]:
    bb = b if strat=="OVERALL" else b[b["FORTYPE"]==strat]
    hh = h if strat=="OVERALL" else h[h["FORTYPE"]==strat]
    tb=tbl(bb,strat); th=tbl(hh,strat)
    print("\n"+"="*118)
    print(f"  STRATUM = {strat}   (cohort n_rows={len(bb)}, htdbh n_rows={len(hh)})")
    print("="*118)
    hdr=f"  {'model':<14}{'n':>5} | "+ " | ".join(f"{nm}_med%(coh/htd/Δ)  {nm}_agg%(coh/htd/Δ)" for nm,_,_ in METRICS)
    print("  metric blocks: cohort / htdbh / delta(htdbh-cohort)")
    for m in MODELS:
        nb,rb=tb[m]; nh,rh=th[m]
        segs=[]
        for nm,_,_ in METRICS:
            _,bm,ba=rb[nm]; _,hm,ha=rh[nm]
            dm=hm-bm; da=ha-ba
            segs.append(f"{nm}: med {bm:+6.1f}/{hm:+6.1f}/{dm:+5.1f}  agg {ba:+6.1f}/{ha:+6.1f}/{da:+5.1f}")
        print(f"  {m:<14}{nh:>5} | "+" | ".join(segs))
