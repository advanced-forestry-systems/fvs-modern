#!/usr/bin/env python3
# aggregate_posterior_ci.py - combine the per-draw FVS posterior files
# (post_<ST>/posterior_<ST>_d<i>.csv: ST,variant,draw,year,mean_density,n_plots)
# into the relative parameter-uncertainty CI per state x year:
#   rel_lo = p2.5(mean_density)/median, rel_hi = p97.5/median, anchored so 2025=1.
# Writes posterior_ci_all.csv (ST,year,mean,rel_lo,rel_hi,ndraw).
import csv, glob, os, statistics as st
FV="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
def pctl(v,p):
    v=sorted(v); k=(len(v)-1)*p; f=int(k);
    return v[f] if f+1>=len(v) else v[f]+(v[f+1]-v[f])*(k-f)
by={}   # (ST,year) -> list of mean_density across draws
for f in glob.glob(f"{FV}/post_*/posterior_*_d*.csv"):
    for r in csv.DictReader(open(f)):
        try: by.setdefault((r["ST"],int(float(r["year"]))),[]).append(float(r["mean_density"]))
        except: pass
rows=[]
states=sorted({k[0] for k in by})
for ST in states:
    yrs=sorted({y for (s,y) in by if s==ST})
    for y in yrs:
        v=by[(ST,y)]
        if len(v)<3: continue
        med=st.median(v); lo=pctl(v,0.025); hi=pctl(v,0.975)
        rows.append((ST,y,round(med,3),round(lo/med,4) if med else 1,round(hi/med,4) if med else 1,len(v)))
with open(f"{FV}/posterior_ci_all.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(["ST","year","mean","rel_lo","rel_hi","ndraw"]); w.writerows(rows)
print(f"{len(states)} states -> posterior_ci_all.csv: {states}")
# widest parameter band per state at its last year (sanity)
for ST in states:
    ys=[r for r in rows if r[0]==ST];
    if ys: last=ys[-1]; print(f"  {ST}: {last[1]} rel [{last[3]},{last[4]}] (+/-{round(50*(last[4]-last[3]),1)}%) n={last[5]}")
