import csv, os, collections
SRC="/fs/scratch/PUOM0008/crsfaaron/FIA/ENTIRE_FVS_STANDINIT_PLOT.csv"
OUT="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_plotlist"
os.makedirs(OUT, exist_ok=True)
import re
writers={}; files={}; counts=collections.Counter()
with open(SRC, newline="") as f:
    r=csv.DictReader(f)
    for row in r:
        var=(row.get("VARIANT") or "").strip().lower()
        if not var: continue
        cn=(row.get("STAND_CN") or "").strip()
        if not cn: continue
        grp=row.get("GROUPS") or ""
        m=re.search(r"State=(\d+)", grp); state=m.group(1) if m else ""
        iy=(row.get("INV_YEAR") or "").strip()
        lat=(row.get("LATITUDE") or "").strip(); lon=(row.get("LONGITUDE") or "").strip()
        if var not in writers:
            fh=open(os.path.join(OUT,f"plotlist_{var}.csv"),"w",newline=""); files[var]=fh
            w=csv.writer(fh); w.writerow(["PLOT","FIRST_PLTCN","STATECD","FIRST_INVYR","PUB_LAT","PUB_LONG","VARIANT"])
            writers[var]=w
        counts[var]+=1
        writers[var].writerow([counts[var], cn, state, iy, lat, lon, var.upper()])
for fh in files.values(): fh.close()
tot=sum(counts.values())
print("variants:", dict(sorted(counts.items(), key=lambda x:-x[1])))
print("TOTAL plots:", tot)
