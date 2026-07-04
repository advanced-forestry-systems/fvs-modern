import sys, numpy as np, pandas as pd
sys.path.insert(0,"/fs/scratch/PUOM0008/crsfaaron/akwork/pylib")
import shapefile
from matplotlib.path import Path
SHP="/fs/scratch/PUOM0008/crsfaaron/akwork/eco_shp/NA_Eco_L3_WGS84.shp"
r=shapefile.Reader(SHP)
flds=[f[0] for f in r.fields[1:]]; print("fields:",flds,flush=True)
def fi(pat):
    for i,f in enumerate(flds):
        if pat.lower() in f.lower(): return i
    return None
iL1=fi("L1NAME") or fi("L1_NAME"); iL3=fi("L3NAME"); iL3c=fi("L3CODE") or fi("L3_KEY")
print("idx L1",iL1,"L3",iL3,"L3code",iL3c,flush=True)
polys=[]
for sr in r.iterShapeRecords():
    sh=sr.shape; rec=sr.record
    if not sh.points: continue
    parts=list(sh.parts)+[len(sh.points)]
    paths=[Path(sh.points[parts[i]:parts[i+1]]) for i in range(len(parts)-1) if parts[i+1]-parts[i]>=3]
    bb=sh.bbox  # xmin,ymin,xmax,ymax
    polys.append((bb,paths,rec[iL1] if iL1 is not None else None,rec[iL3] if iL3 is not None else None,rec[iL3c] if iL3c is not None else None))
print("polygons:",len(polys),flush=True)
s=pd.read_csv("/users/PUOM0008/crsfaaron/magplot/magp_sites.csv")
s=s[(s["province"]=="BC")&np.isfinite(s["latitude"])&np.isfinite(s["longitude"])].copy()
out=[]
for _,row in s.iterrows():
    x,y=row["longitude"],row["latitude"]; hit=(None,None,None)
    for bb,paths,l1,l3,l3c in polys:
        if x<bb[0] or x>bb[2] or y<bb[1] or y>bb[3]: continue
        if any(p.contains_point((x,y)) for p in paths): hit=(l1,l3,l3c); break
    out.append((row["magp_site_id"],hit[0],hit[1],hit[2]))
o=pd.DataFrame(out,columns=["magp_site_id","ecoL1","ecoL3","ecoL3code"])
o.to_csv("/fs/scratch/PUOM0008/crsfaaron/akwork/bc_site_ecoregion.csv",index=False)
print("joined",len(o),"matched",o["ecoL3"].notna().sum(),flush=True)
print("=== by L1 ==="); print(o[o.ecoL1.notna()].groupby("ecoL1").size().sort_values(ascending=False).to_string())
print("=== top L3 ==="); print(o[o.ecoL3.notna()].groupby(["ecoL3code","ecoL3"]).size().sort_values(ascending=False).head(18).to_string())
print("DONE")
