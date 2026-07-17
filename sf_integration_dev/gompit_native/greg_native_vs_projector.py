import sys, csv
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs-modern-greg-consumer")
sys.path.insert(0, "/fs/scratch/PUOM0008/crsfaaron/fvs-modern-greg-consumer/sf_integration_dev")
from config.config_loader import FvsConfigLoader
from conus_greg_projector import GregProjector

L = FvsConfigLoader("ne", version="conus_greg")
gp = GregProjector(L, emt=-15.0, td=25.0)

# identical fixed state to test_gregdghg.f90
dbh,cr,ht,bal,ccfl,cch,elev,td,emt = 8.0,0.5,50.0,80.0,120.0,0.4,1500.0,25.0,-15.0

def load(path):
    d={}
    with open(path) as f:
        next(f)
        for line in f:
            s=line.strip().split(",")
            if len(s)<2: continue
            d[int(s[0])]=float(s[1])
    return d

fdg=load("/fs/scratch/PUOM0008/crsfaaron/fvs-modern-greg-hook/sf_integration_dev/gompit_native/dg_fortran.csv")
fhg=load("/fs/scratch/PUOM0008/crsfaaron/fvs-modern-greg-hook/sf_integration_dev/gompit_native/hg_fortran.csv")

dg_maxad=0.0; dg_n=0; dg_arg=None
for spcd in sorted(set(fdg)&gp.dg_set):
    pv=gp.dg_annual(spcd, dbh, cr, ht, bal, elev, emt=emt)
    ad=abs(pv-fdg[spcd]); dg_n+=1
    if ad>dg_maxad: dg_maxad=ad; dg_arg=spcd

hg_maxad=0.0; hg_n=0; hg_arg=None
for spcd in sorted(set(fhg)&gp.hg_set):
    pv=gp.hg_annual(spcd, ht, cr, ccfl, cch, elev, td=td, emt=emt)
    ad=abs(pv-fhg[spcd]); hg_n+=1
    if ad>hg_maxad: hg_maxad=ad; hg_arg=spcd

print("Native-Fortran vs Python-projector, annual increment, fixed state")
print("DG: n=%d  max|native-projector|=%.3e (spcd %s)" % (dg_n, dg_maxad, dg_arg))
print("HG: n=%d  max|native-projector|=%.3e (spcd %s)" % (hg_n, hg_maxad, hg_arg))
