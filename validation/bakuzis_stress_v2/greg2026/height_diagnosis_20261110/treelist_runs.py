## Re-run the Bakuzis v2 scenarios with tree lists, reproducing the harness stands exactly (same rng order).
import os, sys, numpy as np, pandas as pd, argparse, importlib.util
ap=argparse.ArgumentParser(); ap.add_argument('--wt'); ap.add_argument('--lib'); ap.add_argument('--tag'); ap.add_argument('--out')
ap.add_argument('--groups', default='Northern-Hardwood,Oak-Pine'); ap.add_argument('--sites', default='Medium'); ap.add_argument('--dens', default='Medium')
a=ap.parse_args()
spec=importlib.util.spec_from_file_location('eb', os.path.join(a.wt,'validation/bakuzis_stress_v2/engine_bakuzis_v2.py')); eb=importlib.util.module_from_spec(spec); spec.loader.exec_module(eb)
eb.KEY_TMPL=eb.KEY_TMPL.replace('PROCESS\n','TREELIST           0         3\nPROCESS\n')
sys.path.insert(0, os.path.join(a.wt,'deployment','fvs2py'))
from fvs2py import FVS
import ctypes as _ct, fvs2py._core as _core
_libdl=_ct.CDLL('libdl.so.2')
def _unload(dll):
    f=_libdl.dlclose; f.argtypes=(_ct.c_void_p,); f.restype=_ct.c_int; f(dll._handle)
_core.unload_dll=_unload
cfg=os.path.join(a.wt,'config'); os.makedirs(a.out, exist_ok=True)
rng=np.random.default_rng(42); n=0
G=a.groups.split(','); S=a.sites.split(','); D=a.dens.split(',')
for spk,spp in eb.SPP.items():
  for sk,si in eb.SITE.items():
    for dk,(ba,tpa) in eb.DENS.items():
      n+=1; sid='BK%02d'%n
      if spk in G and sk in S and dk in D:
        rd=os.path.join(a.out,'run_%s_%s'%(a.tag,sid))
        s,sir=eb.run_one(FVS,a.lib,cfg,rd,sid,spk[:12],spp,ba,tpa,si,eb.SITE_SPP[spk],rng)
        print('RAN',sid,spk,sk,dk,sir,flush=True)
      else:
        eb.make_stand(spp,ba,tpa,si,rng)   # advance rng identically
print('DONE')
