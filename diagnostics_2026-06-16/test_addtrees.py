#!/usr/bin/env python3
"""Route A mechanism test: load a stand in memory via fvsAddTrees with NO keyword-file database, then
project. If this produces a summary without the extree segfault, the in-process injection path is open.
Synthetic stand first (isolates the mechanism from FIA species mapping)."""
import os, sys, tempfile, faulthandler
faulthandler.enable()
P="/users/PUOM0008/crsfaaron/fvs-modern"
for p in [P+"/deployment/fvs2py", P]: sys.path.insert(0, p)
import numpy as np
try:
    from fvs2py import FVS
except Exception:
    from fvs2py._base import FVS
LIB=os.environ.get("LIB", P+"/lib-test/FVSne.so")
def pr(*a): print(*a, flush=True)

# database-free keyword file: stand setup + summary to the API, trees come from fvsAddTrees
KEY = """STDIDENT
ADDTREE_TEST
STDINFO        1        11       470         5        50
INVYEAR     2000
NOTREES
TIMEINT        0        10
NUMCYCLE       1
PROCESS
STOP
"""
tmp=tempfile.mkdtemp(); kp=os.path.join(tmp,"t.key"); open(kp,"w").write(KEY)
pr("loading", LIB)
fvs=FVS(LIB); fvs.load_keyfile(kp)
pr("keyfile loaded; itrncd", fvs.itrncd)
fvs.run(stop_point_code=7)
pr("after run(7): restart", fvs.restart_code, "ntrees", fvs.dims["ntrees"])

# synthetic stand: 50 trees, a spread of sizes, species index 1, crown ratio 40, plot 1
n=50
rng=np.random.default_rng(5)
dbh=rng.uniform(5,20,n); ht=8.0+2.5*dbh; cr=np.full(n,40.0); sp=np.ones(n); plot=np.ones(n); tpa=np.full(n,5.0)
added=fvs.add_trees(dbh, sp, ht, cr, plot, tpa)
pr("added", added, "trees; ntrees now", fvs.dims["ntrees"])
tbl=fvs.tree_table(("species","dbh","ht"))
pr("in-engine after add: n", len(tbl), " mean dbh", round(float(tbl["dbh"].mean()),2))

# project to completion
k=0
while fvs.itrncd==0 and k<8:
    fvs.run(); k+=1
    if fvs.restart_code==100 or fvs.itrncd!=0: break
s=fvs.summary
bacol=next((c for c in (s.columns if s is not None else []) if str(c).lower() in ("atba","ba","baa")),None)
volcol=next((c for c in (s.columns if s is not None else []) if str(c).lower() in ("mcuft","tcuft")),None)
pr("summary cols", list(s.columns) if s is not None else None)
if s is not None and len(s):
    pr("summary rows:"); pr(s[[c for c in ("year","age","tpa",bacol,volcol,"attopht") if c in s.columns]].to_string(index=False))
ba=float(s[bacol].iloc[-1]) if (s is not None and len(s) and bacol) else float("nan")
pr("final stand BA (atba):", round(ba,2))
pr("ROUTE_A_INPROCESS_PROJECTION_WORKS" if ba>0 else "RAN_BUT_BA_ZERO")
pr("DONE_ADDTREES")
