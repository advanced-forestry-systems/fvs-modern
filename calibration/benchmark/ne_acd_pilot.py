import os, sys
P = "/users/PUOM0008/crsfaaron/fvs-modern"
os.environ["FIA_DATA_DIR"] = "/fs/scratch/PUOM0008/crsfaaron/FIA"
os.environ["FVS_PROJECT_ROOT"] = P; os.environ["FVS_LIB_DIR"] = P+"/lib"; os.environ["FVS_CONFIG_DIR"] = P+"/config"
for p in [P, P+"/calibration/python", P+"/calibration", P+"/deployment/fvs2py"]: sys.path.insert(0, p)
import fia_stand_generator as G
from perseus_100yr_projection import run_fvs_projection
sc, dc = "Medium", "Medium"
stands = G.generate_real_stand("NE", sc, dc, n_plots=8, seed=1)
print(f"NE {sc}/{dc} stands: {len(stands)}")
M2HA=0.2296; TPHc=2.4710538; CMc=2.54; M3c=0.06997; FTc=0.3048
def g(r,k): 
    try: return float(r.get(k,0) or 0)
    except: return 0.0
def met(s):
    f=s.iloc[0]; l=s.iloc[-1]
    return (g(f,"BA")*M2HA, g(l,"BA")*M2HA, g(l,"Tpa")*TPHc, g(l,"QMD")*CMc, g(l,"TopHt")*FTc, g(l,"TCuFt")*M3c)
runs=[("NE",None,"real_NE "),("ACD",None,"real_ACD"),("NE","calibrated","unified ")]
agg={lab.strip():[] for _,_,lab in runs}
for i,(std,tree,cond) in enumerate(stands):
    sid = str(tree["stand_id"].iloc[0])
    print(f"\n== {sid} ({len(tree)} trees) ==")
    for variant,cfg,label in runs:
        try:
            r = run_fvs_projection(std, tree, sid, variant, config_version=cfg, num_cycles=1, cycle_length=10)
            s = r.get("summary")
            if s is not None and len(s):
                ba0,ba,tph,qmd,tht,vol = met(s)
                print(f"  {label}: BA {ba0:5.1f}->{ba:5.1f} m2/ha | TPH {tph:5.0f} | QMD {qmd:4.1f}cm | TopHt {tht:4.1f}m | Vol {vol:5.1f} m3/ha")
                agg[label.strip()].append((ba0,ba,tph,qmd,vol))
            else: print(f"  {label}: empty (exit {r.get('exit_code')})")
        except Exception as e: print(f"  {label}: FAIL {repr(e)[:110]}")
print("\n== MEAN projected (10yr) across stands ==")
for k,v in agg.items():
    if v:
        import statistics as st
        ba=[x[1] for x in v]; tph=[x[2] for x in v]; vol=[x[4] for x in v]
        print(f"  {k:8s}: meanBA {st.mean(ba):5.1f} | meanTPH {st.mean(tph):5.0f} | meanVol {st.mean(vol):5.1f}  (n={len(v)})")
print("DONE_NE_ACD_PILOT")
