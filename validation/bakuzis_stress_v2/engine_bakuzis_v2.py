#!/usr/bin/env python3
"""engine_bakuzis_v2.py  (2026-10-27)

Engine Bakuzis harness, second version. Derived from validation/bakuzis_engine_verdict/engine_bakuzis.py
(fvs-modern 895f5aa) with three changes that the 2026-10-27 review found necessary:

  1. SITE LEVER. The v1 keyfile wrote a SETSITE line. SETSITE is not an FVS keyword; the engine ignored it
     and every site class ran on the variant default, which is why the v1 verdict recorded
     "SETSITE did not vary SITEAR". This version writes SITECODE (site species, site index), the keyword the
     fvs-calibration synthetic gate used successfully on NE, CS, LS, WC in July 2026, and parses the
     "SITE INDEX" value back out of the .out file as a gate: a scenario whose reported site index is not the
     designed one is marked SITE_FAIL and still written, never silently passed.
  2. SEEDING GATE. Year-0 TPA and basal area from the engine summary are compared with the design cell.
     A cell more than 20 percent off its label is marked SEED_FAIL. The Paper 1 FIA-stand ensemble
     (calibration/output/bakuzis) started every cell at 1.1 to 2.4 m2/ha whatever its label; this gate
     makes that failure impossible to miss in any later run.
  3. TWO LIBRARIES. --lib may be given twice (or --lib-a/--lib-b) so the GOMPIT and EXP-HAZ survival
     transforms run on identical stands in one job; each gets its own tag in the output.

Output: <out>/engine_bakuzis_100yr_<tag>.csv in the layout bakuzis_stress_v2.R --source=engine reads
(year, age, tpa, tcuft, atba, attopht, species_group, site_class, density_class, scenario, plus
seed_ok, site_ok, site_index_reported, tag), and <out>/seed_site_gate_<tag>.csv.

Usage (Cardinal, inside the wt-engine worktree with fvs2py built):
  python3 engine_bakuzis_v2.py --wt /fs/scratch/PUOM0008/crsfaaron/wt-engine \
      --lib $W/lib_gompit/FVSne.so --tag gompit --lib $W/lib_exp/FVSne.so --tag exp --out $W/engine_bakuzis_v2_out
"""
import os, sys, re, faulthandler, argparse
import numpy as np, pandas as pd
faulthandler.enable()

SPP = {
 'Spruce-Fir':       [('BF',0.40),('RS',0.35),('BS',0.15),('YB',0.10)],
 'Northern-Hardwood':[('SM',0.35),('YB',0.25),('AB',0.20),('WA',0.20)],
 'Pine':             [('WP',0.60),('RM',0.15),('BF',0.15),('EH',0.10)],
 'Oak-Pine':         [('RO',0.35),('WP',0.25),('RM',0.20),('WO',0.20)],
}
SITE_SPP  = {'Spruce-Fir':'BF','Northern-Hardwood':'SM','Pine':'WP','Oak-Pine':'RO'}
SITE_SPCD = {'Spruce-Fir':12,'Northern-Hardwood':318,'Pine':129,'Oak-Pine':833}   # FIA codes, used with --site-species-numeric
SITE = {'Low':45,'Medium':60,'High':75}
DENS = {'Low':(60,150),'Medium':(120,250),'High':(180,400)}   # (BA ft2/ac, TPA) design cell

def rec(tid,prob,sp,dbh,ht):
    b=[' ']*62
    def put(s,c0):
        for i,ch in enumerate(s): b[c0-1+i]=ch
    put('%7d'%tid,1); put('%6.0f'%prob,8); put('%1d'%1,14)
    put('%-3s'%sp,15); put('%4.1f'%dbh,18); put('%3.0f'%ht,25)
    return ''.join(b)

def make_stand(spp, ba, tpa, si, rng):
    """Lognormal diameter list per species whose summed TPA and BA hit the design cell exactly.
    v1 drew diameters around the species QMD and let BA fall where it fell; here diameters are rescaled
    so the seeding gate tests the ENGINE, not the generator."""
    rows=[]   # (species, tpa_each, dbh)
    for sp,share in spp:
        sp_tpa=tpa*share; sp_ba=ba*share
        if sp_tpa<=0: continue
        n=max(5,int(sp_tpa/6.0))
        qmd=np.sqrt(sp_ba*576/(np.pi*sp_tpa))
        dias=np.clip(rng.lognormal(np.log(max(qmd,1.5)),0.35,n),1.0,30.0)
        for d in dias: rows.append([sp,sp_tpa/n,float(d)])
    each=np.array([r[1] for r in rows]); dbh=np.array([r[2] for r in rows])
    ba_now=np.sum(0.005454154*dbh**2*each); dbh=np.clip(dbh*np.sqrt(ba/ba_now),1.0,40.0)
    recs=[]
    for tid,(r,d) in enumerate(zip(rows,dbh),start=1):
        ht=4.5+(si*1.1)*(1-np.exp(-0.04*d))**1.2
        ht=max(10,min(ht,120))
        recs.append(rec(tid,r[1],r[0],float(d),float(ht)))
    return recs

KEY_TMPL='''STDIDENT
{sid}  {label}
STDINFO          922                   1
DESIGN            -1         1
INVYEAR       2000.0
SITECODE  {sitesp:>10}{si:>10.1f}
NUMCYCLE        20.0
TIMEINT            0         5
TREEDATA
PROCESS
STOP
'''

def parse_site_index(outfile):
    """Return the site index the engine reports in the .out, or None if no .out / no match."""
    if not os.path.exists(outfile): return None
    txt=open(outfile,errors='replace').read()
    m=re.search(r'SITE CODE FOR LARGE AND SMALL TREES=\s*([0-9]+\.?[0-9]*)',txt)
    return float(m.group(1)) if m else None

def run_one(FVS, lib, cfg, rundir, sid, label, spp, ba, tpa, si, sitesp, rng):
    os.makedirs(rundir, exist_ok=True)
    for f in os.listdir(rundir):
        try: os.remove(os.path.join(rundir,f))
        except Exception: pass
    recs=make_stand(spp, ba, tpa, si, rng)
    base=os.path.join(rundir,'bk')
    open(base+'.key','w').write(KEY_TMPL.format(sid=sid,label=label,sitesp=sitesp,si=float(si)))
    open(base+'.tre','w').write('\n'.join(recs)+'\n')
    os.chdir(rundir)
    fvs=FVS(lib_path=lib, config_version=None, config_dir=cfg)
    fvs.load_keyfile(base+'.key'); fvs.run()
    s=fvs.summary
    si_rep=parse_site_index(base+'.out')
    return (s.copy() if s is not None else None), si_rep

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--wt', required=True, help='worktree holding deployment/fvs2py and config/')
    ap.add_argument('--lib', action='append', required=True, help='shared library; repeat for an A/B')
    ap.add_argument('--tag', action='append', default=[], help='one tag per --lib, same order')
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--site-species-numeric', action='store_true', help='write the FIA code instead of the alpha code on SITECODE')
    ap.add_argument('--seed-tol', type=float, default=0.20)
    a=ap.parse_args()
    tags=a.tag if a.tag else ['lib%d'%i for i in range(len(a.lib))]
    if len(tags)!=len(a.lib): sys.exit('one --tag per --lib')
    sys.path.insert(0, os.path.join(a.wt,'deployment','fvs2py'))
    from fvs2py import FVS
    # glibc < 2.34 keeps dlclose in libdl.so.2, which the FVS .so does not list as NEEDED, so fvs2py's
    # unload (dll.dlclose looked up on the FVS handle) fails on the second scenario. Route it through libdl.
    import ctypes as _ct, fvs2py._core as _core
    try:
        _libdl=_ct.CDLL('libdl.so.2')
        def _unload(dll):
            f=_libdl.dlclose; f.argtypes=(_ct.c_void_p,); f.restype=_ct.c_int; f(dll._handle)
        _core.unload_dll=_unload
    except OSError:
        pass
    cfg=os.path.join(a.wt,'config'); os.makedirs(a.out, exist_ok=True)
    for lib,tag in zip(a.lib,tags):
        rng=np.random.default_rng(a.seed)          # identical stands for every lib
        rows=[]; gate=[]; n=0
        for spk,spp in SPP.items():
            for sk,si in SITE.items():
                for dk,(ba,tpa) in DENS.items():
                    n+=1; sid='BK%02d'%n
                    sitesp = str(SITE_SPCD[spk]) if a.site_species_numeric else SITE_SPP[spk]
                    rundir=os.path.join(a.out,'run_%s_%s'%(tag,sid))
                    try:
                        s,si_rep=run_one(FVS, lib, cfg, rundir, sid, spk[:8], spp, ba, tpa, si, sitesp, rng)
                    except Exception as e:
                        print('FAIL %s %s %s/%s/%s : %s'%(tag,sid,spk,sk,dk,e),flush=True); gate.append(dict(tag=tag,scenario=n,status='RUN_FAIL',err=str(e))); continue
                    if s is None or s.empty:
                        print('EMPTY %s %s'%(tag,sid),flush=True); gate.append(dict(tag=tag,scenario=n,status='EMPTY')); continue
                    tpa0=float(s['tpa'].iloc[0]); ba0=float(s['atba'].iloc[0])
                    seed_ok = abs(tpa0-tpa)/tpa <= a.seed_tol and abs(ba0-ba)/ba <= a.seed_tol
                    site_ok = (si_rep is not None) and abs(si_rep-si) <= 1.0
                    s['species_group']=spk; s['site_class']=sk; s['density_class']=dk; s['scenario']=n
                    s['seed_ok']=seed_ok; s['site_ok']=site_ok; s['site_index_reported']=si_rep; s['tag']=tag
                    rows.append(s)
                    gate.append(dict(tag=tag,scenario=n,species_group=spk,site_class=sk,density_class=dk,design_tpa=tpa,design_ba=ba,
                                     tpa0=tpa0,ba0=ba0,seed_ok=seed_ok,design_si=si,si_reported=si_rep,site_ok=site_ok,
                                     status='OK' if (seed_ok and site_ok) else ('SEED_FAIL ' if not seed_ok else '')+('SITE_FAIL' if not site_ok else '')))
                    print('%s %s %s %s/%s/%s tpa0=%.0f/%d ba0=%.0f/%d si=%s/%d'%('OK  ' if seed_ok and site_ok else 'GATE',tag,sid,spk,sk,dk,tpa0,tpa,ba0,ba,si_rep,si),flush=True)
        g=pd.DataFrame(gate); g.to_csv(os.path.join(a.out,'seed_site_gate_%s.csv'%tag),index=False)
        if rows:
            df=pd.concat(rows,ignore_index=True); out=os.path.join(a.out,'engine_bakuzis_100yr_%s.csv'%tag)
            df.to_csv(out,index=False); print('WROTE',out,'rows',len(df),flush=True)
        ok=int((g['status']=='OK').sum()) if 'status' in g else 0
        print('GATE SUMMARY %s: %d of %d scenarios pass seeding and site gates; %d SEED_FAIL, %d SITE_FAIL'%(
            tag, ok, len(g), int(g.get('seed_ok',pd.Series(dtype=bool)).eq(False).sum()), int(g.get('site_ok',pd.Series(dtype=bool)).eq(False).sum())),flush=True)
    print('BAKUZIS_ENGINE_V2_DONE',flush=True)

if __name__=='__main__':
    main()
