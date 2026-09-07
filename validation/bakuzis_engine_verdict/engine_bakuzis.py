#!/usr/bin/env python3
# Engine Bakuzis harness: 36 scenarios (4 spp-grp x 3 site x 3 density),
# 20 cycles x 5 yr = 100 yr, GOMPIT lib, inline-inventory (.tre) path.
# Reads summaries from ENGINE MEMORY (fvs.summary), no SQLite/DSNIN.
import os, sys, faulthandler, argparse
import numpy as np, pandas as pd
faulthandler.enable()
WT='/fs/scratch/PUOM0008/crsfaaron/wt-engine'
sys.path.insert(0, os.path.join(WT,'deployment','fvs2py'))
from fvs2py import FVS

CFG=os.path.join(WT,'config')

# NE FVS alpha codes for our species
SPP = {
 'Spruce-Fir':      [('BF',0.40),('RS',0.35),('BS',0.15),('YB',0.10)],
 'Northern-Hardwood':[('SM',0.35),('YB',0.25),('AB',0.20),('WA',0.20)],
 'Pine':            [('WP',0.60),('RM',0.15),('BF',0.15),('EH',0.10)],
 'Oak-Pine':        [('RO',0.35),('WP',0.25),('RM',0.20),('WO',0.20)],
}
SITE_SPP = {'Spruce-Fir':'BF','Northern-Hardwood':'SM','Pine':'WP','Oak-Pine':'RO'}
# FVS NE numeric spcd for SETSITE site species
SITE_SPCD = {'Spruce-Fir':12,'Northern-Hardwood':318,'Pine':129,'Oak-Pine':833}
SITE = {'Low':45,'Medium':60,'High':75}
ELEV= {'Low':1800,'Medium':800,'High':400}
DENS= {'Low':(60,150),'Medium':(120,250),'High':(180,400)}

def rec(tid,prob,sp,dbh,ht):
    b=[' ']*62
    def put(s,c0):
        for i,ch in enumerate(s): b[c0-1+i]=ch
    put('%7d'%tid,1); put('%6.0f'%prob,8); put('%1d'%1,14)
    put('%-3s'%sp,15); put('%4.1f'%dbh,18); put('%3.0f'%ht,25)
    return ''.join(b)

def make_stand(spp, ba, tpa, si, rng):
    recs=[]; tid=1
    for sp,share in spp:
        sp_tpa=tpa*share; sp_ba=ba*share
        if sp_tpa<=0: continue
        n=max(5,int(sp_tpa/6.0))
        qmd=np.sqrt(sp_ba*576/(np.pi*sp_tpa)) if sp_tpa>0 else 6.0
        dias=np.clip(rng.lognormal(np.log(max(qmd,1.5)),0.35,n),1.0,30.0)
        each=sp_tpa/n
        for d in dias:
            ht=4.5+(si*1.1)*(1-np.exp(-0.04*d))**1.2
            ht=max(10,min(ht,120))
            recs.append(rec(tid,each,sp,float(d),float(ht))); tid+=1
    return recs

KEY_TMPL='''STDIDENT
{sid}  {label}
STDINFO          922                   1
DESIGN            -1         1
INVYEAR       2000.0
SETSITE            0        0.{spcd:9.0f}.{si:9.1f}
NUMCYCLE        20.0
TIMEINT            0         5
TREEDATA
PROCESS
STOP
'''

def run_one(lib, rundir, sid, label, spp, ba, tpa, si, spcd, rng):
    os.makedirs(rundir, exist_ok=True)
    for f in os.listdir(rundir):
        try: os.remove(os.path.join(rundir,f))
        except: pass
    recs=make_stand(spp, ba, tpa, si, rng)
    base=os.path.join(rundir,'bk')
    open(base+'.key','w').write(KEY_TMPL.format(sid=sid,label=label,spcd=float(spcd),si=float(si)))
    open(base+'.tre','w').write('\n'.join(recs)+'\n')
    os.chdir(rundir)
    fvs=FVS(lib_path=lib, config_version=None, config_dir=CFG)
    fvs.load_keyfile(base+'.key'); fvs.run()
    s=fvs.summary
    return s.copy() if s is not None else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--lib', required=True)
    ap.add_argument('--out', default=os.path.join(WT,'engine_bakuzis_out'))
    a=ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng=np.random.default_rng(42)
    rows=[]; n=0
    for spk,spp in SPP.items():
        for sk,si in SITE.items():
            for dk,(ba,tpa) in DENS.items():
                n+=1; sid='BK%02d'%n
                rundir=os.path.join(a.out,'run_%s'%sid)
                try:
                    s=run_one(a.lib, rundir, sid, spk[:8], spp, ba, tpa, si, SITE_SPCD[spk], rng)
                    if s is not None and not s.empty:
                        s['species_group']=spk; s['site_class']=sk; s['density_class']=dk
                        s['scenario']=n; rows.append(s)
                        print('OK %s %s/%s/%s rows=%d maxyr=%s tpa0=%s'%(sid,spk,sk,dk,len(s),int(s['year'].max()),float(s['tpa'].iloc[0])),flush=True)
                    else:
                        print('EMPTY %s %s/%s/%s'%(sid,spk,sk,dk),flush=True)
                except Exception as e:
                    print('FAIL %s %s/%s/%s : %s'%(sid,spk,sk,dk,e),flush=True)
    if rows:
        df=pd.concat(rows,ignore_index=True)
        out=os.path.join(a.out,'engine_bakuzis_100yr.csv')
        df.to_csv(out,index=False)
        print('WROTE',out,'rows',len(df),flush=True)
    else:
        print('NO_ROWS',flush=True)
    print('BAKUZIS_ENGINE_DONE',flush=True)

if __name__=='__main__':
    main()
