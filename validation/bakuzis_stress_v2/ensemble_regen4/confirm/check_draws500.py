"""D4: per-species draw median vs MAP for the ensemble's 100 draws and for all 500 stored draws (NE, ACD)."""
import sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
F='/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030'; X='/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_20261030'
sys.path[:0]=[F]
from config.uncertainty import UncertaintyEngine
from config.config_loader import FvsConfigLoader
rows=[]
for v in ['ne','acd']:
    e=UncertaintyEngine(v,config_dir=F+'/config',seed=42)
    idx100=set(int(i) for i in np.random.default_rng(42).choice(e.n_draws,size=100,replace=False))
    L=FvsConfigLoader(v,version='calibrated',config_dir=F+'/config'); cm=L.config['calibration_multipliers']
    mapb=np.sqrt(np.clip(np.array(cm['dds_multiplier'],float),.01,100)); mapm=np.array(cm['mort_multiplier'],float)
    sdi=FvsConfigLoader._find_sdi_param(L.config['categories'])
    smap=np.array([np.nan if (sdi is None or i>=len(sdi) or isinstance(sdi[i],str)) else float(sdi[i]) for i in range(108)])
    B=[];M=[];S=[];I=[]
    for di in range(e.n_draws):
        m=e.compute_multipliers_for_draw(e.get_draw(di)); B.append(m['growth_mult']); M.append(m['mort_mult']); S.append(m.get('sdi_values',np.full(108,np.nan))[:108]); I.append(di in idx100)
    B,M,S,I=map(np.array,(B,M,S,I))
    for i in range(108):
        r=dict(variant=v.upper(),fvs_sp=i+1,baimult_map=mapb[i],mortmult_map=mapm[i],sdimax_map=smap[i])
        for lab,sel in (('d100',I),('d500',np.ones_like(I))):
            r[f'baimult_med_{lab}']=np.median(B[sel,i]); r[f'mortmult_med_{lab}']=np.median(M[sel,i])
            r[f'sdimax_med_{lab}']=np.nanmedian(S[sel,i]) if np.isfinite(S[sel,i]).any() else np.nan
        r['baimult_q025_d100']=np.quantile(B[I,i],.025); r['baimult_q975_d100']=np.quantile(B[I,i],.975)
        r['mortmult_q025_d100']=np.quantile(M[I,i],.025); r['mortmult_q975_d100']=np.quantile(M[I,i],.975)
        rows.append(r)
d=pd.DataFrame(rows); d.to_csv(X+'/d4_draw_median_vs_map.csv',index=False)
for v,g in d.groupby('variant'):
    for fam in ('baimult','mortmult','sdimax'):
        cal=g[(g[f'{fam}_map'].notna())&((g[f'{fam}_map']-1).abs()>0.01 if fam!='sdimax' else True)]
        if cal.empty: print(v,fam,'no calibrated species'); continue
        for lab in ('d100','d500'):
            r=cal[f'{fam}_med_{lab}']/cal[f'{fam}_map']
            print(v,fam,lab,'n_sp',len(cal),'med/MAP median %.4f min %.4f max %.4f within5pct %d'%(r.median(),r.min(),r.max(),(abs(r-1)<=.05).sum()))
    print(v,'BAIMULT MAP range %.3f-%.3f, draw-median(d100) range %.3f-%.3f, species with any draw BAIMULT != 1: %d'%(g.baimult_map.min(),g.baimult_map.max(),g.baimult_med_d100.min(),g.baimult_med_d100.max(),((g.baimult_q025_d100!=1)|(g.baimult_q975_d100!=1)).sum()))
print('D4_DONE')
