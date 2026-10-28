#!/usr/bin/env python3
"""diag_fia_seeding.py (2026-10-28): why do FIA-seeded Bakuzis stands not reproduce BALIVE?
For NE conditions in each SI x BA bin, recompute plot basal area from the TREE table the way
fia_stand_generator.load_fia_trees does (live trees on PLT_CN, TPA_UNADJ) and compare with BALIVE.
Also tests three candidate fixes: CONDID match, DIA >= 1.0 filter, and PLT_CN dtype/precision.
Prints only aggregate ratios and counts; no plot identifiers or coordinates."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.environ['W'], 'calibration', 'python'))
import fia_stand_generator as g
fia = g.Path(os.environ.get('FIA_DATA_DIR', os.path.expanduser('~/fia_data')))
cond = g.load_fia_conditions('ne', fia)
print('NE COND rows (forested):', len(cond), 'PLT_CN dtype', cond['PLT_CN'].dtype, 'INVYR range', cond['INVYR'].min(), cond['INVYR'].max())
# a sample across the nine cells
rows = []
for sc in g.SI_BINS:
    for dc in g.BA_BINS:
        s = g.stratified_sample_conditions(cond, sc, dc, n_plots=20, seed=7)
        s = s.assign(site_class=sc, density_class=dc); rows.append(s)
samp = pd.concat(rows, ignore_index=True)
print('sampled conditions:', len(samp))
plt_cns = samp['PLT_CN'].astype(int).tolist()
trees = g.load_fia_trees('ne', plt_cns, fia)
print('trees loaded:', len(trees), 'PLT_CN dtype', trees['PLT_CN'].dtype, 'columns', [c for c in ['CONDID','TPA_UNADJ','DIA','STATUSCD'] if c in trees.columns])
trees['ba_rec'] = 0.005454154 * trees['DIA'].fillna(0) ** 2 * trees['TPA_UNADJ'].fillna(0)
out = []
for _, r in samp.iterrows():
    t = trees.loc[trees['PLT_CN'] == r['PLT_CN']]
    out.append(dict(site=r['site_class'], dens=r['density_class'], balive=r['BALIVE'], n_trees=len(t),
                    tpa=t['TPA_UNADJ'].sum(), ba=t['ba_rec'].sum(), ba_ge5=t.loc[t['DIA'] >= 5, 'ba_rec'].sum(),
                    n_micro=(t['TPA_UNADJ'] > 50).sum(), n_dia_na=t['DIA'].isna().sum()))
o = pd.DataFrame(out); o['ratio'] = o['ba'] / o['balive']
pd.set_option('display.width', 200)
print('\nratio of recomputed plot BA to BALIVE, by cell:')
print(o.groupby(['site', 'dens'])['ratio'].describe()[['count', 'mean', '25%', '50%', '75%']].round(2))
print('\nconditions with no trees found:', int((o['n_trees'] == 0).sum()), 'of', len(o))
print('median n_trees', o['n_trees'].median(), 'median TPA', round(o['tpa'].median()), 'median BA', round(o['ba'].median(), 1), 'median BALIVE', round(o['balive'].median(), 1))
print('share of BA from DIA>=5:', round((o['ba_ge5'].sum() / max(o['ba'].sum(), 1e-9)), 3))
# candidate cause: does the TREE csv cover the sampled INVYRs? compare by year
samp2 = samp.assign(found=[int(x) in set(trees['PLT_CN'].astype(int)) for x in samp['PLT_CN']])
print('\nfound-trees rate by INVYR:'); print(samp2.groupby('INVYR')['found'].agg(['mean', 'count']).round(2).T)
# multi-condition plots: how many sampled PLT_CNs have >1 forested condition?
nc = cond.groupby('PLT_CN').size(); samp2['n_cond'] = samp2['PLT_CN'].map(nc)
print('\nconditions per sampled plot:', samp2['n_cond'].value_counts().to_dict())
o2 = o.assign(n_cond=samp2['n_cond'].values)
print('ratio by n_cond:'); print(o2.groupby('n_cond')['ratio'].median().round(2).to_dict())
