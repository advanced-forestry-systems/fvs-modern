#!/usr/bin/env python3
"""
Volume scorecard for corrected (fix) variant runs: sn, wc, em, pn.
Adapted from vol_scorecard.py — same FIA join logic, no greg_dg config.
EM corrected states are MT (30) and SD (46), not IA/MN/WI as before.
"""

import pandas as pd
import numpy as np
import os
import sys

STATE_ABBR = {
    1:'AL',  5:'AR',  6:'CA',  8:'CO', 12:'FL', 13:'GA',
    16:'ID', 21:'KY', 22:'LA', 23:'ME', 24:'MD', 25:'MA',
    28:'MS', 30:'MT', 32:'NV', 33:'NH', 34:'NJ', 35:'NM',
    36:'NY', 37:'NC', 40:'OK', 41:'OR', 42:'PA', 45:'SC',
    46:'SD', 47:'TN', 48:'TX', 49:'UT', 50:'VT', 51:'VA',
    53:'WA', 54:'WV', 56:'WY'
}

FIA_DIR   = '/fs/scratch/PUOM0008/crsfaaron/FIA'
SC_DIR    = '/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard'

CONFIG_ORDER  = ["default","calibrated","hg_bayes","hg_organon","hg_conus_sp","hg_conus_clim"]
VARIANT_ORDER = ["sn","wc","em","pn"]

PLOT_COLS = ['CN','STATECD','COUNTYCD','PLOT','MEASYEAR']
TREE_COLS = ['PLT_CN','STATUSCD','VOLCFNET','VOLBFNET','TPA_UNADJ']

# ─────────────────────────────────────────────────────────────────────────────
def get_obs_volumes(sc_sub):
    """Given a slice of the scorecard (one or more states), return a DataFrame
    with columns [PLOT, MEASYEAR2, MCuFt_OBS, BdFt_OBS] by running the FIA
    two-step join for every unique state in sc_sub."""
    obs_rows = []
    states_needed = sorted(sc_sub['STATE'].unique())

    for state in states_needed:
        abbr = STATE_ABBR.get(state)
        if not abbr:
            print(f"  WARNING: no abbreviation for STATECD={state}", flush=True)
            continue

        plot_file = os.path.join(FIA_DIR, f'{abbr}_PLOT.csv')
        tree_file = os.path.join(FIA_DIR, f'{abbr}_TREE.csv')
        for f in [plot_file, tree_file]:
            if not os.path.exists(f):
                print(f"  WARNING: {f} not found, skipping state {state}", flush=True)
                continue

        state_sc = sc_sub[sc_sub['STATE'] == state][['PLOT','MEASYEAR2']].drop_duplicates()
        measyear1_cns = set(state_sc['PLOT'].astype(str))
        needed_myr2   = set(state_sc['MEASYEAR2'].astype(int))

        print(f"\n  State {state} ({abbr}): {len(measyear1_cns)} MEASYEAR1 CNs, "
              f"target years={sorted(needed_myr2)}", flush=True)

        # Step A: PLOT table
        try:
            plot_df = pd.read_csv(plot_file,
                                  usecols=lambda c: c in PLOT_COLS,
                                  low_memory=False)
        except Exception as exc:
            print(f"  ERROR reading {plot_file}: {exc}", flush=True)
            continue

        plot_df['CN']       = plot_df['CN'].astype(str)
        plot_df['MEASYEAR'] = pd.to_numeric(plot_df['MEASYEAR'], errors='coerce')
        plot_df['COUNTYCD'] = pd.to_numeric(plot_df['COUNTYCD'], errors='coerce').astype('Int64')
        plot_df['PLOT']     = pd.to_numeric(plot_df['PLOT'],     errors='coerce').astype('Int64')

        m1_rows    = plot_df[plot_df['CN'].isin(measyear1_cns)][['CN','COUNTYCD','PLOT']].copy()
        cn_to_phys = {row.CN: (row.COUNTYCD, row.PLOT) for _, row in m1_rows.iterrows()}
        print(f"  Physical IDs found: {len(cn_to_phys)}/{len(measyear1_cns)}", flush=True)
        if not cn_to_phys:
            continue

        needed_physplots = set(cn_to_phys.values())
        m2_cands = plot_df[plot_df['MEASYEAR'].isin(needed_myr2)].copy()
        m2_cands['phys'] = list(zip(m2_cands['COUNTYCD'], m2_cands['PLOT']))
        m2_cands = m2_cands[m2_cands['phys'].isin(needed_physplots)]
        phys_yr_to_cn2 = {
            (row.COUNTYCD, row.PLOT, int(row.MEASYEAR)): row.CN
            for _, row in m2_cands.iterrows()
            if pd.notna(row.MEASYEAR)
        }
        print(f"  (physplot, MEASYEAR2)->CN2 mappings: {len(phys_yr_to_cn2)}", flush=True)

        sc_to_cn2  = {}
        no_match   = 0
        for _, row in state_sc.iterrows():
            m1_cn = str(row['PLOT'])
            myr2  = int(row['MEASYEAR2'])
            phys  = cn_to_phys.get(m1_cn)
            if phys is None:
                no_match += 1; continue
            cn2 = phys_yr_to_cn2.get((phys[0], phys[1], myr2))
            if cn2 is None:
                no_match += 1; continue
            sc_to_cn2[(m1_cn, myr2)] = cn2
        print(f"  Resolved {len(sc_to_cn2)} plot-periods ({no_match} unresolved)", flush=True)
        if not sc_to_cn2:
            continue

        needed_cn2s = set(sc_to_cn2.values())

        # Step B: TREE table
        try:
            chunks = []
            for chunk in pd.read_csv(tree_file,
                                     usecols=lambda c: c in TREE_COLS,
                                     chunksize=400_000,
                                     low_memory=False):
                chunk['PLT_CN'] = chunk['PLT_CN'].astype(str)
                filt = chunk[chunk['PLT_CN'].isin(needed_cn2s)]
                if len(filt):
                    chunks.append(filt)
        except Exception as exc:
            print(f"  ERROR reading {tree_file}: {exc}", flush=True)
            continue

        if not chunks:
            print(f"  No tree records for state {state}", flush=True)
            continue

        trees = pd.concat(chunks, ignore_index=True)
        trees = trees[trees['STATUSCD'] == 1].copy()
        for col in ['VOLCFNET','VOLBFNET','TPA_UNADJ']:
            trees[col] = pd.to_numeric(trees[col], errors='coerce').fillna(0.0)
        trees['cf_contrib'] = trees['VOLCFNET'] * trees['TPA_UNADJ']
        trees['bf_contrib'] = trees['VOLBFNET'] * trees['TPA_UNADJ']

        agg = (trees.groupby('PLT_CN')
                    .agg(MCuFt_OBS=('cf_contrib','sum'),
                         BdFt_OBS =('bf_contrib','sum'))
                    .reset_index())
        agg['PLT_CN'] = agg['PLT_CN'].astype(str)
        cn2_to_sc = {v: k for k, v in sc_to_cn2.items()}
        agg['sc_PLOT']      = agg['PLT_CN'].map(lambda x: cn2_to_sc.get(x,(None,None))[0])
        agg['sc_MEASYEAR2'] = agg['PLT_CN'].map(lambda x: cn2_to_sc.get(x,(None,None))[1])
        agg = agg.dropna(subset=['sc_PLOT'])
        obs_rows.append(agg[['sc_PLOT','sc_MEASYEAR2','MCuFt_OBS','BdFt_OBS']])
        print(f"  -> {len(agg)} plot-period volume records", flush=True)

    if not obs_rows:
        return pd.DataFrame(columns=['PLOT','MEASYEAR2','MCuFt_OBS','BdFt_OBS'])

    obs = pd.concat(obs_rows, ignore_index=True)
    obs = obs.rename(columns={'sc_PLOT':'PLOT','sc_MEASYEAR2':'MEASYEAR2'})
    obs['PLOT']      = obs['PLOT'].astype(str)
    obs['MEASYEAR2'] = obs['MEASYEAR2'].astype(int)
    return obs

# ─────────────────────────────────────────────────────────────────────────────
def metrics(grp, pred_col, obs_col):
    df = grp[[pred_col, obs_col]].dropna()
    n  = len(df)
    if n == 0:
        return dict(n=0, obs_mean=np.nan, pred_mean=np.nan,
                    bias_pct=np.nan, RMSE=np.nan, MAE=np.nan)
    obs  = df[obs_col].values
    pred = df[pred_col].values
    obs_mean  = obs.mean()
    pred_mean = pred.mean()
    bias_pct  = 100*(pred_mean - obs_mean)/obs_mean if obs_mean else np.nan
    resid     = pred - obs
    return dict(n=n, obs_mean=obs_mean, pred_mean=pred_mean,
                bias_pct=bias_pct,
                RMSE=np.sqrt((resid**2).mean()),
                MAE=np.abs(resid).mean())

# ─────────────────────────────────────────────────────────────────────────────
# Load and join all 4 fix scorecards
all_joined = []
for var in VARIANT_ORDER:
    fix_path = os.path.join(SC_DIR, f'scorecard_{var}_fix.csv')
    print(f"\n{'='*60}")
    print(f"Processing variant: {var}")
    print(f"{'='*60}")
    sc = pd.read_csv(fix_path)
    sc['PLOT']      = sc['PLOT'].astype(str)
    sc['MEASYEAR2'] = sc['MEASYEAR2'].astype(int)
    sc['STATE']     = sc['STATE'].astype(int)
    print(f"  {len(sc)} rows, {sc[['PLOT','MEASYEAR2']].drop_duplicates().shape[0]} plots, "
          f"states={sorted(sc['STATE'].unique())}", flush=True)

    obs_vol = get_obs_volumes(sc)
    print(f"\n  Total obs plot-period records: {len(obs_vol)}", flush=True)

    joined = sc.merge(obs_vol, on=['PLOT','MEASYEAR2'], how='left')
    n_match = joined['MCuFt_OBS'].notna().sum()
    print(f"  Join: {n_match}/{len(joined)} rows have observed volume "
          f"({100*n_match/len(joined):.1f}%)", flush=True)
    all_joined.append(joined)

all_df = pd.concat(all_joined, ignore_index=True)

# Save joined file
out_joined = os.path.join(SC_DIR, 'scorecard_fix_vol_all.csv')
all_df.to_csv(out_joined, index=False)
print(f"\nSaved joined data -> {out_joined}", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Metrics
jv = all_df[all_df['MCuFt_OBS'].notna()].copy()
print(f"\nComputing metrics on {len(jv)} rows...", flush=True)

rows = []
for (variant, config), grp in jv.groupby(['variant','config']):
    cf = metrics(grp, 'MCuFt_PRED', 'MCuFt_OBS')
    bf = metrics(grp, 'BdFt_PRED',  'BdFt_OBS')
    ba = metrics(grp, 'BA_PRED',    'BA_OBS')
    rows.append(dict(
        variant=variant, config=config, n=cf['n'],
        MCuFt_obs_mean=cf['obs_mean'],  MCuFt_pred_mean=cf['pred_mean'],
        MCuFt_bias_pct=cf['bias_pct'],  MCuFt_RMSE=cf['RMSE'], MCuFt_MAE=cf['MAE'],
        BdFt_obs_mean =bf['obs_mean'],  BdFt_pred_mean =bf['pred_mean'],
        BdFt_bias_pct =bf['bias_pct'],  BdFt_RMSE=bf['RMSE'],  BdFt_MAE=bf['MAE'],
        BA_obs_mean   =ba['obs_mean'],  BA_pred_mean   =ba['pred_mean'],
        BA_bias_pct   =ba['bias_pct'],  BA_RMSE=ba['RMSE'],
    ))

summary = pd.DataFrame(rows)
summary['variant'] = pd.Categorical(summary['variant'], categories=VARIANT_ORDER, ordered=True)
summary['config']  = pd.Categorical(summary['config'],  categories=CONFIG_ORDER,  ordered=True)
summary = summary.sort_values(['variant','config']).reset_index(drop=True)

out_summary = os.path.join(SC_DIR, 'scorecard_fix_vol_summary.csv')
summary.to_csv(out_summary, index=False)
print(f"Saved summary -> {out_summary}", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Print pivots
def print_pivot(summary, vcol, title, fmt='{:.1f}'):
    print(f"\n{'='*76}")
    print(f"  {title}")
    print(f"{'='*76}")
    piv = (summary.pivot(index='variant', columns='config', values=vcol)
                  .reindex(index=VARIANT_ORDER, columns=CONFIG_ORDER))
    cw = 14
    header = f"{'variant':<9}" + ''.join(c.rjust(cw) for c in CONFIG_ORDER)
    print(header)
    print('-'*len(header))
    for var in VARIANT_ORDER:
        if var not in piv.index:
            print(f"{var:<9}{'(no data)':>14}")
            continue
        vals = ''.join(
            (fmt.format(v) if pd.notna(v) else 'NA').rjust(cw)
            for v in piv.loc[var]
        )
        print(f"{var:<9}{vals}")

print_pivot(summary, 'MCuFt_bias_pct', 'MCuFt Bias%  (+ = overprediction)')
print_pivot(summary, 'MCuFt_RMSE',     'MCuFt RMSE  (ft3 ac-1)')
print_pivot(summary, 'BdFt_bias_pct',  'BdFt  Bias%  (+ = overprediction)')
print_pivot(summary, 'BdFt_RMSE',      'BdFt  RMSE  (bd ft ac-1)')
print_pivot(summary, 'BA_bias_pct',    'BA    Bias%  (+ = overprediction)')
print_pivot(summary, 'BA_RMSE',        'BA    RMSE  (ft2 ac-1)')

# Mean observed by variant
print(f"\n{'='*76}")
print("  Mean Observed MCuFt ac-1 by Corrected Variant")
print(f"{'='*76}")
obs_ctx = (jv.drop_duplicates(subset=['PLOT','MEASYEAR2','variant'])
             .groupby('variant')['MCuFt_OBS']
             .agg(mean='mean', median='median', sd='std', n='count')
             .reindex(VARIANT_ORDER))
print(obs_ctx.to_string(float_format='{:.1f}'.format))

# EM comparison
print(f"\n{'='*76}")
print("  EM Corrected (MT+SD) vs EM Original (IA/MN/WI) Comparison")
print(f"{'='*76}")
em_data = summary[summary['variant'] == 'em'][
    ['config','n','MCuFt_obs_mean','MCuFt_pred_mean','MCuFt_bias_pct','MCuFt_RMSE',
     'BA_pred_mean','BA_bias_pct']
].copy()
print("Old EM (Midwest): BA bias -16.18%, MCuFt mean pred ~1,603")
print()
print(em_data.to_string(index=False, float_format='{:.1f}'.format))

# Best combos
best_cf = summary.loc[summary['MCuFt_RMSE'].idxmin()]
best_bf = summary.loc[summary['BdFt_RMSE'].idxmin()]
print(f"\n{'='*76}")
print("  Best Combinations (corrected variants only)")
print(f"{'='*76}")
print(f"Best MCuFt RMSE: variant={best_cf['variant']}, config={best_cf['config']}, "
      f"RMSE={best_cf['MCuFt_RMSE']:.1f}  (bias {best_cf['MCuFt_bias_pct']:+.1f}%)")
print(f"Best BdFt  RMSE: variant={best_bf['variant']}, config={best_bf['config']}, "
      f"RMSE={best_bf['BdFt_RMSE']:.1f}  (bias {best_bf['BdFt_bias_pct']:+.1f}%)")

print("\nDone.", flush=True)
