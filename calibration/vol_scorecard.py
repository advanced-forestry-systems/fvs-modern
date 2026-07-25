#!/usr/bin/env python3
"""
Volume accuracy scorecard for FVS fvs-modern 7-config x 5-variant comparison.

Join logic (confirmed from FIA schema):
  scorecard.PLOT  = PLOT.CN at MEASYEAR1 (the first-period plot record)
  We need trees at MEASYEAR2, which have a DIFFERENT PLT_CN.

Correct path:
  1. Look up each scorecard PLOT in the state PLOT table -> get (COUNTYCD, PLOT_number)
  2. Find the PLOT.CN where (STATECD, COUNTYCD, PLOT_number, MEASYEAR == MEASYEAR2)
  3. Filter TREE on those MEASYEAR2 PLT_CNs

Skip EM variant (state assignments still being corrected).
"""

import pandas as pd
import numpy as np
import os
import sys

STATE_ABBR = {
    1:'AL', 5:'AR', 6:'CA', 8:'CO', 12:'FL', 13:'GA',
    16:'ID', 21:'KY', 22:'LA', 23:'ME', 24:'MD', 25:'MA',
    28:'MS', 30:'MT', 32:'NV', 33:'NH', 34:'NJ', 35:'NM',
    36:'NY', 37:'NC', 40:'OK', 41:'OR', 42:'PA', 45:'SC',
    47:'TN', 48:'TX', 49:'UT', 50:'VT', 51:'VA', 53:'WA',
    54:'WV', 56:'WY'
}

FIA_DIR        = '/fs/scratch/PUOM0008/crsfaaron/FIA'
SCORECARD_PATH = '/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard/scorecard_all_variants.csv'
OUT_JOINED     = '/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard/scorecard_vol_all.csv'
OUT_SUMMARY    = '/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard/scorecard_vol_summary.csv'

CONFIG_ORDER  = ["default","calibrated","greg_dg","hg_bayes","hg_organon","hg_conus_sp","hg_conus_clim"]
VARIANT_ORDER = ["ne","sn","wc","pn"]

# ── 1. Load and filter scorecard ─────────────────────────────────────────────
print("Loading scorecard...", flush=True)
sc = pd.read_csv(SCORECARD_PATH)
sc = sc[sc['variant'] != 'em'].copy()
sc['PLOT']      = sc['PLOT'].astype(str)
sc['MEASYEAR2'] = sc['MEASYEAR2'].astype(int)
sc['STATE']     = sc['STATE'].astype(int)
print(f"  {len(sc)} rows (EM dropped), variants: {sorted(sc['variant'].unique())}", flush=True)

states_needed = sorted(sc['STATE'].unique())
print(f"  States: {states_needed}", flush=True)

PLOT_COLS = ['CN','STATECD','COUNTYCD','PLOT','MEASYEAR']
TREE_COLS = ['PLT_CN','INVYR','STATECD','STATUSCD','VOLCFNET','VOLBFNET','TPA_UNADJ']

obs_rows = []

for state in states_needed:
    abbr = STATE_ABBR.get(state)
    if not abbr:
        print(f"WARNING: no abbreviation for STATECD={state}", flush=True)
        continue

    plot_file = os.path.join(FIA_DIR, f'{abbr}_PLOT.csv')
    tree_file = os.path.join(FIA_DIR, f'{abbr}_TREE.csv')
    if not os.path.exists(plot_file):
        print(f"WARNING: {plot_file} not found, skipping state {state}", flush=True)
        continue
    if not os.path.exists(tree_file):
        print(f"WARNING: {tree_file} not found, skipping state {state}", flush=True)
        continue

    state_sc = sc[sc['STATE'] == state][['PLOT','MEASYEAR2']].drop_duplicates()
    measyear1_cns  = set(state_sc['PLOT'].astype(str))
    needed_myr2    = set(state_sc['MEASYEAR2'].astype(int))

    print(f"\nState {state} ({abbr}): {len(measyear1_cns)} MEASYEAR1 plot CNs, "
          f"target MEASYEAR2 = {sorted(needed_myr2)}", flush=True)

    # ── Step A: Read PLOT table to build two maps ─────────────────────────
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

    # Map A: MEASYEAR1 CN -> (COUNTYCD, PLOT_number)
    m1_rows = plot_df[plot_df['CN'].isin(measyear1_cns)][
        ['CN','COUNTYCD','PLOT']].copy()
    cn_to_phys = {row.CN: (row.COUNTYCD, row.PLOT) for _, row in m1_rows.iterrows()}
    print(f"  Found physical IDs for {len(cn_to_phys)}/{len(measyear1_cns)} MEASYEAR1 CNs", flush=True)

    if not cn_to_phys:
        print(f"  No MEASYEAR1 CNs matched in PLOT table for state {state}", flush=True)
        continue

    # Map B: (COUNTYCD, PLOT_number, MEASYEAR) -> CN for MEASYEAR2 years
    needed_physplots = set(cn_to_phys.values())
    m2_candidates = plot_df[
        plot_df['MEASYEAR'].isin(needed_myr2)
    ].copy()
    m2_candidates['phys'] = list(zip(m2_candidates['COUNTYCD'], m2_candidates['PLOT']))
    m2_candidates = m2_candidates[m2_candidates['phys'].isin(needed_physplots)]

    # (COUNTYCD, PLOT_number, MEASYEAR) -> CN
    phys_yr_to_cn2 = {
        (row.COUNTYCD, row.PLOT, int(row.MEASYEAR)): row.CN
        for _, row in m2_candidates.iterrows()
        if pd.notna(row.MEASYEAR)
    }
    print(f"  Built {len(phys_yr_to_cn2)} (physplot, MEASYEAR2) -> CN2 mappings", flush=True)

    # Build the MEASYEAR2 PLT_CNs we need (one per scorecard row in this state)
    # Also keep the scorecard PLOT (MEASYEAR1 CN) for the final join key
    sc_to_cn2 = {}   # scorecard PLOT -> list of (MEASYEAR2, CN2)
    no_match_count = 0
    for _, row in state_sc.iterrows():
        m1_cn   = str(row['PLOT'])
        myr2    = int(row['MEASYEAR2'])
        phys    = cn_to_phys.get(m1_cn)
        if phys is None:
            no_match_count += 1
            continue
        cn2 = phys_yr_to_cn2.get((phys[0], phys[1], myr2))
        if cn2 is None:
            no_match_count += 1
            continue
        sc_to_cn2[(m1_cn, myr2)] = cn2

    print(f"  Resolved {len(sc_to_cn2)} scorecard plot-periods to MEASYEAR2 PLT_CNs "
          f"({no_match_count} unresolved)", flush=True)

    if not sc_to_cn2:
        print(f"  No MEASYEAR2 PLT_CNs found for state {state}", flush=True)
        continue

    needed_cn2s = set(sc_to_cn2.values())

    # ── Step B: Read TREE table, filter to MEASYEAR2 PLT_CNs ─────────────
    try:
        chunks = []
        for chunk in pd.read_csv(tree_file,
                                  usecols=lambda c: c in TREE_COLS,
                                  chunksize=400_000,
                                  low_memory=False):
            chunk['PLT_CN'] = chunk['PLT_CN'].astype(str)
            filt = chunk[chunk['PLT_CN'].isin(needed_cn2s)]
            if len(filt) > 0:
                chunks.append(filt)
    except Exception as exc:
        print(f"  ERROR reading {tree_file}: {exc}", flush=True)
        continue

    if not chunks:
        print(f"  No tree records found for MEASYEAR2 PLT_CNs in state {state}", flush=True)
        continue

    trees = pd.concat(chunks, ignore_index=True)
    print(f"  {len(trees)} raw tree records (STATUSCD all)", flush=True)

    # Live trees only
    trees = trees[trees['STATUSCD'] == 1].copy()

    for col in ['VOLCFNET','VOLBFNET','TPA_UNADJ']:
        trees[col] = pd.to_numeric(trees[col], errors='coerce').fillna(0.0)

    trees['cf_contrib'] = trees['VOLCFNET'] * trees['TPA_UNADJ']
    trees['bf_contrib'] = trees['VOLBFNET'] * trees['TPA_UNADJ']

    # Aggregate by PLT_CN (= MEASYEAR2 PLT_CN)
    agg = (trees.groupby('PLT_CN')
               .agg(MCuFt_OBS=('cf_contrib','sum'),
                    BdFt_OBS =('bf_contrib','sum'))
               .reset_index())
    agg['PLT_CN'] = agg['PLT_CN'].astype(str)

    # Reverse map: CN2 -> (scorecard PLOT, MEASYEAR2)
    cn2_to_sc = {v: k for k, v in sc_to_cn2.items()}  # CN2 -> (m1_cn, myr2)
    agg['sc_PLOT']     = agg['PLT_CN'].map(lambda x: cn2_to_sc.get(x, (None,None))[0])
    agg['sc_MEASYEAR2'] = agg['PLT_CN'].map(lambda x: cn2_to_sc.get(x, (None,None))[1])

    agg = agg.dropna(subset=['sc_PLOT'])
    obs_rows.append(agg[['sc_PLOT','sc_MEASYEAR2','MCuFt_OBS','BdFt_OBS']])
    print(f"  Aggregated to {len(agg)} plot-year records with volume obs", flush=True)

if not obs_rows:
    print("ERROR: no observed volume data found -- aborting", flush=True)
    sys.exit(1)

obs_vol = pd.concat(obs_rows, ignore_index=True)
obs_vol = obs_vol.rename(columns={'sc_PLOT':'PLOT','sc_MEASYEAR2':'MEASYEAR2'})
obs_vol['PLOT']      = obs_vol['PLOT'].astype(str)
obs_vol['MEASYEAR2'] = obs_vol['MEASYEAR2'].astype(int)
print(f"\nTotal observed plot-year records: {len(obs_vol)}", flush=True)

# ── 3. Join to scorecard ──────────────────────────────────────────────────────
joined = sc.merge(
    obs_vol[['PLOT','MEASYEAR2','MCuFt_OBS','BdFt_OBS']],
    on=['PLOT','MEASYEAR2'],
    how='left'
)

n_matched = joined['MCuFt_OBS'].notna().sum()
n_total   = len(joined)
print(f"Join: {n_matched}/{n_total} rows have observed volume ({100*n_matched/n_total:.1f}%)", flush=True)

# Diagnose missing by state
miss = joined[joined['MCuFt_OBS'].isna()].groupby(['STATE','variant']).size().reset_index(name='n_missing')
if len(miss) > 0:
    print("Missing obs by state/variant:")
    print(miss.to_string())

print(f"Saving joined file -> {OUT_JOINED}", flush=True)
joined.to_csv(OUT_JOINED, index=False)

# ── 4. Metrics ────────────────────────────────────────────────────────────────
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
    bias_pct  = 100*(pred_mean-obs_mean)/obs_mean if obs_mean else np.nan
    resid     = pred - obs
    return dict(n=n, obs_mean=obs_mean, pred_mean=pred_mean,
                bias_pct=bias_pct,
                RMSE=np.sqrt((resid**2).mean()),
                MAE =np.abs(resid).mean())

jv = joined[joined['MCuFt_OBS'].notna()].copy()
print(f"\nMetrics on {len(jv)} rows with observed volume", flush=True)

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

print(f"Saving summary -> {OUT_SUMMARY}", flush=True)
summary.to_csv(OUT_SUMMARY, index=False)

# ── 5. Print pivot tables ─────────────────────────────────────────────────────
def print_pivot(summary, vcol, title, fmt='{:.1f}'):
    print(f"\n{'='*76}")
    print(f"  {title}")
    print(f"{'='*76}")
    piv = (summary.pivot(index='variant', columns='config', values=vcol)
                  .reindex(index=VARIANT_ORDER, columns=CONFIG_ORDER))
    cw = 13
    header = f"{'variant':<9}" + ''.join(c.rjust(cw) for c in CONFIG_ORDER)
    print(header)
    print('-'*len(header))
    for var in VARIANT_ORDER:
        if var not in piv.index:
            print(f"{var:<9}{'(no data)':>13}")
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

# Mean observed MCuFt by variant (one row per physical plot-period)
print(f"\n{'='*76}")
print("  Mean Observed MCuFt ac-1 by Variant")
print(f"{'='*76}")
obs_ctx = (jv.drop_duplicates(subset=['PLOT','MEASYEAR2'])
             .groupby('variant')['MCuFt_OBS']
             .agg(mean='mean', median='median', sd='std', n='count')
             .reindex(VARIANT_ORDER))
print(obs_ctx.to_string(float_format='{:.1f}'.format))

# Best combinations
best_cf = summary.loc[summary['MCuFt_RMSE'].idxmin()]
best_bf = summary.loc[summary['BdFt_RMSE'].idxmin()]
print(f"\n{'='*76}")
print("  Best Combinations")
print(f"{'='*76}")
print(f"Best MCuFt RMSE : variant={best_cf['variant']}, config={best_cf['config']}, "
      f"RMSE={best_cf['MCuFt_RMSE']:.1f} ft3 ac-1  (bias {best_cf['MCuFt_bias_pct']:+.1f}%)")
print(f"Best BdFt  RMSE : variant={best_bf['variant']}, config={best_bf['config']}, "
      f"RMSE={best_bf['BdFt_RMSE']:.1f} bd ft ac-1  (bias {best_bf['BdFt_bias_pct']:+.1f}%)")

print("\nDone.", flush=True)
