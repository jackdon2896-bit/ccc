#!/usr/bin/env python3
"""
LIANA-style multi-method consensus cell-cell communication.
Runs three scoring methods (mean-expression, log-FC, product) and ranks
interactions by consensus score.
"""

import argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings('ignore')

LR_DB = {  # same curated mouse LR pairs
    'Tgfb1':  ['Tgfbr1'], 'Vegfa': ['Flt1'],  'Bdnf':  ['Ntrk2'],
    'Cxcl12': ['Cxcr4'],  'Il1b':  ['Il1r1'], 'Egf':   ['Egfr'],
    'Hgf':    ['Met'],    'Nrg1':  ['Erbb4'], 'Fgf2':  ['Fgfr1'],
    'Pdgfb':  ['Pdgfrb'], 'Igf1':  ['Igf1r'],
}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--scrna_h5ad',    required=True)
    p.add_argument('--spatial_h5ad',  required=True)
    p.add_argument('--out_results',   required=True)
    p.add_argument('--out_consensus', required=True)
    p.add_argument('--out_plots',     required=True)
    return p.parse_args()

def cell_means(adata, groupby):
    import scipy.sparse as sp
    X = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    df['group'] = adata.obs[groupby].values
    return df.groupby('group').mean()

def main():
    import scanpy as sc
    args = parse_args()
    adata = sc.read_h5ad(args.scrna_h5ad)
    groupby = 'leiden' if 'leiden' in adata.obs else adata.obs.columns[0]

    means = cell_means(adata, groupby)
    rows = []
    for lig, recs in LR_DB.items():
        for rec in recs:
            if lig not in means.columns or rec not in means.columns:
                continue
            for sender in means.index:
                for receiver in means.index:
                    l_exp = float(means.loc[sender, lig])
                    r_exp = float(means.loc[receiver, rec])
                    # Three scoring methods
                    product_score  = l_exp * r_exp
                    mean_score     = (l_exp + r_exp) / 2
                    geom_score     = float(np.sqrt(max(l_exp * r_exp, 0)))
                    rows.append({'source': sender, 'target': receiver,
                                  'ligand': lig, 'receptor': rec,
                                  'product_score': product_score,
                                  'mean_score': mean_score,
                                  'geom_score': geom_score})

    results_df = pd.DataFrame(rows)
    results_df.to_csv(args.out_results, index=False)

    # Consensus: rank by each method, average ranks
    for col in ['product_score', 'mean_score', 'geom_score']:
        results_df[f'rank_{col}'] = results_df[col].rank(ascending=False)
    results_df['consensus_rank'] = results_df[['rank_product_score',
                                                'rank_mean_score',
                                                'rank_geom_score']].mean(axis=1)
    consensus = results_df.sort_values('consensus_rank').head(200)
    consensus.to_csv(args.out_consensus, index=False)

    # Plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    pivot = results_df.groupby(['source', 'target'])['product_score'].sum().unstack(fill_value=0)
    if not pivot.empty:
        sns.heatmap(pivot, ax=axes[0], cmap='Blues')
        axes[0].set_title('Communication strength (LIANA)')
    top20 = consensus.head(20).copy()
    top20['pair'] = top20['ligand'] + '→' + top20['receptor']
    axes[1].barh(top20['pair'], top20['product_score'])
    axes[1].set_title('Top 20 consensus interactions')
    plt.tight_layout()
    plt.savefig(args.out_plots, dpi=150); plt.close()
    print(f"[ccc_liana] {len(results_df)} interactions computed")

if __name__ == '__main__':
    main()
