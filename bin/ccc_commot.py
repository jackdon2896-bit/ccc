#!/usr/bin/env python3
"""
COMMOT-style spatial cell-cell communication.
Uses spatial coordinates + distance threshold to weight LR interactions
by cell proximity, producing directional communication flows.
"""

import argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

LR_DB = {
    'Tgfb1': ['Tgfbr1'], 'Vegfa': ['Flt1'], 'Bdnf': ['Ntrk2'],
    'Cxcl12': ['Cxcr4'], 'Egf': ['Egfr'],   'Hgf':  ['Met'],
}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--spatial_h5ad',  required=True)
    p.add_argument('--spatial_graph', required=True)
    p.add_argument('--out_h5ad',      required=True)
    p.add_argument('--out_flows',     required=True)
    p.add_argument('--out_plots',     required=True)
    p.add_argument('--distance_thr',  type=float, default=500)
    return p.parse_args()

def main():
    import scanpy as sc
    try:
        import squidpy as sq
        HAS_SQ = True
    except ImportError:
        HAS_SQ = False
    args = parse_args()

    adata = sc.read_h5ad(args.spatial_h5ad)
    print(f"[ccc_commot] {adata.shape}")

    if 'spatial' not in adata.obsm:
        graph_df = pd.read_csv(args.spatial_graph, index_col=0)
        adata.obsm['spatial'] = np.zeros((adata.n_obs, 2))

    coords = adata.obsm['spatial'][:, :2]  # (n_spots, 2)
    n = adata.n_obs
    groupby = 'leiden' if 'leiden' in adata.obs else adata.obs.columns[0]

    # Build distance-gated communication flows
    import scipy.sparse as sp
    flows_rows = []

    for lig, recs in LR_DB.items():
        for rec in recs:
            if lig not in adata.var_names or rec not in adata.var_names:
                continue
            l_exp = np.array(adata[:, lig].X.toarray()).flatten() if hasattr(adata[:, lig].X, 'toarray') else np.array(adata[:, lig].X).flatten()
            r_exp = np.array(adata[:, rec].X.toarray()).flatten() if hasattr(adata[:, rec].X, 'toarray') else np.array(adata[:, rec].X).flatten()

            # Vectorised pairwise distances
            diff = coords[:, None, :] - coords[None, :, :]  # (n,n,2)
            dist = np.linalg.norm(diff, axis=2)              # (n,n)
            w = np.exp(-dist / (args.distance_thr / 3)) * (dist < args.distance_thr)

            # Communication strength: sender l_exp * receiver r_exp * spatial_weight
            strength = l_exp[:, None] * r_exp[None, :] * w  # (n,n)

            # Summarise by cluster pairs
            ct = adata.obs[groupby].values
            for ci in np.unique(ct):
                for cj in np.unique(ct):
                    mask_i = ct == ci; mask_j = ct == cj
                    s = float(strength[np.ix_(mask_i, mask_j)].sum())
                    if s > 0:
                        flows_rows.append({'ligand': lig, 'receptor': rec,
                                           'sender': ci, 'receiver': cj,
                                           'spatial_flow': s})

    flows_df = pd.DataFrame(flows_rows)
    flows_df.to_csv(args.out_flows, index=False)

    # Add aggregate sender/receiver scores to adata
    if not flows_df.empty:
        sender_score = flows_df.groupby('sender')['spatial_flow'].sum()
        receiver_score= flows_df.groupby('receiver')['spatial_flow'].sum()
        adata.obs['commot_sender']   = adata.obs[groupby].map(sender_score).fillna(0)
        adata.obs['commot_receiver'] = adata.obs[groupby].map(receiver_score).fillna(0)

    adata.write_h5ad(args.out_h5ad, compression='gzip')

    # Plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    if 'spatial' in adata.obsm:
        sc.pl.spatial(adata, color='commot_sender',   ax=axes[0], show=False, title='COMMOT senders')
        sc.pl.spatial(adata, color='commot_receiver', ax=axes[1], show=False, title='COMMOT receivers')
    plt.tight_layout()
    import scanpy as sc
    plt.savefig(args.out_plots, dpi=150); plt.close()
    print(f"[ccc_commot] {len(flows_df)} spatial flows computed")

if __name__ == '__main__':
    main()
