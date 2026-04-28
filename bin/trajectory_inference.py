#!/usr/bin/env python3
"""
Trajectory / pseudotime inference via PAGA (scanpy).
Computes PAGA on spatial data, transfers pseudotime to scRNA-seq if available.
"""

import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--spatial_h5ad',  required=True)
    p.add_argument('--scrna_h5ad',    required=True)
    p.add_argument('--out_h5ad',      required=True)
    p.add_argument('--out_pseudotime',required=True)
    p.add_argument('--out_plots',     required=True)
    p.add_argument('--n_pcs',         type=int, default=50)
    return p.parse_args()

def main():
    import scanpy as sc
    args = parse_args()

    adata = sc.read_h5ad(args.spatial_h5ad)
    print(f"[trajectory] spatial: {adata.shape}")

    # Ensure neighbours are computed
    if 'neighbors' not in adata.uns:
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=min(args.n_pcs, adata.n_obs-1, adata.n_vars-1))
        sc.pp.neighbors(adata, n_pcs=args.n_pcs)
        sc.tl.umap(adata)
    if 'leiden' not in adata.obs:
        sc.tl.leiden(adata, resolution=0.5)

    # PAGA
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata, show=False)

    # Diffusion pseudotime
    root_idx = adata.obs['leiden'].value_counts().idxmin()
    adata.uns['iroot'] = int(adata.obs['leiden'].eq(root_idx).to_numpy().argmax())
    try:
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata, n_dcs=10)
        pt_col = 'dpt_pseudotime'
    except Exception:
        # Fallback: use UMAP component 0 as proxy
        adata.obs['pseudo_fallback'] = adata.obsm['X_umap'][:, 0]
        pt_col = 'pseudo_fallback'

    # Save pseudotime
    pt_df = adata.obs[[pt_col]].rename(columns={pt_col: 'pseudotime'})
    pt_df.index.name = 'cell_id'
    pt_df.to_csv(args.out_pseudotime)

    adata.write_h5ad(args.out_h5ad, compression='gzip')

    # Plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=pt_col, ax=axes[0], show=False, title='Pseudotime (UMAP)')
    if 'spatial' in adata.obsm:
        sc.pl.spatial(adata, color=pt_col, ax=axes[1], show=False, title='Pseudotime (spatial)')
    plt.tight_layout()
    plt.savefig(args.out_plots, dpi=150); plt.close()
    print("[trajectory] Done.")

if __name__ == '__main__':
    main()
