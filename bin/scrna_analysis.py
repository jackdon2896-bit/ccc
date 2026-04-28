#!/usr/bin/env python3
"""Integrate multiple scRNA-seq H5AD files, run full analysis, export markers."""

import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ads',       nargs='+', required=True)
    p.add_argument('--out_h5ad',    required=True)
    p.add_argument('--out_plots',   required=True)
    p.add_argument('--out_metrics', required=True)
    p.add_argument('--out_markers', required=True)
    p.add_argument('--n_top_genes', type=int, default=2000)
    p.add_argument('--n_pcs',       type=int, default=50)
    p.add_argument('--resolution',  type=float, default=0.5)
    return p.parse_args()

def main():
    import scanpy as sc
    args = parse_args()

    adatas = [sc.read_h5ad(f) for f in args.h5ads]
    print(f"[scrna_analysis] Concatenating {len(adatas)} datasets")

    if len(adatas) > 1:
        import anndata
        adata = anndata.concat(adatas, label='batch',
                               keys=[f'batch{i}' for i in range(len(adatas))],
                               merge='same')
    else:
        adata = adatas[0].copy()
    adata.var_names_make_unique()

    # Normalise (might already be done, guard with raw check)
    if adata.X.max() > 100:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, flavor='seurat_v3',
                                  batch_key='batch' if 'batch' in adata.obs else None)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(args.n_pcs, adata.n_obs-1, adata.n_vars-1))

    # Harmony batch correction if multiple batches
    if 'batch' in adata.obs and adata.obs['batch'].nunique() > 1:
        try:
            sc.external.pp.harmony_integrate(adata, 'batch')
            use_rep = 'X_pca_harmony'
        except Exception:
            use_rep = 'X_pca'
    else:
        use_rep = 'X_pca'

    sc.pp.neighbors(adata, use_rep=use_rep, n_pcs=args.n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution)
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Export marker genes
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(args.out_markers, index=False)

    # Plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Clusters')
    sc.pl.umap(adata, color='brain_region' if 'brain_region' in adata.obs else 'batch',
               ax=axes[1], show=False, title='Brain region / batch')
    sc.pl.umap(adata, color='total_counts', ax=axes[2], show=False, title='Counts')
    plt.tight_layout()
    plt.savefig(args.out_plots, dpi=150); plt.close()

    adata.write_h5ad(args.out_h5ad, compression='gzip')
    metrics = dict(n_cells=int(adata.n_obs), n_genes=int(adata.n_vars),
                   n_clusters=int(adata.obs['leiden'].nunique()),
                   n_batches=len(adatas))
    Path(args.out_metrics).write_text(json.dumps(metrics, indent=2))
    print(f"[scrna_analysis] Done: {adata.shape}")

if __name__ == '__main__':
    main()
