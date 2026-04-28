#!/usr/bin/env python3
"""Full spatial transcriptomics analysis: PCA, UMAP, clustering, Squidpy spatial stats."""

import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad');        p.add_argument('--mask')
    p.add_argument('--out_h5ad');    p.add_argument('--out_graph')
    p.add_argument('--out_plots');   p.add_argument('--out_metrics')
    p.add_argument('--n_top_genes', type=int, default=2000)
    p.add_argument('--n_pcs',       type=int, default=50)
    p.add_argument('--resolution',  type=float, default=0.5)
    p.add_argument('--n_neighbors', type=int, default=15)
    return p.parse_args()

def main():
    import scanpy as sc
    try:
        import squidpy as sq
        HAS_SQ = True
    except ImportError:
        HAS_SQ = False
    args = parse_args()

    adata = sc.read_h5ad(args.h5ad)
    print(f"[spatial_analysis] {adata.shape}")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(args.n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    sc.pp.neighbors(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution)

    # Spatial graph via squidpy if available
    if HAS_SQ and 'spatial' in adata.obsm:
        sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
        sq.gr.nhood_enrichment(adata, cluster_key='leiden')
        sq.gr.spatial_autocorr(adata, mode='moran')
        adj = adata.obsp['spatial_connectivities']
        import scipy.sparse as sp
        graph_df = pd.DataFrame.sparse.from_spmatrix(adj,
                    index=adata.obs_names, columns=adata.obs_names)
    else:
        # Fallback: k-NN distance graph from spatial coords
        if 'spatial' in adata.obsm:
            from sklearn.neighbors import kneighbors_graph
            G = kneighbors_graph(adata.obsm['spatial'], n_neighbors=args.n_neighbors)
            graph_df = pd.DataFrame(G.toarray(), index=adata.obs_names,
                                    columns=adata.obs_names)
        else:
            graph_df = pd.DataFrame(index=adata.obs_names, columns=adata.obs_names)
    graph_df.to_csv(args.out_graph)

    # Marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    sc.pl.umap(adata, color='leiden', ax=axes[0,0], show=False, title='Clusters')
    sc.pl.umap(adata, color='total_counts', ax=axes[0,1], show=False, title='Counts')
    if 'spatial' in adata.obsm:
        sc.pl.spatial(adata, color='leiden', ax=axes[1,0], show=False, title='Spatial clusters')
        sc.pl.spatial(adata, color='total_counts', ax=axes[1,1], show=False, title='Spatial counts')
    plt.tight_layout()
    plt.savefig(args.out_plots, dpi=150); plt.close()

    adata.write_h5ad(args.out_h5ad, compression='gzip')
    n_clusters = int(adata.obs['leiden'].nunique())
    metrics = dict(n_spots=int(adata.n_obs), n_genes=int(adata.n_vars),
                   n_clusters=n_clusters)
    Path(args.out_metrics).write_text(json.dumps(metrics, indent=2))
    print(f"[spatial_analysis] {n_clusters} clusters")

if __name__ == '__main__':
    main()
