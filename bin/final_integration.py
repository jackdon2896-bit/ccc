#!/usr/bin/env python3
"""
Golden Standard Multi-Modal Integration.
Integrates spatial, scRNA-seq, ML predictions, GNN embeddings,
pseudotime and CCC networks into a single annotated AnnData.

Integration strategy:
  spatial ↔ scRNA   : correlation-based transfer (Squidpy OT when available)
  ML / GNN           : direct obs column injection
  CCC                : sender/receiver scores per cluster → per cell
  Trajectory         : pseudotime column transfer
"""

import argparse, json, glob, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--result_files',    nargs='+', required=True,
                   help='flat list of result files emitted by upstream processes')
    p.add_argument('--original_image',  required=True)
    p.add_argument('--out_h5ad',        required=True)
    p.add_argument('--out_plots_dir',   required=True)
    p.add_argument('--out_metrics',     required=True)
    p.add_argument('--genome',          default='mm10')
    return p.parse_args()

def classify_file(f):
    fl = f.lower()
    if 'spatial_analyzed' in fl:   return 'spatial'
    if 'scrna_integrated' in fl:   return 'scrna'
    if 'ml_results'        in fl:   return 'ml_h5ad'
    if 'ml_predictions'    in fl:   return 'ml_csv'
    if 'gnn_results'       in fl:   return 'gnn_h5ad'
    if 'gnn_embeddings'    in fl:   return 'gnn_csv'
    if 'trajectory_results' in fl:  return 'traj_h5ad'
    if 'pseudotime'        in fl:   return 'traj_csv'
    if 'cellchat_networks' in fl:   return 'cellchat'
    if 'liana_results'     in fl:   return 'liana'
    if 'commot_flows'      in fl:   return 'commot'
    return 'unknown'

def add_col_safe(adata, col, series):
    common = series.index.intersection(adata.obs_names)
    if len(common):
        adata.obs.loc[common, col] = series.loc[common]

def transfer_cell_types(adata_spatial, adata_scrna):
    """Correlation-based scRNA → spatial cell type transfer."""
    try:
        common_genes = adata_spatial.var_names.intersection(adata_scrna.var_names)
        if len(common_genes) < 50:
            return adata_spatial
        import scipy.sparse as sp
        Xs = adata_spatial[:, common_genes].X
        Xr = adata_scrna[:, common_genes].X
        if sp.issparse(Xs): Xs = Xs.toarray()
        if sp.issparse(Xr): Xr = Xr.toarray()
        # Pearson-like: row-normalise then dot
        Xs_n = Xs / (np.linalg.norm(Xs, axis=1, keepdims=True) + 1e-8)
        Xr_n = Xr / (np.linalg.norm(Xr, axis=1, keepdims=True) + 1e-8)
        sim = Xs_n @ Xr_n.T  # (n_spatial, n_scrna)
        best = sim.argmax(axis=1)
        groupby = 'leiden' if 'leiden' in adata_scrna.obs else adata_scrna.obs.columns[0]
        labels = adata_scrna.obs[groupby].values
        adata_spatial.obs['scrna_mapped_cluster'] = [str(labels[b]) for b in best]
        adata_spatial.obs['scrna_mapping_score']  = sim.max(axis=1)
    except Exception as e:
        print(f"  [transfer] Warning: {e}")
    return adata_spatial

def main():
    import scanpy as sc
    args = parse_args()
    Path(args.out_plots_dir).mkdir(parents=True, exist_ok=True)

    print("[final_integration] Classifying result files...")
    buckets = {}
    for f in args.result_files:
        k = classify_file(f)
        buckets.setdefault(k, []).append(f)

    # ── Load base spatial AnnData ──────────────────────────────────────────
    spatial_files = buckets.get('spatial', [])
    if not spatial_files:
        # fallback: use any h5ad
        spatial_files = [f for f in args.result_files if f.endswith('.h5ad')]
    if not spatial_files:
        raise RuntimeError("No spatial AnnData found among result files.")
    adata = sc.read_h5ad(spatial_files[0])
    print(f"  Base spatial: {adata.shape}")

    # ── Transfer scRNA-seq labels ──────────────────────────────────────────
    for f in buckets.get('scrna', []):
        print(f"  Integrating scRNA: {f}")
        adata_scrna = sc.read_h5ad(f)
        adata = transfer_cell_types(adata, adata_scrna)

    # ── Inject ML predictions ──────────────────────────────────────────────
    for f in buckets.get('ml_csv', []):
        ml = pd.read_csv(f)
        if 'cell_id' in ml.columns:
            ml = ml.set_index('cell_id')
        for col in ['predicted_label', 'ml_confidence']:
            if col in ml.columns:
                add_col_safe(adata, f'ml_{col}', ml[col])

    # ── GNN embeddings ─────────────────────────────────────────────────────
    for f in buckets.get('gnn_csv', []):
        emb = pd.read_csv(f, index_col=0)
        emb_cols = [c for c in emb.columns if c.startswith('gnn_')]
        if emb_cols:
            common = emb.index.intersection(adata.obs_names)
            mat = np.zeros((adata.n_obs, len(emb_cols)))
            idx_map = {c: i for i, c in enumerate(adata.obs_names)}
            for c in common:
                mat[idx_map[c]] = emb.loc[c, emb_cols].values
            adata.obsm['X_gnn'] = mat

    # ── Pseudotime ─────────────────────────────────────────────────────────
    for f in buckets.get('traj_csv', []):
        pt = pd.read_csv(f, index_col=0)
        if 'pseudotime' in pt.columns:
            add_col_safe(adata, 'pseudotime', pt['pseudotime'])

    # ── CCC scores ─────────────────────────────────────────────────────────
    groupby = 'leiden' if 'leiden' in adata.obs else None
    for method, key in [('cellchat', 'cellchat'), ('liana', 'liana'), ('commot', 'commot')]:
        for f in buckets.get(key, []):
            df = pd.read_csv(f)
            if 'sender' in df.columns and 'receiver' in df.columns:
                scol = 'strength' if 'strength' in df.columns else (
                       'product_score' if 'product_score' in df.columns else
                       'spatial_flow' if 'spatial_flow' in df.columns else None)
                if scol and groupby:
                    sender_s = df.groupby('sender')[scol].sum()
                    receiver_s= df.groupby('receiver')[scol].sum()
                    cct = adata.obs[groupby]
                    adata.obs[f'{method}_sender']   = cct.map(sender_s).fillna(0)
                    adata.obs[f'{method}_receiver'] = cct.map(receiver_s).fillna(0)

    # ── Integration metrics ────────────────────────────────────────────────
    metrics = dict(
        n_spots   = int(adata.n_obs),
        n_genes   = int(adata.n_vars),
        has_scrna_mapping  = 'scrna_mapped_cluster' in adata.obs,
        has_ml_predictions = 'ml_predicted_label' in adata.obs,
        has_gnn_embeddings = 'X_gnn' in adata.obsm,
        has_pseudotime     = 'pseudotime' in adata.obs,
        ccc_methods = sum(1 for col in adata.obs if col.endswith('_sender')),
        genome = args.genome,
    )

    # ── Save ──────────────────────────────────────────────────────────────
    adata.write_h5ad(args.out_h5ad, compression='gzip')
    Path(args.out_metrics).write_text(json.dumps(metrics, indent=2))

    # ── Master plot ────────────────────────────────────────────────────────
    n_cols = min(3, 1 + sum([
        'scrna_mapped_cluster' in adata.obs,
        'ml_ml_predicted_label' in adata.obs or 'ml_predicted_label' in adata.obs,
        'pseudotime' in adata.obs,
    ]))
    fig, axes = plt.subplots(2, n_cols, figsize=(6*n_cols, 10))
    if axes.ndim == 1:
        axes = axes[None, :]

    plot_cols = ['leiden'] if 'leiden' in adata.obs else [adata.obs.columns[0]]
    if 'scrna_mapped_cluster' in adata.obs:   plot_cols.append('scrna_mapped_cluster')
    if 'pseudotime' in adata.obs:             plot_cols.append('pseudotime')
    plot_cols = plot_cols[:n_cols]

    if 'spatial' in adata.obsm:
        for i, col in enumerate(plot_cols):
            sc.pl.spatial(adata, color=col, ax=axes[0, i], show=False, title=col)
    if 'X_umap' in adata.obsm:
        for i, col in enumerate(plot_cols):
            sc.pl.umap(adata, color=col, ax=axes[1, i], show=False, title=f'UMAP {col}')

    plt.suptitle('Golden Standard Multi-Modal Integration', fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{args.out_plots_dir}/integration_summary.png", dpi=150)
    plt.close()

    print(f"[final_integration] Done. Success rate: {sum(metrics[k] for k in metrics if isinstance(metrics[k], bool))}/5")

if __name__ == '__main__':
    main()
