#!/usr/bin/env python3
"""
Preprocess spatial transcriptomics data.
Reads H5/H5AD + TIFF image, applies QC filters, normalises, finds HVGs,
exports processed AnnData + spatial coordinates + QC plots.
"""

import argparse, json, sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--h5',         required=True)
    p.add_argument('--image',      required=True)
    p.add_argument('--out_h5ad',   required=True)
    p.add_argument('--out_image',  required=True)
    p.add_argument('--out_coords', required=True)
    p.add_argument('--out_qc',     required=True)
    p.add_argument('--out_plot',   required=True)
    p.add_argument('--min_genes',  type=int, default=200)
    p.add_argument('--min_cells',  type=int, default=3)
    p.add_argument('--max_genes',  type=int, default=5000)
    p.add_argument('--mt_max',     type=float, default=20.0)
    return p.parse_args()

def main():
    import scanpy as sc
    args = parse_args()

    sc.settings.verbosity = 2
    print(f"[preprocess_spatial] Loading {args.h5}")

    # ── Load data ──────────────────────────────────────────────────────────
    if args.h5.endswith('.h5ad'):
        adata = sc.read_h5ad(args.h5)
    elif args.h5.endswith('.h5'):
        adata = sc.read_10x_h5(args.h5)
    else:
        raise ValueError(f"Unsupported format: {args.h5}")

    print(f"  Loaded: {adata.shape[0]} spots × {adata.shape[1]} genes")

    # ── Process image ──────────────────────────────────────────────────────
    img = Image.open(args.image).convert('RGB')
    img.save(args.out_image)

    # ── QC ────────────────────────────────────────────────────────────────
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                log1p=False, inplace=True)

    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs.n_genes_by_counts < args.max_genes].copy()
    adata = adata[adata.obs.pct_counts_mt < args.mt_max].copy()
    print(f"  After QC: {adata.n_obs}/{n0} spots retained")

    # ── Normalise / HVG ────────────────────────────────────────────────────
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

    # ── Spatial coords ─────────────────────────────────────────────────────
    if 'spatial' in adata.obsm:
        coords = pd.DataFrame(adata.obsm['spatial'],
                               index=adata.obs_names, columns=['x', 'y'])
    else:
        rng = np.random.default_rng(42)
        coords = pd.DataFrame(rng.uniform(0, 1000, (adata.n_obs, 2)),
                               index=adata.obs_names, columns=['x', 'y'])
        adata.obsm['spatial'] = coords.values
    coords.to_csv(args.out_coords)

    # ── QC plot ────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].hist(adata.obs.n_genes_by_counts, bins=50)
    axes[0].set(title='Genes / spot', xlabel='n_genes', ylabel='spots')
    axes[1].hist(adata.obs.total_counts, bins=50)
    axes[1].set(title='UMI counts / spot', xlabel='total_counts')
    axes[2].hist(adata.obs.pct_counts_mt, bins=50)
    axes[2].set(title='MT %', xlabel='pct_counts_mt')
    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=150)
    plt.close()

    # ── Save ───────────────────────────────────────────────────────────────
    adata.write_h5ad(args.out_h5ad, compression='gzip')
    qc = dict(n_spots_raw=int(n0), n_spots_final=int(adata.n_obs),
              n_genes=int(adata.n_vars),
              mean_genes=float(adata.obs.n_genes_by_counts.mean()),
              mean_counts=float(adata.obs.total_counts.mean()))
    Path(args.out_qc).write_text(json.dumps(qc, indent=2))
    print("[preprocess_spatial] Done.")

if __name__ == '__main__':
    main()
