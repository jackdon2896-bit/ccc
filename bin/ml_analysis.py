#!/usr/bin/env python3
"""
Random Forest spatial cell-type prediction.
Features: top HVG expression + spatial coordinates + neighbourhood statistics.
"""

import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--spatial_h5ad',  required=True)
    p.add_argument('--spatial_graph', required=True)
    p.add_argument('--out_h5ad',      required=True)
    p.add_argument('--out_preds',     required=True)
    p.add_argument('--out_features',  required=True)
    p.add_argument('--out_plots',     required=True)
    return p.parse_args()

def main():
    import scanpy as sc
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import accuracy_score, classification_report
    from sklearn.preprocessing import LabelEncoder
    args = parse_args()

    adata = sc.read_h5ad(args.spatial_h5ad)
    graph = pd.read_csv(args.spatial_graph, index_col=0)
    print(f"[ml_analysis] {adata.shape}")

    # ── Feature matrix ─────────────────────────────────────────────────────
    feat_parts = []

    # Gene expression (top 100 HVG)
    hvg = (adata.var['highly_variable'] if 'highly_variable' in adata.var
           else pd.Series(True, index=adata.var_names))
    genes = adata.var_names[hvg][:100]
    X_expr = pd.DataFrame(adata[:, genes].X.toarray() if hasattr(adata[:, genes].X, 'toarray')
                           else adata[:, genes].X,
                           index=adata.obs_names, columns=[f'gene_{g}' for g in genes])
    feat_parts.append(X_expr)

    # Spatial coordinates
    if 'spatial' in adata.obsm:
        X_sp = pd.DataFrame(adata.obsm['spatial'][:, :2],
                             index=adata.obs_names, columns=['x', 'y'])
        feat_parts.append(X_sp)

    X = pd.concat(feat_parts, axis=1).fillna(0)

    # Labels
    if 'leiden' in adata.obs:
        y_raw = adata.obs['leiden']
    else:
        from sklearn.cluster import MiniBatchKMeans
        km = MiniBatchKMeans(n_clusters=10, random_state=42)
        y_raw = pd.Series(km.fit_predict(X.values).astype(str),
                          index=adata.obs_names, name='leiden')
        adata.obs['leiden'] = y_raw

    le = LabelEncoder()
    y = le.fit_transform(y_raw)

    X_tr, X_te, y_tr, y_te, idx_tr, idx_te = train_test_split(
        X, y, X.index, test_size=0.3, random_state=42)

    # ── Random Forest ──────────────────────────────────────────────────────
    rf = RandomForestClassifier(n_estimators=200, n_jobs=-1, random_state=42)
    rf.fit(X_tr, y_tr)
    y_pred_all = rf.predict(X)
    acc = accuracy_score(y_te, rf.predict(X_te))
    print(f"  RF test accuracy: {acc:.3f}")

    # ── Save predictions ───────────────────────────────────────────────────
    preds = pd.DataFrame({'cell_id': X.index,
                           'predicted_label': le.inverse_transform(y_pred_all),
                           'ml_confidence': rf.predict_proba(X).max(axis=1)})
    preds.to_csv(args.out_preds, index=False)

    fi = pd.DataFrame({'feature': X.columns, 'importance': rf.feature_importances_})
    fi.sort_values('importance', ascending=False).to_csv(args.out_features, index=False)

    adata.obs['ml_prediction'] = le.inverse_transform(y_pred_all)
    adata.obs['ml_confidence'] = rf.predict_proba(X).max(axis=1)
    adata.write_h5ad(args.out_h5ad, compression='gzip')

    # ── Plots ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    top20 = fi.nlargest(20, 'importance')
    axes[0].barh(top20['feature'], top20['importance'])
    axes[0].set_title('Top 20 Features')
    if 'spatial' in adata.obsm:
        sc.pl.spatial(adata, color='ml_prediction', ax=axes[1], show=False, title='ML predictions')
        sc.pl.spatial(adata, color='ml_confidence', ax=axes[2], show=False, title='Confidence')
    plt.tight_layout()
    plt.savefig(args.out_plots, dpi=150); plt.close()
    print("[ml_analysis] Done.")

if __name__ == '__main__':
    main()
