#!/usr/bin/env python3
"""
Graph Neural Network spatial analysis.
Uses PyTorch + torch_geometric (GCNConv) for cell-type classification on
the spatial graph. Falls back to a pure-numpy spectral approach if
torch_geometric is unavailable.
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
    p.add_argument('--spatial_h5ad',    required=True)
    p.add_argument('--spatial_graph',   required=True)
    p.add_argument('--out_h5ad',        required=True)
    p.add_argument('--out_embeddings',  required=True)
    p.add_argument('--out_plots',       required=True)
    p.add_argument('--out_metrics',     required=True)
    p.add_argument('--n_epochs',  type=int,   default=200)
    p.add_argument('--lr',        type=float, default=0.01)
    p.add_argument('--hidden_dim',type=int,   default=128)
    return p.parse_args()

def spectral_embed(X, adj, hidden_dim=64):
    """Simple spectral embedding fallback (no torch needed)."""
    from sklearn.decomposition import TruncatedSVD
    # Graph-smoothed features: X_smooth = D^{-1} A X
    deg = np.array(adj.sum(axis=1)).flatten() + 1e-8
    X_smooth = adj.dot(X) / deg[:, None]
    svd = TruncatedSVD(n_components=min(hidden_dim, X_smooth.shape[1]-1), random_state=42)
    return svd.fit_transform(X_smooth)

def main():
    import scanpy as sc
    from sklearn.preprocessing import LabelEncoder
    from sklearn.metrics import accuracy_score
    from sklearn.model_selection import train_test_split
    args = parse_args()

    adata = sc.read_h5ad(args.spatial_h5ad)
    graph_df = pd.read_csv(args.spatial_graph, index_col=0)
    print(f"[gnn_analysis] {adata.shape}")

    # Feature matrix (HVG expression)
    hvg = (adata.var['highly_variable'] if 'highly_variable' in adata.var
           else pd.Series(True, index=adata.var_names))
    genes = adata.var_names[hvg][:500]
    X_raw = (adata[:, genes].X.toarray() if hasattr(adata[:, genes].X, 'toarray')
             else np.array(adata[:, genes].X))

    # Adjacency
    import scipy.sparse as sp
    common = [c for c in adata.obs_names if c in graph_df.index and c in graph_df.columns]
    if common:
        adj_np = graph_df.loc[common, common].values.astype(float)
        X_feat = X_raw[[list(adata.obs_names).index(c) for c in common]]
    else:
        adj_np = np.eye(len(adata.obs_names))
        X_feat = X_raw
        common = list(adata.obs_names)
    adj_sp = sp.csr_matrix(adj_np)

    # Try torch_geometric; fall back to spectral
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch_geometric.nn import GCNConv
        from torch_geometric.utils import from_scipy_sparse_matrix

        edge_index, edge_weight = from_scipy_sparse_matrix(adj_sp)
        x_t = torch.tensor(X_feat, dtype=torch.float32)

        le = LabelEncoder()
        labels_raw = adata.obs['leiden'].values if 'leiden' in adata.obs else np.zeros(len(common), dtype=int)
        y_t = torch.tensor(le.fit_transform(labels_raw), dtype=torch.long)
        n_classes = int(y_t.max()) + 1

        class GNN(nn.Module):
            def __init__(self, in_f, hid, out):
                super().__init__()
                self.c1 = GCNConv(in_f, hid)
                self.c2 = GCNConv(hid, out)
            def forward(self, x, ei, ew=None):
                h = F.relu(self.c1(x, ei, ew))
                h = F.dropout(h, 0.5, training=self.training)
                return self.c2(h, ei, ew), h  # logits, embeddings

        model = GNN(X_feat.shape[1], args.hidden_dim, n_classes)
        opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=5e-4)
        idx = np.arange(len(common))
        tr_idx, te_idx = train_test_split(idx, test_size=0.3, random_state=42)
        tr_mask = torch.zeros(len(common), dtype=torch.bool)
        te_mask = torch.zeros(len(common), dtype=torch.bool)
        tr_mask[tr_idx] = True; te_mask[te_idx] = True

        model.train()
        for ep in range(args.n_epochs):
            opt.zero_grad()
            out, _ = model(x_t, edge_index, edge_weight)
            loss = F.cross_entropy(out[tr_mask], y_t[tr_mask])
            loss.backward(); opt.step()

        model.eval()
        with torch.no_grad():
            logits, embeddings = model(x_t, edge_index, edge_weight)
            preds = logits.argmax(dim=1).numpy()
            emb_np = embeddings.numpy()
        acc = accuracy_score(y_t[te_mask].numpy(), preds[te_mask])
        print(f"  GNN test accuracy: {acc:.3f}")
        method = 'GCN'

    except Exception as exc:
        print(f"  torch_geometric unavailable ({exc}), using spectral fallback")
        emb_np = spectral_embed(X_feat, adj_sp, args.hidden_dim)
        from sklearn.ensemble import RandomForestClassifier
        le = LabelEncoder()
        labels_raw = adata.obs['leiden'].values if 'leiden' in adata.obs else np.zeros(len(common), dtype=int)
        y_all = le.fit_transform(labels_raw)
        tr_idx, te_idx = train_test_split(np.arange(len(common)), test_size=0.3, random_state=42)
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(emb_np[tr_idx], y_all[tr_idx])
        preds = clf.predict(emb_np)
        acc = accuracy_score(y_all[te_idx], clf.predict(emb_np[te_idx]))
        print(f"  Spectral+RF test accuracy: {acc:.3f}")
        method = 'Spectral'

    # Save embeddings
    embed_df = pd.DataFrame(emb_np, index=common,
                             columns=[f'gnn_{i}' for i in range(emb_np.shape[1])])
    embed_df.to_csv(args.out_embeddings)

    adata.obs.loc[common, 'gnn_prediction'] = le.inverse_transform(preds) if 'le' in dir() else preds.astype(str)
    adata.obsm['X_gnn'] = np.zeros((adata.n_obs, emb_np.shape[1]))
    for i, c in enumerate(common):
        adata.obsm['X_gnn'][list(adata.obs_names).index(c)] = emb_np[i]
    adata.write_h5ad(args.out_h5ad, compression='gzip')

    # Plots
    from sklearn.manifold import TSNE
    tsne_2d = TSNE(n_components=2, random_state=42).fit_transform(emb_np[:2000] if len(emb_np) > 2000 else emb_np)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(tsne_2d[:, 0], tsne_2d[:, 1], c=preds[:len(tsne_2d)], cmap='tab20', s=4)
    axes[0].set_title(f'{method} embeddings (t-SNE)')
    if 'spatial' in adata.obsm:
        sc.pl.spatial(adata, color='gnn_prediction', ax=axes[1], show=False, title='GNN predictions')
    plt.tight_layout()

    import scanpy as sc
    plt.savefig(args.out_plots, dpi=150); plt.close()

    metrics = dict(method=method, test_accuracy=float(acc),
                   n_cells=int(adata.n_obs), embedding_dim=int(emb_np.shape[1]))
    Path(args.out_metrics).write_text(json.dumps(metrics, indent=2))
    print("[gnn_analysis] Done.")

if __name__ == '__main__':
    main()
