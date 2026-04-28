#!/usr/bin/env python3
"""
CellChat-style ligand-receptor analysis.
Loads curated mouse LR pairs, computes sender/receiver scores per cluster,
builds communication matrix and exports to CSV.
"""

import argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings('ignore')

# Curated mouse brain LR pairs (ligand → receptor list)
MOUSE_LR = {
    'Tgfb1':  ['Tgfbr1', 'Tgfbr2'],
    'Vegfa':  ['Flt1', 'Kdr'],
    'Bdnf':   ['Ntrk2'],
    'Ntf3':   ['Ntrk3'],
    'Cxcl12': ['Cxcr4', 'Cxcr7'],
    'Ccl2':   ['Ccr2'],
    'Il1b':   ['Il1r1'],
    'Tnf':    ['Tnfrsf1a', 'Tnfrsf1b'],
    'Egf':    ['Egfr'],
    'Hgf':    ['Met'],
    'Sema3a': ['Nrp1'],
    'Nrg1':   ['Erbb4'],
    'Fgf2':   ['Fgfr1'],
    'Pdgfb':  ['Pdgfrb'],
    'Igf1':   ['Igf1r'],
}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--scrna_h5ad',    required=True)
    p.add_argument('--spatial_h5ad',  required=True)
    p.add_argument('--out_networks',  required=True)
    p.add_argument('--out_pathways',  required=True)
    p.add_argument('--out_matrix',    required=True)
    p.add_argument('--out_plots',     required=True)
    return p.parse_args()

def mean_expr(adata, gene, groupby):
    if gene not in adata.var_names:
        return pd.Series(dtype=float)
    X = adata[:, gene].X
    x = X.toarray().flatten() if hasattr(X, 'toarray') else np.array(X).flatten()
    df = pd.DataFrame({'expr': x, 'group': adata.obs[groupby].values})
    return df.groupby('group')['expr'].mean()

def main():
    import scanpy as sc
    args = parse_args()

    adata = sc.read_h5ad(args.scrna_h5ad)
    groupby = 'leiden' if 'leiden' in adata.obs else adata.obs.columns[0]
    cell_types = adata.obs[groupby].unique().tolist()

    rows, comm_mat = [], pd.DataFrame(0.0, index=cell_types, columns=cell_types)
    pathways = set()

    for ligand, receptors in MOUSE_LR.items():
        lig_mean = mean_expr(adata, ligand, groupby)
        for receptor in receptors:
            rec_mean = mean_expr(adata, receptor, groupby)
            for sender in cell_types:
                for receiver in cell_types:
                    ls = float(lig_mean.get(sender, 0))
                    rr = float(rec_mean.get(receiver, 0))
                    strength = ls * rr
                    if strength > 0:
                        rows.append({'sender': sender, 'receiver': receiver,
                                     'ligand': ligand, 'receptor': receptor,
                                     'strength': strength})
                        comm_mat.loc[sender, receiver] += strength
                        pathways.add(ligand[:4])   # simplified pathway name

    net_df = pd.DataFrame(rows)
    net_df.to_csv(args.out_networks, index=False)
    comm_mat.to_csv(args.out_matrix)
    pd.DataFrame({'pathway': sorted(pathways)}).to_csv(args.out_pathways, index=False)

    # Plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    if not comm_mat.empty:
        sns.heatmap(comm_mat, ax=axes[0], cmap='viridis', annot=False)
        axes[0].set_title('Communication matrix')
    if not net_df.empty:
        top = net_df.nlargest(20, 'strength')
        top['pair'] = top['ligand'] + '→' + top['receptor']
        axes[1].barh(top['pair'], top['strength'])
        axes[1].set_title('Top 20 LR pairs')
    plt.tight_layout()
    plt.savefig(args.out_plots, dpi=150); plt.close()
    print(f"[ccc_cellchat] {len(rows)} interactions, {len(pathways)} pathways")

if __name__ == '__main__':
    main()
