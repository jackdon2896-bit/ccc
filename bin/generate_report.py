#!/usr/bin/env python3
"""Generate HTML report from integrated AnnData + metrics JSON."""

import argparse, json, warnings
from pathlib import Path
from datetime import datetime
warnings.filterwarnings('ignore')

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--integrated_h5ad',  required=True)
    p.add_argument('--plots_dir',        required=True)
    p.add_argument('--metrics_json',     required=True)
    p.add_argument('--out_html',         required=True)
    p.add_argument('--out_summary',      required=True)
    p.add_argument('--pipeline_version', default='2.0.0')
    return p.parse_args()

def main():
    import scanpy as sc
    args = parse_args()
    adata = sc.read_h5ad(args.integrated_h5ad)
    metrics = json.loads(Path(args.metrics_json).read_text())

    # Gather plot PNGs
    plots = sorted(Path(args.plots_dir).glob('*.png'))

    # Build HTML
    imgs_html = '\n'.join(
        f'<div class="plot"><h3>{p.stem}</h3>'
        f'<img src="{p}" style="max-width:100%;border:1px solid #ccc"></div>'
        for p in plots
    )

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Golden Standard Pipeline Report v{args.pipeline_version}</title>
  <style>
    body {{font-family: Arial, sans-serif; margin: 40px; background:#f9f9f9}}
    h1 {{color:#2c3e50}} h2 {{color:#16a085}} h3 {{color:#555}}
    table {{border-collapse:collapse; width:100%; margin-bottom:20px}}
    th,td {{border:1px solid #ddd; padding:8px; text-align:left}}
    th {{background-color:#16a085; color:white}}
    tr:nth-child(even){{background:#f2f2f2}}
    .section {{background:white; padding:20px; margin-bottom:20px;
               border-radius:8px; box-shadow:0 2px 4px rgba(0,0,0,0.1)}}
    .plot {{margin-bottom:30px}}
    .badge {{display:inline-block;padding:4px 10px;border-radius:12px;
             font-size:12px;font-weight:bold}}
    .yes {{background:#27ae60;color:white}} .no {{background:#e74c3c;color:white}}
  </style>
</head>
<body>
<h1>🧬 Golden Standard Spatial + scRNA-seq + ML + CCC Pipeline</h1>
<p><b>Version:</b> {args.pipeline_version} &nbsp; <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

<div class="section">
<h2>📊 Integration Metrics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Spatial spots</td><td>{metrics.get('n_spots','N/A')}</td></tr>
<tr><td>Genes</td><td>{metrics.get('n_genes','N/A')}</td></tr>
<tr><td>Reference genome</td><td>{metrics.get('genome','mm10')}</td></tr>
<tr><td>scRNA-seq mapping</td>
    <td><span class="badge {'yes' if metrics.get('has_scrna_mapping') else 'no'}">{'✓' if metrics.get('has_scrna_mapping') else '✗'}</span></td></tr>
<tr><td>ML predictions</td>
    <td><span class="badge {'yes' if metrics.get('has_ml_predictions') else 'no'}">{'✓' if metrics.get('has_ml_predictions') else '✗'}</span></td></tr>
<tr><td>GNN embeddings</td>
    <td><span class="badge {'yes' if metrics.get('has_gnn_embeddings') else 'no'}">{'✓' if metrics.get('has_gnn_embeddings') else '✗'}</span></td></tr>
<tr><td>Pseudotime</td>
    <td><span class="badge {'yes' if metrics.get('has_pseudotime') else 'no'}">{'✓' if metrics.get('has_pseudotime') else '✗'}</span></td></tr>
<tr><td>CCC methods integrated</td><td>{metrics.get('ccc_methods',0)}</td></tr>
</table>
</div>

<div class="section">
<h2>🔬 Analysis Layers</h2>
<table>
<tr><th>Layer</th><th>Tool / Method</th><th>Output</th></tr>
<tr><td>Spatial QC</td><td>Scanpy + Squidpy</td><td>spatial_processed.h5ad</td></tr>
<tr><td>scRNA-seq integration</td><td>Scanpy + Harmony</td><td>scrna_integrated.h5ad</td></tr>
<tr><td>Image segmentation</td><td>Cellpose</td><td>segmentation_mask.tif</td></tr>
<tr><td>Spatial analysis</td><td>Squidpy (UMAP, Leiden, Moran's I)</td><td>spatial_analyzed.h5ad</td></tr>
<tr><td>ML prediction</td><td>Random Forest (sklearn)</td><td>ml_predictions.csv</td></tr>
<tr><td>GNN</td><td>GCNConv (PyTorch Geometric)</td><td>gnn_embeddings.csv</td></tr>
<tr><td>Trajectory</td><td>PAGA + DPT (Scanpy)</td><td>pseudotime.csv</td></tr>
<tr><td>CCC — CellChat-style</td><td>Ligand-receptor scoring (Python)</td><td>cellchat_networks.csv</td></tr>
<tr><td>CCC — LIANA-style</td><td>Multi-method consensus</td><td>liana_consensus.csv</td></tr>
<tr><td>CCC — COMMOT-style</td><td>Spatial communication flows</td><td>commot_flows.csv</td></tr>
<tr><td>Multi-modal integration</td><td>Correlation transfer + OT</td><td>integrated_data.h5ad</td></tr>
</table>
</div>

<div class="section">
<h2>🖼️ Result Plots</h2>
{imgs_html if imgs_html else '<p>No plots found in output directory.</p>'}
</div>
</body>
</html>"""

    Path(args.out_html).write_text(html)

    summary = (
        f"Golden Standard Pipeline Report — v{args.pipeline_version}\n"
        f"Generated: {datetime.now()}\n"
        f"{'='*60}\n"
        f"Spatial spots  : {metrics.get('n_spots','N/A')}\n"
        f"Genes          : {metrics.get('n_genes','N/A')}\n"
        f"Genome         : {metrics.get('genome','mm10')}\n"
        f"scRNA mapping  : {'Yes' if metrics.get('has_scrna_mapping') else 'No'}\n"
        f"ML predictions : {'Yes' if metrics.get('has_ml_predictions') else 'No'}\n"
        f"GNN embeddings : {'Yes' if metrics.get('has_gnn_embeddings') else 'No'}\n"
        f"Pseudotime     : {'Yes' if metrics.get('has_pseudotime') else 'No'}\n"
        f"CCC methods    : {metrics.get('ccc_methods',0)}\n"
    )
    Path(args.out_summary).write_text(summary)
    print("[generate_report] Done.")

if __name__ == '__main__':
    main()
