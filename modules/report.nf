process FINAL_REPORT {
    tag "Final Report"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path integrated_data
    path communication_results
    path ml_results
    
    output:
    path "PIPELINE_REPORT.html", emit: html
    path "summary_statistics.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    from datetime import datetime
    
    # Load all results
    integrated = sc.read_h5ad("${integrated_data}")
    comm = sc.read_h5ad("${communication_results}")
    ml = sc.read_h5ad("${ml_results}")
    
    # Generate HTML report
    html_content = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <title>Cell-Cell Communication Pipeline Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 40px;
                background-color: #f5f5f5;
            }}
            .container {{
                background-color: white;
                padding: 30px;
                border-radius: 10px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            }}
            h1 {{
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }}
            h2 {{
                color: #34495e;
                margin-top: 30px;
            }}
            .metric {{
                display: inline-block;
                background-color: #ecf0f1;
                padding: 15px 25px;
                margin: 10px;
                border-radius: 5px;
                border-left: 4px solid #3498db;
            }}
            .metric-value {{
                font-size: 24px;
                font-weight: bold;
                color: #2c3e50;
            }}
            .metric-label {{
                font-size: 12px;
                color: #7f8c8d;
                text-transform: uppercase;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
            }}
            th, td {{
                border: 1px solid #ddd;
                padding: 12px;
                text-align: left;
            }}
            th {{
                background-color: #3498db;
                color: white;
            }}
            tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            .success {{
                color: #27ae60;
                font-weight: bold;
            }}
            .footer {{
                margin-top: 40px;
                padding-top: 20px;
                border-top: 1px solid #ddd;
                color: #7f8c8d;
                font-size: 12px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧬 Cell-Cell Communication Spatial Pipeline Report</h1>
            <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            
            <h2>📊 Summary Statistics</h2>
            <div>
                <div class="metric">
                    <div class="metric-value">{integrated.n_obs}</div>
                    <div class="metric-label">Total Cells</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{integrated.n_vars}</div>
                    <div class="metric-label">Genes Analyzed</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{len(integrated.obs['leiden'].unique())}</div>
                    <div class="metric-label">Clusters Identified</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{len(comm.uns['liana_res'])}</div>
                    <div class="metric-label">Interactions Tested</div>
                </div>
            </div>
            
            <h2>🔬 Analysis Modules Completed</h2>
            <table>
                <tr>
                    <th>Module</th>
                    <th>Status</th>
                    <th>Key Output</th>
                </tr>
                <tr>
                    <td>Image Preprocessing</td>
                    <td class="success">✓ Complete</td>
                    <td>Normalized TIFF image</td>
                </tr>
                <tr>
                    <td>Cell Segmentation</td>
                    <td class="success">✓ Complete</td>
                    <td>Cellpose masks generated</td>
                </tr>
                <tr>
                    <td>Quality Control</td>
                    <td class="success">✓ Complete</td>
                    <td>{integrated.n_obs} cells passed filters</td>
                </tr>
                <tr>
                    <td>Clustering</td>
                    <td class="success">✓ Complete</td>
                    <td>{len(integrated.obs['leiden'].unique())} clusters</td>
                </tr>
                <tr>
                    <td>Cell Communication (LIANA)</td>
                    <td class="success">✓ Complete</td>
                    <td>{len(comm.uns['liana_res'])} interactions</td>
                </tr>
                <tr>
                    <td>ML Classification</td>
                    <td class="success">✓ Complete</td>
                    <td>Random Forest trained</td>
                </tr>
            </table>
            
            <h2>🎯 Top Cell-Cell Interactions</h2>
            <table>
                <tr>
                    <th>Source → Target</th>
                    <th>Ligand - Receptor</th>
                    <th>Significance Rank</th>
                </tr>
    '''
    
    # Add top interactions
    liana_res = comm.uns['liana_res']
    top_interactions = liana_res.nsmallest(10, 'magnitude_rank')
    
    for _, row in top_interactions.iterrows():
        html_content += f'''
                <tr>
                    <td>{row['source']} → {row['target']}</td>
                    <td>{row['ligand']} - {row['receptor']}</td>
                    <td>{row['magnitude_rank']:.4f}</td>
                </tr>
        '''
    
    html_content += f'''
            </table>
            
            <h2>📈 ML Classification Results</h2>
            <p><strong>Model:</strong> Random Forest Classifier</p>
            <p><strong>Mean Prediction Confidence:</strong> {ml.obs['ml_confidence'].mean():.3f}</p>
            <p><strong>Features Used:</strong> Highly variable genes</p>
            
            <h2>📁 Output Files</h2>
            <ul>
                <li><strong>clustering/</strong> - UMAP plots and cluster assignments</li>
                <li><strong>communication/</strong> - Cell-cell interaction networks</li>
                <li><strong>ml_classification/</strong> - ML model predictions and feature importance</li>
                <li><strong>visualizations/</strong> - All visualization plots</li>
                <li><strong>integration/</strong> - Integrated spatial + scRNA-seq data</li>
            </ul>
            
            <div class="footer">
                <p>Pipeline: Cell-Cell Communication Spatial Analysis v1.0.0</p>
                <p>Generated by Seqera AI | Nextflow DSL2</p>
            </div>
        </div>
    </body>
    </html>
    '''
    
    # Save HTML
    with open("PIPELINE_REPORT.html", "w") as f:
        f.write(html_content)
    
    # Generate text summary
    with open("summary_statistics.txt", "w") as f:
        f.write("="*70 + "\\n")
        f.write("CELL-CELL COMMUNICATION SPATIAL PIPELINE - FINAL SUMMARY\\n")
        f.write("="*70 + "\\n\\n")
        
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n\\n")
        
        f.write("DATASET STATISTICS:\\n")
        f.write(f"  Total cells: {integrated.n_obs}\\n")
        f.write(f"  Total genes: {integrated.n_vars}\\n")
        f.write(f"  Clusters identified: {len(integrated.obs['leiden'].unique())}\\n")
        f.write(f"  Mean genes/cell: {integrated.obs['n_genes_by_counts'].mean():.1f}\\n\\n")
        
        f.write("CELL-CELL COMMUNICATION:\\n")
        f.write(f"  Interactions tested: {len(liana_res)}\\n")
        significant = liana_res[liana_res['magnitude_rank'] < 0.05]
        f.write(f"  Significant interactions: {len(significant)}\\n")
        f.write(f"  Unique ligands: {liana_res['ligand'].nunique()}\\n")
        f.write(f"  Unique receptors: {liana_res['receptor'].nunique()}\\n\\n")
        
        f.write("MACHINE LEARNING:\\n")
        f.write(f"  Model: Random Forest Classifier\\n")
        f.write(f"  Mean prediction confidence: {ml.obs['ml_confidence'].mean():.3f}\\n")
        f.write(f"  Min confidence: {ml.obs['ml_confidence'].min():.3f}\\n")
        f.write(f"  Max confidence: {ml.obs['ml_confidence'].max():.3f}\\n\\n")
        
        f.write("OUTPUT DIRECTORIES:\\n")
        f.write(f"  Main results: ${params.outdir}/\\n")
        f.write(f"  Report location: ${params.outdir}/PIPELINE_REPORT.html\\n")
        
        f.write("\\n" + "="*70 + "\\n")
        f.write("Pipeline execution completed successfully!\\n")
        f.write("="*70 + "\\n")
    
    print("✅ Final report generated")
    """
}
