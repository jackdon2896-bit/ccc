process SPATIAL_QC {
    tag "Quality Control"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    path h5_file
    path roi_coords
    
    output:
    path "spatial_adata.h5ad", emit: adata
    path "qc_metrics.png", emit: qc_plot
    path "qc_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Load spatial data
    adata = sc.read_10x_h5("${h5_file}")
    
    # Load ROI coordinates
    coords = pd.read_csv("${roi_coords}")
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filtering
    adata = adata[adata.obs.n_genes_by_counts > ${params.min_genes}, :]
    adata = adata[:, adata.var.n_cells_by_counts > ${params.min_cells}]
    adata = adata[adata.obs.pct_counts_mt < ${params.max_mito_pct}, :]
    
    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Store raw counts
    adata.raw = adata
    
    # Save
    adata.write("spatial_adata.h5ad")
    
    # QC plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=False, ax=axes[0, 0], show=False)
    
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[0, 1], show=False)
    
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1, 0], show=False)
    
    axes[1, 1].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black')
    axes[1, 1].set_xlabel('Number of Genes')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Gene Count Distribution')
    
    plt.tight_layout()
    plt.savefig("qc_metrics.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Summary
    with open("qc_summary.txt", "w") as f:
        f.write(f"Cells after QC: {adata.n_obs}\\n")
        f.write(f"Genes after QC: {adata.n_vars}\\n")
        f.write(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.1f}\\n")
        f.write(f"Mean counts per cell: {adata.obs['total_counts'].mean():.1f}\\n")
        f.write(f"Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}\\n")
    
    print(f"✅ QC complete: {adata.n_obs} cells, {adata.n_vars} genes")
    """
}
