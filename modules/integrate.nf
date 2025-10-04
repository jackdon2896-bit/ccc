process SCRNA_INTEGRATION {
    tag "scRNA-seq Integration"
    label 'process_high'
    publishDir "${params.outdir}/integration", mode: 'copy'
    
    input:
    path spatial_adata
    path scrna_ref
    
    output:
    path "integrated_adata.h5ad", emit: integrated
    path "integration_plot.png", emit: plot
    path "integration_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Load data
    spatial = sc.read_h5ad("${spatial_adata}")
    reference = sc.read_h5ad("${scrna_ref}")
    
    print(f"Spatial data: {spatial.n_obs} cells, {spatial.n_vars} genes")
    print(f"Reference data: {reference.n_obs} cells, {reference.n_vars} genes")
    
    # Find common genes
    common_genes = spatial.var_names.intersection(reference.var_names)
    print(f"Common genes: {len(common_genes)}")
    
    if len(common_genes) < 100:
        raise ValueError(f"Too few common genes ({len(common_genes)}) for integration")
    
    # Subset to common genes
    spatial_subset = spatial[:, common_genes].copy()
    reference_subset = reference[:, common_genes].copy()
    
    # Label datasets
    spatial_subset.obs['dataset'] = 'spatial'
    reference_subset.obs['dataset'] = 'reference'
    
    # Concatenate
    combined = sc.concat([spatial_subset, reference_subset], join='inner')
    
    # Normalize and scale
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    sc.pp.scale(combined, max_value=10)
    
    # PCA
    sc.tl.pca(combined, svd_solver='arpack')
    
    # Batch correction with Harmony
    try:
        import harmonypy
        ho = harmonypy.run_harmony(
            combined.obsm['X_pca'],
            combined.obs,
            'dataset',
            max_iter_harmony=20
        )
        combined.obsm['X_pca_harmony'] = ho.Z_corr.T
        
        # Use harmony-corrected PCs for UMAP
        sc.pp.neighbors(combined, use_rep='X_pca_harmony')
        sc.tl.umap(combined)
        
        integration_method = "Harmony"
    except ImportError:
        print("Harmony not available, using standard integration")
        sc.pp.neighbors(combined)
        sc.tl.umap(combined)
        integration_method = "Standard PCA"
    
    # Transfer annotations
    spatial_integrated = combined[combined.obs['dataset'] == 'spatial'].copy()
    
    # Save
    spatial_integrated.write("integrated_adata.h5ad")
    
    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sc.pl.umap(combined, color='dataset', ax=axes[0], show=False, title='Dataset')
    
    if 'leiden' in combined.obs.columns:
        sc.pl.umap(combined, color='leiden', ax=axes[1], show=False, title='Clusters')
    
    # Show integration quality
    spatial_cells = combined.obs['dataset'] == 'spatial'
    axes[2].scatter(
        combined.obsm['X_umap'][~spatial_cells, 0],
        combined.obsm['X_umap'][~spatial_cells, 1],
        c='lightgray', s=1, alpha=0.3, label='Reference'
    )
    axes[2].scatter(
        combined.obsm['X_umap'][spatial_cells, 0],
        combined.obsm['X_umap'][spatial_cells, 1],
        c='red', s=5, alpha=0.6, label='Spatial'
    )
    axes[2].set_title('Integration Quality')
    axes[2].legend()
    axes[2].axis('off')
    
    plt.tight_layout()
    plt.savefig("integration_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Summary
    with open("integration_summary.txt", "w") as f:
        f.write(f"scRNA-seq Integration Summary\\n")
        f.write(f"{'='*50}\\n\\n")
        f.write(f"Integration method: {integration_method}\\n")
        f.write(f"Spatial cells: {spatial.n_obs}\\n")
        f.write(f"Reference cells: {reference.n_obs}\\n")
        f.write(f"Common genes: {len(common_genes)}\\n")
        f.write(f"Final integrated cells: {spatial_integrated.n_obs}\\n")
    
    print(f"✅ Integration complete: {spatial_integrated.n_obs} cells")
    """
}
