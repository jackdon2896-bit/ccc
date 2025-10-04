#!/usr/bin/env nextflow

/*
========================================================================================
    Spatial + scRNA-seq Integration Pipeline - Mode A (Full)
========================================================================================
    Production-ready pipeline for academic research
    Optimized for Seqera Cloud Community Showcase
    Author: Seqera AI
    Version: 1.0.0
========================================================================================
*/

nextflow.enable.dsl = 2

// ========================================================================================
// PARAMETERS
// ========================================================================================

params.scrna_data = null
params.spatial_data = null
params.output_dir = null
params.mode = 'full'

// Integration parameters
params.n_top_genes = 2000
params.n_pcs = 50
params.resolution = 0.5

// Spatial parameters
params.spot_diameter = 100
params.n_neighbors = 15

// Resource parameters
params.max_cpus = 8
params.max_memory = '64.GB'
params.max_time = '6.h'

// ========================================================================================
// INPUT VALIDATION
// ========================================================================================

if (!params.scrna_data) {
    error "ERROR: --scrna_data is required! Please provide path to scRNA-seq data."
}

if (!params.spatial_data) {
    error "ERROR: --spatial_data is required! Please provide path to spatial data."
}

if (!params.output_dir) {
    error "ERROR: --output_dir is required! Please provide output directory path."
}

log.info """\
========================================================================================
    SPATIAL + scRNA-seq INTEGRATION PIPELINE - MODE A
========================================================================================
    Mode              : ${params.mode}
    scRNA-seq data    : ${params.scrna_data}
    Spatial data      : ${params.spatial_data}
    Output directory  : ${params.output_dir}
    
    Integration Settings:
    - Top genes       : ${params.n_top_genes}
    - PCs             : ${params.n_pcs}
    - Resolution      : ${params.resolution}
    
    Spatial Settings:
    - Spot diameter   : ${params.spot_diameter}
    - Neighbors       : ${params.n_neighbors}
    
    Resources:
    - Max CPUs        : ${params.max_cpus}
    - Max memory      : ${params.max_memory}
    - Max time        : ${params.max_time}
========================================================================================
"""

// ========================================================================================
// PROCESS: INTEGRATE SCRNA DATA
// ========================================================================================

process INTEGRATE_SCRNA_DATA {
    tag "scRNA Integration"
    publishDir "${params.output_dir}/integration", mode: 'copy'
    
    cpus 8
    memory '64.GB'
    time '4.h'
    
    container 'quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0'
    
    input:
    path scrna_data
    
    output:
    path 'integrated_image.h5ad', emit: integrated
    path 'integration_qc.png', emit: qc
    path 'integration_metrics.txt', emit: metrics
    path 'versions.txt', emit: versions
    
    script:
    """
    #!/usr/bin/env python3
    
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from datetime import datetime
    
    print("=" * 80)
    print("STARTING scRNA-seq INTEGRATION")
    print("=" * 80)
    print(f"Start time: {datetime.now()}")
    print(f"Input data: ${scrna_data}")
    print(f"Top genes: ${params.n_top_genes}")
    print(f"PCs: ${params.n_pcs}")
    print(f"Resolution: ${params.resolution}")
    print("")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Configure scanpy
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white', figsize=(8, 6))
    
    print("Step 1: Loading scRNA-seq data...")
    # Load data (handle different formats)
    try:
        if '${scrna_data}'.endswith('.h5ad'):
            adata = sc.read_h5ad('${scrna_data}')
        elif '${scrna_data}'.endswith('.h5'):
            adata = sc.read_10x_h5('${scrna_data}')
        else:
            adata = sc.read_10x_mtx('${scrna_data}')
        print(f"Loaded data: {adata.shape[0]} cells x {adata.shape[1]} genes")
    except Exception as e:
        print(f"ERROR loading data: {e}")
        raise
    
    print("\\nStep 2: Quality control filtering...")
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    # Filter cells
    n_cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"Filtered cells: {n_cells_before} -> {adata.n_obs} ({adata.n_obs/n_cells_before*100:.1f}% retained)")
    
    print("\\nStep 3: Normalization and log transformation...")
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print("\\nStep 4: Feature selection...")
    # Identify highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=${params.n_top_genes},
        flavor='seurat_v3',
        batch_key=None
    )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Selected {n_hvg} highly variable genes")
    
    print("\\nStep 5: Scaling and PCA...")
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # Run PCA
    sc.tl.pca(adata, n_comps=${params.n_pcs}, svd_solver='arpack')
    print(f"PCA completed with {${params.n_pcs}} components")
    
    print("\\nStep 6: Computing neighborhood graph...")
    # Compute neighbors
    sc.pp.neighbors(adata, n_pcs=${params.n_pcs}, n_neighbors=${params.n_neighbors})
    
    print("\\nStep 7: UMAP embedding...")
    # Compute UMAP
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)
    
    print("\\nStep 8: Clustering...")
    # Leiden clustering
    sc.tl.leiden(adata, resolution=${params.resolution})
    n_clusters = len(adata.obs['leiden'].unique())
    print(f"Identified {n_clusters} clusters")
    
    print("\\nStep 9: Finding marker genes...")
    # Find marker genes for each cluster
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    print("\\nStep 10: Generating QC plots...")
    # Create comprehensive QC figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # UMAP colored by cluster
    sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False, title='Clusters')
    
    # UMAP colored by n_genes
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[0, 1], show=False, title='N genes')
    
    # UMAP colored by total_counts
    sc.pl.umap(adata, color='total_counts', ax=axes[0, 2], show=False, title='Total counts')
    
    # PCA variance ratio
    axes[1, 0].plot(np.cumsum(adata.uns['pca']['variance_ratio']))
    axes[1, 0].set_xlabel('PC')
    axes[1, 0].set_ylabel('Cumulative variance explained')
    axes[1, 0].set_title('PCA Variance')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Cluster sizes
    cluster_sizes = adata.obs['leiden'].value_counts().sort_index()
    axes[1, 1].bar(range(len(cluster_sizes)), cluster_sizes.values)
    axes[1, 1].set_xlabel('Cluster')
    axes[1, 1].set_ylabel('N cells')
    axes[1, 1].set_title('Cluster Sizes')
    axes[1, 1].grid(True, alpha=0.3, axis='y')
    
    # Gene expression distribution
    axes[1, 2].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black')
    axes[1, 2].set_xlabel('N genes per cell')
    axes[1, 2].set_ylabel('N cells')
    axes[1, 2].set_title('Gene Expression Distribution')
    axes[1, 2].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('integration_qc.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("\\nStep 11: Saving integration metrics...")
    # Save metrics
    with open('integration_metrics.txt', 'w') as f:
        f.write("=" * 80 + "\\n")
        f.write("scRNA-seq INTEGRATION METRICS\\n")
        f.write("=" * 80 + "\\n\\n")
        f.write(f"Analysis date: {datetime.now()}\\n\\n")
        f.write("INPUT DATA:\\n")
        f.write(f"  Original cells: {n_cells_before}\\n")
        f.write(f"  Final cells: {adata.n_obs}\\n")
        f.write(f"  Total genes: {adata.n_vars}\\n\\n")
        f.write("PROCESSING:\\n")
        f.write(f"  Highly variable genes: {n_hvg}\\n")
        f.write(f"  PCs computed: ${params.n_pcs}\\n")
        f.write(f"  Clustering resolution: ${params.resolution}\\n\\n")
        f.write("RESULTS:\\n")
        f.write(f"  Number of clusters: {n_clusters}\\n")
        f.write(f"  Variance explained (top 50 PCs): {np.sum(adata.uns['pca']['variance_ratio'][:50]):.3f}\\n\\n")
        f.write("CLUSTER SIZES:\\n")
        for cluster, size in cluster_sizes.items():
            f.write(f"  Cluster {cluster}: {size} cells ({size/adata.n_obs*100:.1f}%)\\n")
        f.write("\\n" + "=" * 80 + "\\n")
    
    print("\\nStep 12: Saving integrated data...")
    # Save integrated data
    adata.write('integrated_image.h5ad', compression='gzip')
    print(f"Saved integrated data: integrated_image.h5ad ({adata.n_obs} cells x {adata.n_vars} genes)")
    
    # Save versions
    with open('versions.txt', 'w') as f:
        f.write(f"scanpy=={sc.__version__}\\n")
        f.write(f"numpy=={np.__version__}\\n")
        f.write(f"pandas=={pd.__version__}\\n")
    
    print("\\n" + "=" * 80)
    print("INTEGRATION COMPLETED SUCCESSFULLY!")
    print("=" * 80)
    print(f"End time: {datetime.now()}")
    print(f"Output files:")
    print(f"  - integrated_image.h5ad")
    print(f"  - integration_qc.png")
    print(f"  - integration_metrics.txt")
    print("")
    """
}

// ========================================================================================
// PROCESS: SPATIAL ANALYSIS
// ========================================================================================

process SPATIAL_ANALYSIS {
    tag "Spatial Analysis"
    publishDir "${params.output_dir}/spatial", mode: 'copy'
    
    cpus 4
    memory '32.GB'
    time '3.h'
    
    container 'quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0'
    
    input:
    path integrated_data
    path spatial_data
    
    output:
    path 'spatial_analysis_results.h5ad', emit: results
    path 'spatial_plots.png', emit: plots
    path 'spatial_metrics.txt', emit: metrics
    path 'deconvolution_results.csv', emit: deconvolution
    
    script:
    """
    #!/usr/bin/env python3
    
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from datetime import datetime
    
    print("=" * 80)
    print("STARTING SPATIAL ANALYSIS")
    print("=" * 80)
    print(f"Start time: {datetime.now()}")
    print(f"Integrated data: ${integrated_data}")
    print(f"Spatial data: ${spatial_data}")
    print("")
    
    # Set random seed
    np.random.seed(42)
    
    # Configure plotting
    sc.settings.set_figure_params(dpi=150, facecolor='white', figsize=(8, 6))
    
    print("Step 1: Loading integrated scRNA-seq data...")
    adata_ref = sc.read_h5ad('${integrated_data}')
    print(f"Reference data: {adata_ref.shape[0]} cells x {adata_ref.shape[1]} genes")
    
    print("\\nStep 2: Loading spatial data...")
    # Load spatial data (handle different formats)
    try:
        if '${spatial_data}'.endswith('.h5ad'):
            adata_spatial = sc.read_h5ad('${spatial_data}')
        elif '${spatial_data}'.endswith('.h5'):
            adata_spatial = sc.read_visium('${spatial_data}')
        else:
            adata_spatial = sc.read_10x_h5('${spatial_data}')
        print(f"Spatial data: {adata_spatial.shape[0]} spots x {adata_spatial.shape[1]} genes")
    except Exception as e:
        print(f"ERROR loading spatial data: {e}")
        raise
    
    print("\\nStep 3: Preprocessing spatial data...")
    # QC filtering for spatial data
    sc.pp.calculate_qc_metrics(adata_spatial, inplace=True)
    sc.pp.filter_genes(adata_spatial, min_cells=10)
    
    # Normalize
    sc.pp.normalize_total(adata_spatial, target_sum=1e4)
    sc.pp.log1p(adata_spatial)
    
    print("\\nStep 4: Finding common genes...")
    # Get common genes between reference and spatial
    common_genes = adata_ref.var_names.intersection(adata_spatial.var_names)
    print(f"Common genes: {len(common_genes)}")
    
    # Subset to common genes
    adata_ref_subset = adata_ref[:, common_genes].copy()
    adata_spatial_subset = adata_spatial[:, common_genes].copy()
    
    print("\\nStep 5: Computing spatial neighborhoods...")
    # Compute spatial neighbors
    if 'spatial' in adata_spatial.obsm:
        sc.pp.neighbors(adata_spatial_subset, n_neighbors=${params.n_neighbors}, use_rep='spatial')
    else:
        # If no spatial coordinates, use PCA
        sc.pp.highly_variable_genes(adata_spatial_subset, n_top_genes=2000)
        sc.pp.pca(adata_spatial_subset, n_comps=30)
        sc.pp.neighbors(adata_spatial_subset, n_neighbors=${params.n_neighbors})
    
    print("\\nStep 6: Clustering spatial data...")
    # Cluster spatial data
    sc.tl.leiden(adata_spatial_subset, resolution=0.5)
    n_spatial_clusters = len(adata_spatial_subset.obs['leiden'].unique())
    print(f"Identified {n_spatial_clusters} spatial clusters")
    
    print("\\nStep 7: Cell type deconvolution...")
    # Simple deconvolution using correlation
    # Get cluster centroids from reference
    cluster_means = pd.DataFrame()
    for cluster in adata_ref.obs['leiden'].unique():
        cluster_cells = adata_ref[adata_ref.obs['leiden'] == cluster]
        cluster_mean = cluster_cells[:, common_genes].X.mean(axis=0)
        cluster_means[f'cluster_{cluster}'] = np.asarray(cluster_mean).flatten()
    
    cluster_means.index = common_genes
    
    # Correlate each spatial spot with cluster means
    deconv_results = pd.DataFrame(index=adata_spatial_subset.obs_names)
    
    for i, spot in enumerate(adata_spatial_subset.obs_names):
        spot_expr = np.asarray(adata_spatial_subset[i, :].X.toarray()).flatten()
        
        # Compute correlation with each cluster
        correlations = {}
        for cluster_name in cluster_means.columns:
            cluster_expr = cluster_means[cluster_name].values
            corr = np.corrcoef(spot_expr, cluster_expr)[0, 1]
            correlations[cluster_name] = max(0, corr)  # Remove negative correlations
        
        # Normalize to sum to 1
        total = sum(correlations.values())
        if total > 0:
            for k in correlations:
                correlations[k] /= total
        
        for cluster_name, prop in correlations.items():
            deconv_results.loc[spot, cluster_name] = prop
    
    # Add dominant cell type
    adata_spatial_subset.obs['dominant_celltype'] = deconv_results.idxmax(axis=1)
    adata_spatial_subset.obs['celltype_confidence'] = deconv_results.max(axis=1)
    
    print("\\nStep 8: Finding spatially variable genes...")
    # Find spatially variable genes (if spatial coordinates available)
    if 'spatial' in adata_spatial.obsm:
        sc.pp.highly_variable_genes(adata_spatial_subset, flavor='seurat_v3')
    
    print("\\nStep 9: Generating spatial plots...")
    # Create comprehensive spatial analysis figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Spatial clusters
    if 'spatial' in adata_spatial_subset.obsm:
        sc.pl.spatial(adata_spatial_subset, color='leiden', ax=axes[0, 0], show=False, title='Spatial Clusters')
    else:
        axes[0, 0].text(0.5, 0.5, 'No spatial coordinates', ha='center', va='center')
        axes[0, 0].set_title('Spatial Clusters')
    
    # Plot 2: Dominant cell type
    if 'spatial' in adata_spatial_subset.obsm:
        sc.pl.spatial(adata_spatial_subset, color='dominant_celltype', ax=axes[0, 1], show=False, title='Dominant Cell Type')
    else:
        axes[0, 1].text(0.5, 0.5, 'No spatial coordinates', ha='center', va='center')
        axes[0, 1].set_title('Dominant Cell Type')
    
    # Plot 3: Cell type confidence
    if 'spatial' in adata_spatial_subset.obsm:
        sc.pl.spatial(adata_spatial_subset, color='celltype_confidence', ax=axes[0, 2], show=False, title='Confidence', cmap='viridis')
    else:
        axes[0, 2].text(0.5, 0.5, 'No spatial coordinates', ha='center', va='center')
        axes[0, 2].set_title('Confidence')
    
    # Plot 4: Cluster composition heatmap
    composition = pd.crosstab(adata_spatial_subset.obs['leiden'], adata_spatial_subset.obs['dominant_celltype'])
    sns.heatmap(composition, annot=True, fmt='d', cmap='YlOrRd', ax=axes[1, 0])
    axes[1, 0].set_title('Cluster Composition')
    axes[1, 0].set_xlabel('Cell Type')
    axes[1, 0].set_ylabel('Spatial Cluster')
    
    # Plot 5: Deconvolution confidence distribution
    axes[1, 1].hist(adata_spatial_subset.obs['celltype_confidence'], bins=30, edgecolor='black')
    axes[1, 1].set_xlabel('Confidence')
    axes[1, 1].set_ylabel('N spots')
    axes[1, 1].set_title('Deconvolution Confidence')
    axes[1, 1].grid(True, alpha=0.3, axis='y')
    
    # Plot 6: Cell type proportions
    celltype_counts = adata_spatial_subset.obs['dominant_celltype'].value_counts()
    axes[1, 2].bar(range(len(celltype_counts)), celltype_counts.values)
    axes[1, 2].set_xticks(range(len(celltype_counts)))
    axes[1, 2].set_xticklabels(celltype_counts.index, rotation=45, ha='right')
    axes[1, 2].set_ylabel('N spots')
    axes[1, 2].set_title('Cell Type Distribution')
    axes[1, 2].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('spatial_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("\\nStep 10: Saving metrics...")
    # Save metrics
    with open('spatial_metrics.txt', 'w') as f:
        f.write("=" * 80 + "\\n")
        f.write("SPATIAL ANALYSIS METRICS\\n")
        f.write("=" * 80 + "\\n\\n")
        f.write(f"Analysis date: {datetime.now()}\\n\\n")
        f.write("INPUT DATA:\\n")
        f.write(f"  Reference cells: {adata_ref.n_obs}\\n")
        f.write(f"  Spatial spots: {adata_spatial_subset.n_obs}\\n")
        f.write(f"  Common genes: {len(common_genes)}\\n\\n")
        f.write("RESULTS:\\n")
        f.write(f"  Spatial clusters: {n_spatial_clusters}\\n")
        f.write(f"  Cell types detected: {len(celltype_counts)}\\n")
        f.write(f"  Mean confidence: {adata_spatial_subset.obs['celltype_confidence'].mean():.3f}\\n\\n")
        f.write("CELL TYPE DISTRIBUTION:\\n")
        for celltype, count in celltype_counts.items():
            f.write(f"  {celltype}: {count} spots ({count/len(adata_spatial_subset)*100:.1f}%)\\n")
        f.write("\\n" + "=" * 80 + "\\n")
    
    print("\\nStep 11: Saving deconvolution results...")
    # Save deconvolution results
    deconv_results.to_csv('deconvolution_results.csv')
    
    print("\\nStep 12: Saving spatial analysis results...")
    # Save spatial data with annotations
    adata_spatial_subset.write('spatial_analysis_results.h5ad', compression='gzip')
    
    print("\\n" + "=" * 80)
    print("SPATIAL ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)
    print(f"End time: {datetime.now()}")
    print(f"Output files:")
    print(f"  - spatial_analysis_results.h5ad")
    print(f"  - spatial_plots.png")
    print(f"  - spatial_metrics.txt")
    print(f"  - deconvolution_results.csv")
    print("")
    """
}

// ========================================================================================
// WORKFLOW
// ========================================================================================

workflow {
    // Print workflow info
    log.info """
    
    Starting Mode A workflow...
    
    """
    
    // Create input channels
    scrna_ch = channel.fromPath(params.scrna_data)
    spatial_ch = channel.fromPath(params.spatial_data)
    
    // Run integration
    INTEGRATE_SCRNA_DATA(scrna_ch)
    
    // Run spatial analysis
    SPATIAL_ANALYSIS(
        INTEGRATE_SCRNA_DATA.out.integrated,
        spatial_ch
    )
    
    // Workflow completion handler
    workflow.onComplete {
        log.info """
        
        ========================================================================================
        Pipeline completed at: ${workflow.complete}
        Duration            : ${workflow.duration}
        Success             : ${workflow.success}
        Exit status         : ${workflow.exitStatus}
        Error report        : ${workflow.errorReport ?: 'No errors'}
        ========================================================================================
        
        Results saved to: ${params.output_dir}
        
        Output files:
          Integration:
            - integrated_image.h5ad
            - integration_qc.png
            - integration_metrics.txt
          
          Spatial Analysis:
            - spatial_analysis_results.h5ad
            - spatial_plots.png
            - spatial_metrics.txt
            - deconvolution_results.csv
        
        ${workflow.success ? '✓ PIPELINE COMPLETED SUCCESSFULLY!' : '✗ Pipeline failed - check error report'}
        ========================================================================================
        """
    }
    
    workflow.onError {
        log.error """
        ========================================================================================
        ERROR: Pipeline execution failed!
        ========================================================================================
        Error message: ${workflow.errorMessage}
        Error report : ${workflow.errorReport}
        ========================================================================================
        """
    }
}
