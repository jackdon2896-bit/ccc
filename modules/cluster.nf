process CLUSTERING {
    tag "Clustering & UMAP"
    label 'process_medium'
    publishDir "${params.outdir}/clustering", mode: 'copy'
    
    input:
    path adata
    
    output:
    path "clustered_adata.h5ad", emit: adata
    path "umap_clusters.png", emit: umap
    path "cluster_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    # Load data
    adata = sc.read_h5ad("${adata}")
    
    # Feature selection
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=${params.n_neighbors}, n_pcs=30)
    sc.tl.umap(adata)
    
    # Leiden clustering
    sc.tl.leiden(adata, resolution=${params.resolution})
    
    # Save
    adata.write("clustered_adata.h5ad")
    
    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Clusters')
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1], show=False, title='Gene Counts')
    sc.pl.umap(adata, color='total_counts', ax=axes[2], show=False, title='Total Counts')
    
    plt.tight_layout()
    plt.savefig("umap_clusters.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Cluster statistics
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    
    with open("cluster_stats.txt", "w") as f:
        f.write(f"Total clusters: {len(cluster_counts)}\\n")
        f.write(f"\\nCells per cluster:\\n")
        for cluster, count in cluster_counts.items():
            f.write(f"  Cluster {cluster}: {count} cells\\n")
        
        # Find marker genes for each cluster
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
        f.write(f"\\nTop marker genes per cluster:\\n")
        for i in range(min(5, len(cluster_counts))):
            markers = sc.get.rank_genes_groups_df(adata, group=str(i))['names'].head(5).tolist()
            f.write(f"  Cluster {i}: {', '.join(markers)}\\n")
    
    print(f"✅ Clustering complete: {len(cluster_counts)} clusters identified")
    """
}
