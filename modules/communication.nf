process CELL_COMMUNICATION {
    tag "Cell-Cell Communication"
    label 'process_high'
    publishDir "${params.outdir}/communication", mode: 'copy'
    
    input:
    path adata
    
    output:
    path "communication_results.h5ad", emit: results
    path "ligand_receptor_network.csv", emit: network
    path "communication_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import liana as li
    import pandas as pd
    
    # Load data
    adata = sc.read_h5ad("${adata}")
    
    # Run LIANA for cell-cell communication
    print("Running LIANA cell-cell communication analysis...")
    
    # Ensure we have cluster information
    if 'leiden' not in adata.obs.columns:
        raise ValueError("Clustering must be performed before communication analysis")
    
    # Run LIANA with CellPhoneDB resource
    li.mt.rank_aggregate(
        adata,
        groupby='leiden',
        resource_name='cellphonedb',
        expr_prop=0.1,
        min_cells=5,
        n_perms=100,
        seed=42,
        verbose=True,
        use_raw=False
    )
    
    # Extract results
    liana_res = adata.uns['liana_res']
    
    # Filter significant interactions
    significant = liana_res[liana_res['magnitude_rank'] < 0.05]
    
    # Save results
    adata.write("communication_results.h5ad")
    significant.to_csv("ligand_receptor_network.csv", index=False)
    
    # Summary statistics
    with open("communication_summary.txt", "w") as f:
        f.write(f"Cell-Cell Communication Analysis Results\\n")
        f.write(f"{'='*50}\\n\\n")
        f.write(f"Total interactions tested: {len(liana_res)}\\n")
        f.write(f"Significant interactions: {len(significant)}\\n")
        f.write(f"Unique ligands: {liana_res['ligand'].nunique()}\\n")
        f.write(f"Unique receptors: {liana_res['receptor'].nunique()}\\n")
        f.write(f"Cluster pairs analyzed: {liana_res[['source', 'target']].drop_duplicates().shape[0]}\\n\\n")
        
        f.write(f"Top 10 interactions by magnitude:\\n")
        top_interactions = significant.nsmallest(10, 'magnitude_rank')
        for idx, row in top_interactions.iterrows():
            f.write(f"  {row['source']} -> {row['target']}: {row['ligand']} - {row['receptor']} (rank: {row['magnitude_rank']:.3f})\\n")
    
    print(f"✅ Communication analysis complete: {len(significant)} significant interactions")
    """
}
