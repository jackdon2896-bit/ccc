process COMMUNICATION_VIZ {
    tag "Communication Visualization"
    label 'process_medium'
    publishDir "${params.outdir}/visualizations", mode: 'copy'
    
    input:
    path communication_results
    
    output:
    path "network_plot.png", emit: network
    path "dotplot_interactions.png", emit: dotplot
    path "chord_diagram.png", emit: chord
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import liana as li
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    # Load results
    adata = sc.read_h5ad("${communication_results}")
    liana_res = adata.uns['liana_res']
    
    # Filter significant interactions
    significant = liana_res[liana_res['magnitude_rank'] < 0.05]
    
    if len(significant) == 0:
        print("Warning: No significant interactions found, using top 50 by rank")
        significant = liana_res.nsmallest(50, 'magnitude_rank')
    
    # Network plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create adjacency matrix
    clusters = sorted(set(significant['source'].unique()) | set(significant['target'].unique()))
    adj_matrix = pd.DataFrame(0, index=clusters, columns=clusters)
    
    for _, row in significant.iterrows():
        adj_matrix.loc[row['source'], row['target']] += 1
    
    # Plot heatmap
    import seaborn as sns
    sns.heatmap(adj_matrix, annot=True, fmt='d', cmap='YlOrRd', 
                square=True, linewidths=0.5, ax=ax)
    ax.set_title('Cell-Cell Communication Network\\n(Number of Interactions)')
    ax.set_xlabel('Target Cluster')
    ax.set_ylabel('Source Cluster')
    plt.tight_layout()
    plt.savefig("network_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Dotplot of top interactions
    top_interactions = significant.nsmallest(30, 'magnitude_rank')
    
    fig, ax = plt.subplots(figsize=(10, 12))
    
    # Create labels
    top_interactions['interaction'] = (
        top_interactions['ligand'] + ' - ' + top_interactions['receptor']
    )
    top_interactions['pair'] = (
        top_interactions['source'] + ' → ' + top_interactions['target']
    )
    
    # Plot
    scatter = ax.scatter(
        range(len(top_interactions)),
        top_interactions['pair'],
        s=200,
        c=top_interactions['magnitude_rank'],
        cmap='viridis_r',
        alpha=0.8,
        edgecolors='black',
        linewidth=0.5
    )
    
    ax.set_yticks(range(len(top_interactions)))
    ax.set_yticklabels(top_interactions['interaction'], fontsize=8)
    ax.set_xlabel('Interaction Index')
    ax.set_title('Top 30 Ligand-Receptor Interactions')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Magnitude Rank (lower = stronger)')
    
    plt.tight_layout()
    plt.savefig("dotplot_interactions.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Chord diagram (simplified circular plot)
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))
    
    # Count interactions per cluster pair
    interaction_counts = significant.groupby(['source', 'target']).size().reset_index(name='count')
    
    # Get unique clusters
    all_clusters = sorted(set(interaction_counts['source']) | set(interaction_counts['target']))
    n_clusters = len(all_clusters)
    
    # Assign angles
    angles = np.linspace(0, 2 * np.pi, n_clusters, endpoint=False)
    cluster_angles = dict(zip(all_clusters, angles))
    
    # Plot arcs for interactions
    for _, row in interaction_counts.iterrows():
        source_angle = cluster_angles[row['source']]
        target_angle = cluster_angles[row['target']]
        
        # Draw arc
        theta = np.linspace(source_angle, target_angle, 50)
        r = np.ones(50) * row['count'] / interaction_counts['count'].max()
        ax.plot(theta, r, alpha=0.3, linewidth=2)
    
    # Plot cluster positions
    ax.scatter(angles, np.ones(n_clusters), s=500, c=range(n_clusters), cmap='tab20')
    
    # Add labels
    for cluster, angle in cluster_angles.items():
        ax.text(angle, 1.15, f'C{cluster}', ha='center', va='center', fontsize=10)
    
    ax.set_ylim(0, 1.3)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title('Communication Chord Diagram', y=1.08, fontsize=14)
    
    plt.tight_layout()
    plt.savefig("chord_diagram.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Generated {len(significant)} interaction visualizations")
    """
}
