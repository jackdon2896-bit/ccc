process ML_FEATURE_VIZ {
    tag "ML Feature Visualization"
    label 'process_low'
    publishDir "${params.outdir}/visualizations", mode: 'copy'
    
    input:
    path ml_results
    
    output:
    path "feature_importance_plot.png", emit: importance
    path "prediction_confidence.png", emit: confidence
    path "cluster_comparison.png", emit: comparison
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    
    # Load results
    adata = sc.read_h5ad("${ml_results}")
    
    # Feature importance plot
    feature_importance = pd.read_csv("${params.outdir}/ml_classification/feature_importance.csv")
    top_features = feature_importance.head(20)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(range(len(top_features)), top_features['importance'])
    ax.set_yticks(range(len(top_features)))
    ax.set_yticklabels(top_features['feature'])
    ax.set_xlabel('Feature Importance')
    ax.set_title('Top 20 Most Important Features for Classification')
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig("feature_importance_plot.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Prediction confidence plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram of confidence scores
    axes[0].hist(adata.obs['ml_confidence'], bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Prediction Confidence')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title('Distribution of ML Prediction Confidence')
    axes[0].axvline(adata.obs['ml_confidence'].mean(), color='red', 
                    linestyle='--', label=f'Mean: {adata.obs["ml_confidence"].mean():.2f}')
    axes[0].legend()
    
    # Confidence by cluster
    confidence_by_cluster = adata.obs.groupby('leiden')['ml_confidence'].mean().sort_values()
    axes[1].barh(range(len(confidence_by_cluster)), confidence_by_cluster.values)
    axes[1].set_yticks(range(len(confidence_by_cluster)))
    axes[1].set_yticklabels([f'Cluster {i}' for i in confidence_by_cluster.index])
    axes[1].set_xlabel('Mean Confidence')
    axes[1].set_title('Average Prediction Confidence by Cluster')
    
    plt.tight_layout()
    plt.savefig("prediction_confidence.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Cluster comparison (true vs predicted)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='True Clusters')
    sc.pl.umap(adata, color='ml_prediction', ax=axes[1], show=False, title='ML Predictions')
    
    plt.tight_layout()
    plt.savefig("cluster_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print("✅ ML feature visualizations complete")
    """
}
