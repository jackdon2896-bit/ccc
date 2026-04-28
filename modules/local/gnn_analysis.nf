process GNN_ANALYSIS {
    tag "gnn"
    label 'process_high'
    publishDir "${params.outdir}/machine_learning/gnn", mode: 'copy'

    input:
    path spatial_h5ad
    path spatial_graph

    output:
    path "gnn_results.h5ad",   emit: results
    path "gnn_embeddings.csv", emit: embeddings
    path "gnn_plots.png",      emit: plots
    path "gnn_metrics.json",   emit: metrics

    script:
    """
    python3 ${projectDir}/bin/gnn_analysis.py \\
        --spatial_h5ad  ${spatial_h5ad} \\
        --spatial_graph ${spatial_graph} \\
        --out_h5ad      gnn_results.h5ad \\
        --out_embeddings gnn_embeddings.csv \\
        --out_plots     gnn_plots.png \\
        --out_metrics   gnn_metrics.json \\
        --n_epochs      ${params.n_epochs} \\
        --lr            ${params.learning_rate} \\
        --hidden_dim    ${params.hidden_dim}
    """
}
