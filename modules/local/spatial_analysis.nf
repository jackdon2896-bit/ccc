process SPATIAL_ANALYSIS {
    tag "spatial_analysis"
    label 'process_high'
    publishDir "${params.outdir}/spatial_analysis", mode: 'copy'

    input:
    path h5ad_file
    path mask_file

    output:
    path "spatial_analyzed.h5ad",  emit: processed
    path "spatial_graph.csv",      emit: spatial_graph
    path "spatial_plots.png",      emit: plots
    path "spatial_metrics.json",   emit: metrics

    script:
    """
    python3 ${projectDir}/bin/spatial_analysis.py \\
        --h5ad           ${h5ad_file} \\
        --mask           ${mask_file} \\
        --out_h5ad       spatial_analyzed.h5ad \\
        --out_graph      spatial_graph.csv \\
        --out_plots      spatial_plots.png \\
        --out_metrics    spatial_metrics.json \\
        --n_top_genes    ${params.n_top_genes} \\
        --n_pcs          ${params.n_pcs} \\
        --resolution     ${params.clustering_resolution} \\
        --n_neighbors    ${params.n_neighbors}
    """
}
