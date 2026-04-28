process SCRNA_ANALYSIS {
    tag "scrna_analysis"
    label 'process_high'
    publishDir "${params.outdir}/scrna_analysis", mode: 'copy'

    input:
    path h5ad_files   // collected list from all brain regions

    output:
    path "scrna_integrated.h5ad", emit: processed
    path "scrna_plots.png",       emit: plots
    path "scrna_metrics.json",    emit: metrics
    path "marker_genes.csv",      emit: markers

    script:
    """
    python3 ${projectDir}/bin/scrna_analysis.py \\
        --h5ads       ${h5ad_files} \\
        --out_h5ad    scrna_integrated.h5ad \\
        --out_plots   scrna_plots.png \\
        --out_metrics scrna_metrics.json \\
        --out_markers marker_genes.csv \\
        --n_top_genes ${params.n_top_genes} \\
        --n_pcs       ${params.n_pcs} \\
        --resolution  ${params.clustering_resolution}
    """
}
