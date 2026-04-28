process TRAJECTORY_INFERENCE {
    tag "trajectory"
    label 'process_high'
    publishDir "${params.outdir}/machine_learning/trajectory", mode: 'copy'

    input:
    path spatial_h5ad
    path scrna_h5ad

    output:
    path "trajectory_results.h5ad", emit: results
    path "pseudotime.csv",          emit: pseudotime
    path "trajectory_plots.png",    emit: plots

    script:
    """
    python3 ${projectDir}/bin/trajectory_inference.py \\
        --spatial_h5ad  ${spatial_h5ad} \\
        --scrna_h5ad    ${scrna_h5ad} \\
        --out_h5ad      trajectory_results.h5ad \\
        --out_pseudotime pseudotime.csv \\
        --out_plots     trajectory_plots.png \\
        --n_pcs         ${params.n_pcs}
    """
}
