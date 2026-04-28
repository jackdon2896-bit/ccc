process FINAL_INTEGRATION {
    tag "integration"
    label 'process_high'
    publishDir "${params.outdir}/final_integration", mode: 'copy'

    input:
    path result_files   // collected list: [type, path] tuples as flat files
    path original_image

    output:
    path "integrated_data.h5ad",       emit: data
    path "integration_plots/",         emit: plots
    path "integration_metrics.json",   emit: metrics

    script:
    """
    mkdir -p integration_plots

    python3 ${projectDir}/bin/final_integration.py \\
        --result_files   ${result_files} \\
        --original_image ${original_image} \\
        --out_h5ad       integrated_data.h5ad \\
        --out_plots_dir  integration_plots \\
        --out_metrics    integration_metrics.json \\
        --genome         ${params.genome}
    """
}
