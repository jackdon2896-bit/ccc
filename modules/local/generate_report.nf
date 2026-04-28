process GENERATE_REPORT {
    tag "report"
    label 'process_medium'
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path integrated_h5ad
    path plots_dir
    path metrics_json

    output:
    path "report.html",         emit: html
    path "executive_summary.txt", emit: summary

    script:
    """
    python3 ${projectDir}/bin/generate_report.py \\
        --integrated_h5ad  ${integrated_h5ad} \\
        --plots_dir        ${plots_dir} \\
        --metrics_json     ${metrics_json} \\
        --out_html         report.html \\
        --out_summary      executive_summary.txt \\
        --pipeline_version "2.0.0"
    """
}
