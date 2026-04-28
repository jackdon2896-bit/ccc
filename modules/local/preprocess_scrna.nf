process PREPROCESS_SCRNA {
    tag "${sra_id}"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessed/scrna/${sra_id}", mode: 'copy'

    input:
    tuple val(sra_id), val(brain_region), path(fastq_files)

    output:
    tuple val(sra_id), val(brain_region), path("${sra_id}_processed.h5ad"), emit: h5ad
    path "${sra_id}_qc.json", emit: qc_metrics

    script:
    """
    python3 ${projectDir}/bin/preprocess_scrna.py \\
        --sra_id       ${sra_id} \\
        --brain_region ${brain_region} \\
        --fastq        ${fastq_files} \\
        --out_h5ad     ${sra_id}_processed.h5ad \\
        --out_qc       ${sra_id}_qc.json \\
        --min_genes    ${params.min_genes} \\
        --min_cells    ${params.min_cells} \\
        --max_genes    ${params.max_genes} \\
        --mt_max       ${params.mt_percent_max}
    """
}
