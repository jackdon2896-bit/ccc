process SRA_DOWNLOAD {
    tag "${sra_id}"
    label 'process_medium'
    publishDir "${params.outdir}/raw_data/sra/${sra_id}", mode: 'copy'

    input:
    tuple val(sra_id), val(brain_region)

    output:
    tuple val(sra_id), val(brain_region), path("${sra_id}*.fastq.gz"), emit: fastq
    path "${sra_id}_download.log", emit: log

    script:
    """
    echo "Downloading ${sra_id} (${brain_region})" | tee ${sra_id}_download.log

    prefetch ${sra_id} --max-size 50G 2>&1 | tee -a ${sra_id}_download.log

    fasterq-dump ${sra_id} \\
        --split-files \\
        --threads ${task.cpus} \\
        --progress \\
        --temp ./tmp 2>&1 | tee -a ${sra_id}_download.log

    pigz -p ${task.cpus} ${sra_id}*.fastq

    rm -rf ${sra_id} ./tmp

    echo "Download complete: \$(date)" >> ${sra_id}_download.log
    """
}
