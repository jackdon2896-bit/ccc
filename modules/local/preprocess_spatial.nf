process PREPROCESS_SPATIAL {
    tag "spatial"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessed/spatial", mode: 'copy'

    input:
    path h5_file
    path image_file

    output:
    path "spatial_processed.h5ad", emit: h5ad
    path "spatial_qc.png",         emit: qc_plot
    path "spatial_image_proc.tif", emit: image
    path "spatial_coords.csv",     emit: coords
    path "spatial_qc.json",        emit: qc_metrics

    script:
    """
    python3 ${projectDir}/bin/preprocess_spatial.py \\
        --h5         ${h5_file} \\
        --image      ${image_file} \\
        --out_h5ad   spatial_processed.h5ad \\
        --out_image  spatial_image_proc.tif \\
        --out_coords spatial_coords.csv \\
        --out_qc     spatial_qc.json \\
        --out_plot   spatial_qc.png \\
        --min_genes  ${params.min_genes} \\
        --min_cells  ${params.min_cells} \\
        --max_genes  ${params.max_genes} \\
        --mt_max     ${params.mt_percent_max}
    """
}
