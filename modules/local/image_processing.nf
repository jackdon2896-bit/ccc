process CELLPOSE_SEGMENT {
    tag "cellpose"
    label 'process_high'
    publishDir "${params.outdir}/image_processing", mode: 'copy'

    input:
    path image_file

    output:
    path "segmentation_mask.tif",    emit: mask
    path "segmentation_overlay.png", emit: overlay
    path "cell_boundaries.csv",      emit: boundaries
    path "segmentation_stats.json",  emit: stats

    script:
    """
    python3 ${projectDir}/bin/cellpose_segment.py \\
        --image         ${image_file} \\
        --out_mask      segmentation_mask.tif \\
        --out_overlay   segmentation_overlay.png \\
        --out_bounds    cell_boundaries.csv \\
        --out_stats     segmentation_stats.json
    """
}
