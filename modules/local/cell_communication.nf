/*
    Cell-Cell Communication modules
    CCC_CELLCHAT  — ligand-receptor analysis (CellChat-style)
    CCC_LIANA     — multi-method consensus (LIANA-style)
    CCC_COMMOT    — spatial communication (COMMOT-style)
*/

process CCC_CELLCHAT {
    tag "cellchat"
    label 'process_high'
    publishDir "${params.outdir}/cell_communication/cellchat", mode: 'copy'

    input:
    path scrna_h5ad
    path spatial_h5ad

    output:
    path "cellchat_networks.csv",  emit: networks
    path "cellchat_pathways.csv",  emit: pathways
    path "cellchat_matrix.csv",    emit: matrix
    path "cellchat_plots.png",     emit: plots

    script:
    """
    python3 ${projectDir}/bin/ccc_cellchat.py \\
        --scrna_h5ad      ${scrna_h5ad} \\
        --spatial_h5ad    ${spatial_h5ad} \\
        --out_networks    cellchat_networks.csv \\
        --out_pathways    cellchat_pathways.csv \\
        --out_matrix      cellchat_matrix.csv \\
        --out_plots       cellchat_plots.png
    """
}

process CCC_LIANA {
    tag "liana"
    label 'process_high'
    publishDir "${params.outdir}/cell_communication/liana", mode: 'copy'

    input:
    path scrna_h5ad
    path spatial_h5ad

    output:
    path "liana_results.csv",    emit: interactions
    path "liana_consensus.csv",  emit: consensus
    path "liana_plots.png",      emit: plots

    script:
    """
    python3 ${projectDir}/bin/ccc_liana.py \\
        --scrna_h5ad    ${scrna_h5ad} \\
        --spatial_h5ad  ${spatial_h5ad} \\
        --out_results   liana_results.csv \\
        --out_consensus liana_consensus.csv \\
        --out_plots     liana_plots.png
    """
}

process CCC_COMMOT {
    tag "commot"
    label 'process_high'
    publishDir "${params.outdir}/cell_communication/commot", mode: 'copy'

    input:
    path spatial_h5ad
    path spatial_graph

    output:
    path "commot_results.h5ad",  emit: results
    path "commot_flows.csv",     emit: flows
    path "commot_plots.png",     emit: plots

    script:
    """
    python3 ${projectDir}/bin/ccc_commot.py \\
        --spatial_h5ad  ${spatial_h5ad} \\
        --spatial_graph ${spatial_graph} \\
        --out_h5ad      commot_results.h5ad \\
        --out_flows     commot_flows.csv \\
        --out_plots     commot_plots.png \\
        --distance_thr  ${params.ccc_distance_threshold}
    """
}
