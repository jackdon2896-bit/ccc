#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    Cell-Cell Communication Spatial Transcriptomics Pipeline
========================================================================================
    Advanced spatial omics pipeline with:
    - Cell-cell communication analysis (LIANA/COMMOT)
    - Machine learning tissue classification
    - Multi-modal integration (spatial + scRNA-seq)
    - AWS Batch deployment ready
    
    Author: Seqera AI
    Version: 1.0.0
========================================================================================
*/

// Parameters
params.tiff = null
params.h5 = null
params.sra_ids = null
params.scrna_ref = null
params.outdir = "results"
params.min_genes = 200
params.min_cells = 3
params.max_mito_pct = 20

// Validation
if (!params.tiff) {
    error "ERROR: --tiff parameter is required (spatial TIFF image)"
}

if (!params.h5 && !params.sra_ids) {
    error "ERROR: Either --h5 or --sra_ids parameter is required"
}

log.info """
========================================================================================
    Cell-Cell Communication Spatial Pipeline
========================================================================================
    TIFF Image     : ${params.tiff}
    H5 Matrix      : ${params.h5 ?: 'N/A'}
    SRA IDs        : ${params.sra_ids ?: 'N/A'}
    scRNA Ref      : ${params.scrna_ref ?: 'None (skipping integration)'}
    Output Dir     : ${params.outdir}
========================================================================================
""".stripIndent()

// Import modules
include { PREPROCESS_IMAGE } from './modules/preprocess'
include { CELLPOSE_SEGMENT } from './modules/segment'
include { EXTRACT_ROI } from './modules/roi'
include { SPATIAL_QC } from './modules/qc'
include { CLUSTERING } from './modules/cluster'
include { CELL_COMMUNICATION } from './modules/communication'
include { ML_CLASSIFIER } from './modules/ml_classify'
include { SCRNA_INTEGRATION } from './modules/integrate'
include { COMMUNICATION_VIZ } from './modules/viz_communication'
include { ML_FEATURE_VIZ } from './modules/viz_ml'
include { FINAL_REPORT } from './modules/report'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Stage 1: Image Processing
    tiff_ch = Channel.fromPath(params.tiff, checkIfExists: true)
    preprocessed_ch = PREPROCESS_IMAGE(tiff_ch)
    mask_ch = CELLPOSE_SEGMENT(preprocessed_ch)
    roi_ch = EXTRACT_ROI(preprocessed_ch, mask_ch)
    
    // Stage 2: Load Spatial Data
    if (params.h5) {
        h5_ch = Channel.fromPath(params.h5, checkIfExists: true)
        spatial_data = h5_ch
    } else {
        // TODO: Add SRA download module if needed
        error "SRA download not yet implemented - please provide --h5"
    }
    
    // Stage 3: Spatial Analysis
    qc_data = SPATIAL_QC(spatial_data, roi_ch.coords)
    clustered_data = CLUSTERING(qc_data)
    
    // Stage 4: Cell-Cell Communication
    communication_results = CELL_COMMUNICATION(clustered_data)
    
    // Stage 5: Machine Learning Classification
    ml_results = ML_CLASSIFIER(clustered_data)
    
    // Stage 6: Optional scRNA-seq Integration
    if (params.scrna_ref) {
        ref_ch = Channel.fromPath(params.scrna_ref, checkIfExists: true)
        integrated_data = SCRNA_INTEGRATION(clustered_data, ref_ch)
    } else {
        integrated_data = clustered_data
    }
    
    // Stage 7: Visualization
    COMMUNICATION_VIZ(communication_results)
    ML_FEATURE_VIZ(ml_results)
    
    // Stage 8: Final Report
    FINAL_REPORT(
        integrated_data,
        communication_results,
        ml_results
    )
}

/*
========================================================================================
    WORKFLOW EVENTS
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================================================================
    Pipeline Execution Complete!
    ========================================================================================
    Status          : ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}
    Duration        : ${workflow.duration}
    Results         : ${params.outdir}
    
    Key Outputs:
    - Cell-cell communication networks
    - ML tissue classification (${params.outdir}/ml_classification/)
    - Integrated spatial analysis
    - Interactive visualizations
    ========================================================================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    ========================================================================================
    Pipeline Execution Failed
    ========================================================================================
    Error Message   : ${workflow.errorMessage}
    Error Report    : ${workflow.errorReport}
    ========================================================================================
    """.stripIndent()
}
