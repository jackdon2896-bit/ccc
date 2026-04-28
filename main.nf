#!/usr/bin/env nextflow

/*
========================================================================================
    Golden Standard Spatial + scRNA-seq + ML + CCC Pipeline
========================================================================================
    Integrates:
      - Spatial transcriptomics (TIFF + H5)
      - scRNA-seq (SRA download: SRR6470906-925)
      - Machine learning (GNN, RandomForest, trajectory)
      - Cell-cell communication (CellChat-style, LIANA-style, COMMOT-style)
    Author: jackdon2896-bit
    Version: 2.0.0
========================================================================================
*/

nextflow.enable.dsl = 2

include { COMPLETE_PIPELINE } from './workflows/complete_pipeline'

workflow {
    COMPLETE_PIPELINE()
}
