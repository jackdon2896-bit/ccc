/*
========================================================================================
    Golden Standard Complete Workflow
    Spatial + scRNA-seq + ML + Cell-Cell Communication
========================================================================================
*/

include { SRA_DOWNLOAD         } from '../modules/local/sra_download'
include { PREPROCESS_SPATIAL   } from '../modules/local/preprocess_spatial'
include { PREPROCESS_SCRNA     } from '../modules/local/preprocess_scrna'
include { CELLPOSE_SEGMENT     } from '../modules/local/image_processing'
include { SCRNA_ANALYSIS       } from '../modules/local/scrna_analysis'
include { SPATIAL_ANALYSIS     } from '../modules/local/spatial_analysis'
include { ML_ANALYSIS          } from '../modules/local/ml_analysis'
include { GNN_ANALYSIS         } from '../modules/local/gnn_analysis'
include { TRAJECTORY_INFERENCE } from '../modules/local/trajectory_inference'
include { CCC_CELLCHAT         } from '../modules/local/cell_communication'
include { CCC_LIANA            } from '../modules/local/cell_communication'
include { CCC_COMMOT           } from '../modules/local/cell_communication'
include { FINAL_INTEGRATION    } from '../modules/local/final_integration'
include { GENERATE_REPORT      } from '../modules/local/generate_report'

// Brain region map
def BRAIN_REGIONS = [
    'SRR6470906': 'primary_motor_cortex',
    'SRR6470907': 'primary_motor_cortex',
    'SRR6470908': 'primary_motor_cortex',
    'SRR6470910': 'olfactory_bulb',
    'SRR6470911': 'olfactory_bulb',
    'SRR6470912': 'olfactory_bulb',
    'SRR6470915': 'hippocampus',
    'SRR6470916': 'hippocampus',
    'SRR6470917': 'hippocampus',
    'SRR6470923': 'striatum',
    'SRR6470924': 'striatum',
    'SRR6470925': 'striatum'
]

workflow COMPLETE_PIPELINE {

    // ── Input channels ──────────────────────────────────────────────────────
    def ch_spatial_h5    = channel.fromPath(params.spatial_h5,    checkIfExists: true)
    def ch_spatial_image = channel.fromPath(params.spatial_image, checkIfExists: true)

    def ch_sra = channel.of(params.sra_accessions.tokenize(','))
        .flatten()
        .map { sra_id ->
            def region = BRAIN_REGIONS.getOrDefault(sra_id, 'unknown')
            tuple(sra_id, region)
        }

    // ── Step 1: SRA download ────────────────────────────────────────────────
    def ch_fastq = SRA_DOWNLOAD(ch_sra)

    // ── Step 2: Preprocessing ───────────────────────────────────────────────
    def ch_spatial_processed = PREPROCESS_SPATIAL(ch_spatial_h5, ch_spatial_image)
    def ch_scrna_processed   = PREPROCESS_SCRNA(ch_fastq.fastq)

    // ── Step 3: Image segmentation ──────────────────────────────────────────
    def ch_segmented = CELLPOSE_SEGMENT(ch_spatial_processed.image)

    // ── Step 4: Core analysis ───────────────────────────────────────────────
    def ch_spatial_results = SPATIAL_ANALYSIS(
        ch_spatial_processed.h5ad,
        ch_segmented.mask
    )

    def ch_scrna_results = SCRNA_ANALYSIS(
        ch_scrna_processed.h5ad.map { sra_id, region, h5ad -> h5ad }.collect()
    )

    // ── Step 5: Machine learning ────────────────────────────────────────────
    def ch_ml = ML_ANALYSIS(
        ch_spatial_results.processed,
        ch_spatial_results.spatial_graph
    )

    def ch_gnn = GNN_ANALYSIS(
        ch_spatial_results.processed,
        ch_spatial_results.spatial_graph
    )

    def ch_traj = TRAJECTORY_INFERENCE(
        ch_spatial_results.processed,
        ch_scrna_results.processed
    )

    // ── Step 6: Cell-cell communication ────────────────────────────────────
    def ch_cellchat = CCC_CELLCHAT(
        ch_scrna_results.processed,
        ch_spatial_results.processed
    )

    def ch_liana = CCC_LIANA(
        ch_scrna_results.processed,
        ch_spatial_results.processed
    )

    def ch_commot = CCC_COMMOT(
        ch_spatial_results.processed,
        ch_spatial_results.spatial_graph
    )

    // ── Step 7: Final integration ───────────────────────────────────────────
    def ch_all = channel.empty()
        .mix(ch_spatial_results.processed.map  { f -> tuple('spatial',   f) })
        .mix(ch_scrna_results.processed.map     { f -> tuple('scrna',     f) })
        .mix(ch_ml.results.map                  { f -> tuple('ml',        f) })
        .mix(ch_gnn.embeddings.map              { f -> tuple('gnn',       f) })
        .mix(ch_traj.pseudotime.map             { f -> tuple('traj',      f) })
        .mix(ch_cellchat.networks.map           { f -> tuple('cellchat',  f) })
        .mix(ch_liana.interactions.map          { f -> tuple('liana',     f) })
        .mix(ch_commot.flows.map                { f -> tuple('commot',    f) })
        .collect()

    def ch_integrated = FINAL_INTEGRATION(
        ch_all,
        ch_spatial_processed.image
    )

    // ── Step 8: Report ──────────────────────────────────────────────────────
    GENERATE_REPORT(
        ch_integrated.data,
        ch_integrated.plots,
        ch_integrated.metrics
    )

    emit:
    integrated_data = ch_integrated.data
    report          = GENERATE_REPORT.out.html
}
