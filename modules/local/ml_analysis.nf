process ML_ANALYSIS {
    tag "ml_random_forest"
    label 'process_high'
    publishDir "${params.outdir}/machine_learning/rf", mode: 'copy'

    input:
    path spatial_h5ad
    path spatial_graph

    output:
    path "ml_results.h5ad",      emit: results
    path "ml_predictions.csv",   emit: predictions
    path "feature_importance.csv", emit: features
    path "ml_plots.png",         emit: plots

    script:
    """
    python3 ${projectDir}/bin/ml_analysis.py \\
        --spatial_h5ad   ${spatial_h5ad} \\
        --spatial_graph  ${spatial_graph} \\
        --out_h5ad       ml_results.h5ad \\
        --out_preds      ml_predictions.csv \\
        --out_features   feature_importance.csv \\
        --out_plots      ml_plots.png
    """
}
