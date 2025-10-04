process ML_CLASSIFIER {
    tag "ML Tissue Classification"
    label 'process_ml'
    publishDir "${params.outdir}/ml_classification", mode: 'copy'
    
    input:
    path adata
    
    output:
    path "ml_predictions.h5ad", emit: predictions
    path "feature_importance.csv", emit: features
    path "classification_report.txt", emit: report
    path "confusion_matrix.png", emit: confusion
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score, train_test_split
    from sklearn.metrics import classification_report, confusion_matrix
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Load data
    adata = sc.read_h5ad("${adata}")
    
    # Prepare features (use highly variable genes)
    if 'highly_variable' in adata.var.columns:
        X = adata[:, adata.var['highly_variable']].X
        feature_names = adata.var_names[adata.var['highly_variable']].tolist()
    else:
        X = adata.X
        feature_names = adata.var_names.tolist()
    
    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    # Target variable (cluster labels)
    y = adata.obs['leiden'].values
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # Train Random Forest
    print("Training Random Forest classifier...")
    clf = RandomForestClassifier(
        n_estimators=100,
        max_depth=10,
        random_state=42,
        n_jobs=-1
    )
    clf.fit(X_train, y_train)
    
    # Predictions
    y_pred = clf.predict(X_test)
    adata.obs['ml_prediction'] = clf.predict(X)
    adata.obs['ml_confidence'] = clf.predict_proba(X).max(axis=1)
    
    # Save predictions
    adata.write("ml_predictions.h5ad")
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': clf.feature_importances_
    }).sort_values('importance', ascending=False)
    
    feature_importance.to_csv("feature_importance.csv", index=False)
    
    # Classification report
    with open("classification_report.txt", "w") as f:
        f.write("Machine Learning Tissue Classification Report\\n")
        f.write("=" * 60 + "\\n\\n")
        f.write(f"Model: Random Forest Classifier\\n")
        f.write(f"Features: {X.shape[1]}\\n")
        f.write(f"Training samples: {len(X_train)}\\n")
        f.write(f"Test samples: {len(X_test)}\\n\\n")
        
        # Cross-validation scores
        cv_scores = cross_val_score(clf, X_train, y_train, cv=${params.cv_folds})
        f.write(f"Cross-validation scores (${params.cv_folds}-fold):\\n")
        f.write(f"  Mean accuracy: {cv_scores.mean():.3f} (+/- {cv_scores.std():.3f})\\n\\n")
        
        # Classification report
        f.write("Classification Report:\\n")
        f.write(classification_report(y_test, y_pred))
        
        f.write("\\nTop 20 Most Important Features:\\n")
        for idx, row in feature_importance.head(20).iterrows():
            f.write(f"  {row['feature']}: {row['importance']:.4f}\\n")
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=np.unique(y), yticklabels=np.unique(y))
    plt.title('Confusion Matrix - ML Classification')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig("confusion_matrix.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✅ ML classification complete - Test accuracy: {(y_pred == y_test).mean():.3f}")
    """
}
