process PREPROCESS_IMAGE {
    tag "Preprocessing"
    label 'process_low'
    publishDir "${params.outdir}/preprocessing", mode: 'copy'
    
    input:
    path tiff
    
    output:
    path "preprocessed.tif", emit: image
    path "preprocessing_qc.png", emit: qc_plot
    
    script:
    """
    #!/usr/bin/env python3
    import numpy as np
    from PIL import Image
    import matplotlib.pyplot as plt
    
    # Load image
    img = Image.open("${tiff}")
    img_array = np.array(img)
    
    # Normalize
    if img_array.ndim == 3:
        img_normalized = (img_array - img_array.min()) / (img_array.max() - img_array.min())
    else:
        img_normalized = (img_array - np.min(img_array)) / (np.max(img_array) - np.min(img_array))
    
    # Convert to 8-bit
    img_8bit = (img_normalized * 255).astype(np.uint8)
    
    # Save preprocessed image
    Image.fromarray(img_8bit).save("preprocessed.tif")
    
    # QC plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    if img_array.ndim == 3:
        axes[0].imshow(img_array)
        axes[1].imshow(img_8bit)
    else:
        axes[0].imshow(img_array, cmap='gray')
        axes[1].imshow(img_8bit, cmap='gray')
    
    axes[0].set_title("Original")
    axes[1].set_title("Preprocessed")
    axes[0].axis('off')
    axes[1].axis('off')
    plt.tight_layout()
    plt.savefig("preprocessing_qc.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print("✅ Image preprocessing complete")
    """
}
