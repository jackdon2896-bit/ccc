process CELLPOSE_SEGMENT {
    tag "Cellpose Segmentation"
    label 'process_medium'
    publishDir "${params.outdir}/segmentation", mode: 'copy'
    
    input:
    path image
    
    output:
    path "cell_masks.npy", emit: masks
    path "segmentation_overlay.png", emit: overlay
    path "segmentation_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import numpy as np
    from PIL import Image
    from cellpose import models
    import matplotlib.pyplot as plt
    
    # Load image
    img = np.array(Image.open("${image}"))
    
    # Run Cellpose
    model = models.Cellpose(model_type='${params.cellpose_model}')
    masks, flows, styles, diams = model.eval(
        img,
        diameter=${params.cell_diameter},
        channels=[0, 0]
    )
    
    # Save masks
    np.save("cell_masks.npy", masks)
    
    # Create overlay
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    if img.ndim == 3:
        axes[0].imshow(img)
    else:
        axes[0].imshow(img, cmap='gray')
    axes[0].set_title(f"Original Image")
    axes[0].axis('off')
    
    axes[1].imshow(masks, cmap='tab20')
    axes[1].set_title(f"Segmentation ({masks.max()} cells)")
    axes[1].axis('off')
    
    plt.tight_layout()
    plt.savefig("segmentation_overlay.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Statistics
    with open("segmentation_stats.txt", "w") as f:
        f.write(f"Total cells detected: {masks.max()}\\n")
        f.write(f"Image dimensions: {img.shape}\\n")
        f.write(f"Model: ${params.cellpose_model}\\n")
        f.write(f"Diameter: ${params.cell_diameter}\\n")
    
    print(f"✅ Detected {masks.max()} cells")
    """
}
