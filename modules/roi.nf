process EXTRACT_ROI {
    tag "ROI Extraction"
    label 'process_low'
    publishDir "${params.outdir}/roi", mode: 'copy'
    
    input:
    path image
    path masks
    
    output:
    path "roi_coordinates.csv", emit: coords
    path "roi_visualization.png", emit: viz
    
    script:
    """
    #!/usr/bin/env python3
    import numpy as np
    import pandas as pd
    from PIL import Image
    import matplotlib.pyplot as plt
    from scipy import ndimage
    
    # Load data
    img = np.array(Image.open("${image}"))
    masks = np.load("${masks}")
    
    # Extract centroids for each cell
    cell_ids = []
    x_coords = []
    y_coords = []
    
    for cell_id in range(1, masks.max() + 1):
        cell_mask = (masks == cell_id)
        if cell_mask.sum() > 0:
            centroid = ndimage.center_of_mass(cell_mask)
            cell_ids.append(cell_id)
            y_coords.append(centroid[0])
            x_coords.append(centroid[1])
    
    # Create DataFrame
    df = pd.DataFrame({
        'cell_id': cell_ids,
        'x': x_coords,
        'y': y_coords
    })
    
    df.to_csv("roi_coordinates.csv", index=False)
    
    # Visualization
    fig, ax = plt.subplots(figsize=(10, 10))
    if img.ndim == 3:
        ax.imshow(img, alpha=0.5)
    else:
        ax.imshow(img, cmap='gray', alpha=0.5)
    
    ax.scatter(df['x'], df['y'], c='red', s=5, alpha=0.6)
    ax.set_title(f"ROI Centroids (n={len(df)})")
    ax.axis('off')
    plt.tight_layout()
    plt.savefig("roi_visualization.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Extracted {len(df)} ROI coordinates")
    """
}
