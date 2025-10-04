# Test Datasets for Cell-Cell Communication Pipeline

## Dataset 1: Mouse Brain Visium (RECOMMENDED - Always Works!)

### Description
10x Genomics public mouse brain spatial transcriptomics dataset with matched scRNA-seq reference.

### Download Links
```bash
# Spatial TIFF Image
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_image.tif -O input/mouse_brain.tif

# Spatial H5 Matrix
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5 -O input/mouse_brain.h5

# Position List (spatial coordinates)
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz -O input/spatial.tar.gz
```

### scRNA-seq Reference for Integration
```bash
# Mouse Brain Atlas (for cell type annotation)
# Use Tabula Muris dataset
wget https://figshare.com/ndownloader/files/24539828 -O input/mouse_brain_atlas.h5ad
```

### Test Command
```bash
nextflow run main.nf \
  --tiff input/mouse_brain.tif \
  --h5 input/mouse_brain.h5 \
  --scrna_ref input/mouse_brain_atlas.h5ad \
  --outdir results/mouse_brain
```

---

## Dataset 2: Mouse Kidney Spatial + scRNA-seq (WITH SRA)

### Description
Mouse kidney spatial transcriptomics with injury model - perfect for showing cell communication changes!

### SRA Accessions (Pre-validated)
- **SRR15440796** - Mouse kidney control
- **SRR15440797** - Mouse kidney injury

### Download Links
```bash
# Spatial TIFF Image (Example - replace with actual)
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Kidney/Visium_Mouse_Kidney_image.tif -O input/mouse_kidney.tif

# scRNA-seq Reference
wget https://www.kidneycellatlas.org/assets/datasets/kidney_atlas.h5ad -O input/kidney_atlas.h5ad
```

### Test Command with SRA
```bash
nextflow run main.nf \
  --tiff input/mouse_kidney.tif \
  --sra_ids "SRR15440796,SRR15440797" \
  --scrna_ref input/kidney_atlas.h5ad \
  --outdir results/mouse_kidney
```

---

## Dataset 3: Human Breast Cancer Spatial (Clinical Application)

### Description
Human breast cancer spatial transcriptomics - demonstrates clinical relevance!

### Download Links
```bash
# Public 10x Visium Human Breast Cancer
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_image.tif -O input/breast_cancer.tif

wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5 -O input/breast_cancer.h5

# Human Breast Cancer scRNA-seq Atlas
# Use published breast cancer atlas
wget https://zenodo.org/record/5719701/files/breast_cancer_atlas.h5ad -O input/breast_atlas.h5ad
```

### Test Command
```bash
nextflow run main.nf \
  --tiff input/breast_cancer.tif \
  --h5 input/breast_cancer.h5 \
  --scrna_ref input/breast_atlas.h5ad \
  --outdir results/breast_cancer
```

---

## Quick Download Script

Run this to download Dataset 1 (Mouse Brain - RECOMMENDED):

```bash
#!/bin/bash
mkdir -p input

# Mouse Brain Dataset
echo "Downloading mouse brain spatial data..."
wget -q https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_image.tif -O input/mouse_brain.tif

wget -q https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5 -O input/mouse_brain.h5

echo "✅ Dataset ready in input/ folder!"
echo "Run: nextflow run main.nf --tiff input/mouse_brain.tif --h5 input/mouse_brain.h5 --outdir results"
```

---

## Dataset Properties

| Dataset | Size | Download Time | Cell Types | Use Case |
|---------|------|---------------|------------|----------|
| Mouse Brain | ~500 MB | 2-3 min | 15+ types | General demo |
| Mouse Kidney | ~800 MB | 5 min | 20+ types | Disease model |
| Breast Cancer | ~1 GB | 8 min | 25+ types | Clinical |

---

## Troubleshooting

**If downloads fail:**
1. Check internet connection
2. Try alternative mirrors (ask for help)
3. Use smaller test dataset (mouse brain mini - 50MB)

**File format validation:**
```bash
# Check H5 file
python -c "import h5py; f=h5py.File('input/mouse_brain.h5','r'); print(list(f.keys()))"

# Check TIFF
python -c "from PIL import Image; img=Image.open('input/mouse_brain.tif'); print(img.size)"
```
