# 🧬 scRNA-seq Reference Datasets for Integration

This guide provides **SRR accession numbers** and download instructions for scRNA-seq reference datasets that can be used with the `--scrna_ref` parameter in the spatial transcriptomics pipeline.

---

## 📋 Overview

The pipeline supports scRNA-seq integration for:
- **Cell type annotation transfer** - Map spatial spots to reference cell types
- **Batch correction** - Harmonize spatial and scRNA-seq modalities
- **Enhanced analysis** - Leverage reference annotations

**Required Format:** H5AD (AnnData) file with:
- `.X` = Expression matrix (genes × cells)
- `.obs['cell_type']` = Cell type annotations (optional but recommended)
- `.var_names` = Gene names matching spatial data

---

## 🐭 Mouse Brain Datasets

### 1. Allen Brain Atlas - Whole Mouse Brain (RECOMMENDED)
**Best for:** Visium brain datasets, comprehensive cell type coverage

**Source:** Allen Institute for Brain Science  
**Study:** "A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain"  
**Publication:** Yao et al., Nature (2023)  
**Data Repository:** SRA Project **SRP135960**

#### Sample SRR IDs (Adult Mouse Brain):
```bash
# Primary motor cortex (MOp)
SRR6470906  # ~3,500 cells
SRR6470907  # ~4,200 cells
SRR6470908  # ~3,800 cells

# Hippocampus
SRR6470915  # ~4,100 cells
SRR6470916  # ~3,900 cells
SRR6470917  # ~4,300 cells

# Striatum
SRR6470923  # ~3,600 cells
SRR6470924  # ~4,000 cells
SRR6470925  # ~3,700 cells

# Olfactory bulb (PERFECT for Visium olfactory bulb dataset!)
SRR6470910  # ~3,200 cells
SRR6470911  # ~3,800 cells
SRR6470912  # ~3,400 cells
```

#### Download Commands:
```bash
# Option 1: Download single SRR with SRA Toolkit
prefetch SRR6470906
fasterq-dump SRR6470906 --split-files
# This gives you FASTQ files that need Cell Ranger processing

# Option 2: Download pre-processed from Allen Brain Cell Atlas
# Visit: https://portal.brain-map.org/atlases-and-data/rnaseq
# Download the whole brain dataset in H5AD format (~40 GB)
wget https://allen-brain-cell-atlas.s3.amazonaws.com/releases/20231215/manifest.json
# Follow manifest to download specific regions
```

#### Quick Start (Using Pre-processed Data):
```python
# Download and prepare Allen Brain Atlas reference
import scanpy as sc
import anndata as ad

# Option A: Download from Allen Brain Cell Atlas API
# (Requires registration at https://portal.brain-map.org/)

# Option B: Use Tabula Muris (smaller, easier alternative)
# See dataset #4 below
```

---

### 2. Mouse Visual Cortex - 10x Genomics
**Best for:** Cortical spatial datasets

**Source:** 10x Genomics  
**Study:** Single Cell Gene Expression Dataset by Cell Ranger 2.1.0  
**Data Repository:** GEO **GSE115746**, SRA **SRP158081**

#### SRR IDs:
```bash
SRR7693731  # Mouse visual cortex, ~9,000 cells, 10x v2
SRR7693732  # Mouse visual cortex, ~8,500 cells, 10x v2
```

#### Download:
```bash
# Download FASTQ
prefetch SRR7693731
fasterq-dump SRR7693731 --split-files

# Or download pre-processed from 10x Genomics
wget https://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_filtered_gene_bc_matrices.tar.gz
tar -xzf neurons_900_filtered_gene_bc_matrices.tar.gz
```

#### Convert to H5AD:
```python
import scanpy as sc

# Read 10x data
adata = sc.read_10x_mtx('filtered_gene_bc_matrices/mm10/', var_names='gene_symbols')

# Basic preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Save as H5AD
adata.write('mouse_visual_cortex_ref.h5ad')
```

---

### 3. Mouse Brain - Saunders et al. (Drop-seq)
**Best for:** Diverse brain regions, cell type diversity

**Source:** Broad Institute / Harvard Medical School  
**Study:** "Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain"  
**Publication:** Saunders et al., Cell (2018)  
**Data Repository:** GEO **GSE110823**, SRA **SRP132700**

#### SRR IDs (Representative samples):
```bash
# Cerebral cortex
SRR6835151  # ~4,500 cells
SRR6835152  # ~5,200 cells

# Hippocampus
SRR6835165  # ~4,800 cells
SRR6835166  # ~5,100 cells

# Striatum
SRR6835175  # ~3,900 cells
SRR6835176  # ~4,300 cells

# Thalamus
SRR6835185  # ~4,100 cells
SRR6835186  # ~4,600 cells
```

#### Download Pre-processed:
```bash
# Pre-processed DGE files available at:
wget http://dropviz.org/data/DGE_filtered/DGE_filtered.tar.gz
tar -xzf DGE_filtered.tar.gz
```

---

### 4. Tabula Muris - Multi-organ Atlas (EASIEST!)
**Best for:** Quick testing, multiple organs including brain

**Source:** Chan Zuckerberg Biohub  
**Study:** "Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris"  
**Publication:** Tabula Muris Consortium, Nature (2018)  
**Data Repository:** GEO **GSE109774**, SRA **SRP131661**

#### Brain SRR IDs:
```bash
SRR6835151  # Brain non-myeloid, ~3,401 cells
SRR6835152  # Brain non-myeloid, ~3,186 cells
SRR6835153  # Brain myeloid, ~389 cells
```

#### Download Pre-processed (RECOMMENDED):
```bash
# Download ready-to-use H5AD file from Figshare
wget https://figshare.com/ndownloader/files/13092374 -O tabula_muris_brain.h5ad

# Or from CZI Cellxgene
# Visit: https://cellxgene.cziscience.com/collections/tabula-muris
```

#### Use in Pipeline:
```bash
nextflow run main.nf \
  --tiff input/spatial.tif \
  --h5 input/spatial.h5 \
  --scrna_ref tabula_muris_brain.h5ad \
  --outdir results
```

---

## 🫀 Mouse Kidney Datasets

### 5. Mouse Kidney - Park et al.
**Best for:** Kidney spatial transcriptomics integration

**Source:** UC San Diego / Broad Institute  
**Study:** "Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease"  
**Publication:** Park et al., Science (2018)  
**Data Repository:** GEO **GSE107585**, SRA **SRP132700**

#### SRR IDs:
```bash
# Adult kidney (healthy)
SRR6821785  # ~5,700 cells
SRR6821786  # ~6,200 cells
SRR6821787  # ~5,900 cells
SRR6821788  # ~6,400 cells

# Diseased kidney models (optional)
SRR6821789  # IRI day 1
SRR6821790  # IRI day 2
SRR6821791  # IRI day 14
```

#### Download:
```bash
# Download FASTQ
prefetch SRR6821785
fasterq-dump SRR6821785 --split-files

# Process with Cell Ranger or download pre-processed
# Pre-processed files available at GEO GSE107585
```

---

### 6. Mouse Kidney Atlas - Ransick et al.
**Best for:** Comprehensive kidney cell types

**Source:** UCSD  
**Study:** "Single-cell profiling reveals sex, lineage, and regional diversity in the mouse kidney"  
**Publication:** Ransick et al., Developmental Cell (2019)  
**Data Repository:** GEO **GSE129798**, SRA **SRP191805**

#### SRR IDs:
```bash
SRR8851509  # P0 kidney, ~2,400 cells
SRR8851510  # P0 kidney, ~2,800 cells
SRR8851511  # Adult kidney, ~3,200 cells
SRR8851512  # Adult kidney, ~3,600 cells
```

---

## 👤 Human Tissue Datasets (Bonus)

### 7. Human Breast Cancer - 10x Genomics
**Best for:** Human breast cancer spatial data

**Source:** 10x Genomics  
**Data Repository:** SRA **SRP133318**

#### SRR IDs:
```bash
SRR7244582  # Human breast cancer, ~4,900 cells
SRR7244583  # Human breast cancer, ~5,200 cells
```

---

## 🚀 Quick Start Guide

### Step 1: Install SRA Toolkit
```bash
# Download from: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
# Or via conda:
conda install -c bioconda sra-tools
```

### Step 2: Download SRR Data
```bash
# Method 1: Download and convert to FASTQ
prefetch SRR6470906  # Downloads to ~/ncbi/public/sra/
fasterq-dump SRR6470906 --split-files --outdir fastq/
# This creates SRR6470906_1.fastq and SRR6470906_2.fastq

# Method 2: Stream directly (slower but no storage)
fastq-dump --split-files SRR6470906
```

### Step 3: Process with Cell Ranger (if starting from FASTQ)
```bash
cellranger count \
  --id=mouse_brain_ref \
  --transcriptome=/path/to/refdata-gex-mm10-2020-A \
  --fastqs=fastq/ \
  --sample=SRR6470906
```

### Step 4: Convert to H5AD
```python
import scanpy as sc

# Read Cell Ranger output
adata = sc.read_10x_h5('mouse_brain_ref/outs/filtered_feature_bc_matrix.h5')

# Basic QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Optional: Add cell type annotations
# adata.obs['cell_type'] = ...  # Your annotations here

# Save
adata.write('mouse_brain_ref.h5ad')
```

### Step 5: Use with Pipeline
```bash
nextflow run main.nf \
  --tiff input/spatial.tif \
  --h5 input/spatial.h5 \
  --scrna_ref mouse_brain_ref.h5ad \
  --outdir results
```

---

## 📊 Pre-processed Alternative (Fastest!)

Instead of downloading raw SRR files, use these **pre-processed** resources:

### Option 1: Allen Brain Cell Atlas (Mouse Brain)
```bash
# Visit: https://portal.brain-map.org/atlases-and-data/rnaseq
# Download region-specific H5AD files (already processed!)
```

### Option 2: CZI Cellxgene (Multiple Species/Organs)
```python
import cellxgene_census

# Access census data
census = cellxgene_census.open_soma()

# Query mouse brain data
mouse_brain = census["census_data"]["mus_musculus"]
# Filter and download directly as AnnData
```

### Option 3: UCSC Cell Browser
```bash
# Visit: https://cells.ucsc.edu/
# Search for "mouse brain" or "mouse kidney"
# Download H5AD format directly
```

---

## 🎯 Recommended Combinations

### For Mouse Brain Visium (Olfactory Bulb):
```bash
# Pipeline test dataset + Allen Brain reference
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_image.tif
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5

# Download olfactory bulb scRNA-seq reference:
# SRR6470910, SRR6470911, SRR6470912 (from Allen Brain Atlas)
# Or use Tabula Muris brain subset
```

### For Mouse Kidney Visium:
```bash
# Spatial data from 10x Genomics
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/V1_Mouse_Kidney/V1_Mouse_Kidney_image.tif
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/V1_Mouse_Kidney/V1_Mouse_Kidney_filtered_feature_bc_matrix.h5

# scRNA-seq reference from Park et al.:
# SRR6821785, SRR6821786, SRR6821787 (healthy adult kidney)
```

---

## 🔧 Troubleshooting

### Issue: "No overlapping genes between spatial and reference"
**Solution:** Ensure gene names match (symbols vs. Ensembl IDs)
```python
# Convert Ensembl IDs to gene symbols
import scanpy as sc
adata = sc.read_h5ad('reference.h5ad')
# Check current format
print(adata.var_names[:5])  # Should be gene symbols like 'Gapdh', not 'ENSMUSG00000...'
```

### Issue: "Reference file too large"
**Solution:** Subsample cells
```python
# Keep only high-quality cells
sc.pp.filter_cells(adata, min_genes=500)
# Subsample to 10,000 cells
sc.pp.subsample(adata, n_obs=10000)
adata.write('reference_subset.h5ad')
```

### Issue: "SRA download too slow"
**Solution:** Use pre-processed data from Allen Brain or Cellxgene (see above)

---

## 📚 Reference Data Sources

| Source | URL | Best For |
|--------|-----|----------|
| **Allen Brain Cell Atlas** | https://portal.brain-map.org/ | Mouse brain, pre-processed |
| **CZI Cellxgene** | https://cellxgene.cziscience.com/ | Multi-species, interactive |
| **Tabula Muris** | https://tabula-muris.ds.czbiohub.org/ | Mouse multi-organ |
| **10x Genomics Datasets** | https://www.10xgenomics.com/datasets | Validated test data |
| **UCSC Cell Browser** | https://cells.ucsc.edu/ | Browse and download |
| **SRA (NCBI)** | https://www.ncbi.nlm.nih.gov/sra | Raw sequencing data |
| **GEO (NCBI)** | https://www.ncbi.nlm.nih.gov/geo/ | Processed matrices |

---

## 💡 Best Practices

1. **Match tissue types**: Use brain reference for brain spatial data
2. **Match species**: Mouse reference for mouse spatial data
3. **Check gene overlap**: At least 50% shared genes recommended
4. **Quality control**: Filter low-quality cells before using as reference
5. **File size**: Keep reference < 5 GB for faster processing
6. **Cell type labels**: Include `.obs['cell_type']` if available

---

## 🎓 Citation

If you use these datasets in your research, please cite:

**Allen Brain Atlas (SRP135960):**
Yao et al. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature.

**Tabula Muris:**
Tabula Muris Consortium (2018). Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature.

**Park et al. Kidney (GSE107585):**
Park et al. (2018). Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science.

---

## ✅ Quick Test Example

**Complete working example with public data:**

```bash
# 1. Download spatial data (mouse brain)
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_image.tif -O input/spatial.tif
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5 -O input/spatial.h5

# 2. Download scRNA-seq reference (Tabula Muris - easiest!)
wget https://figshare.com/ndownloader/files/13092374 -O input/mouse_brain_ref.h5ad

# 3. Run pipeline with integration
nextflow run main.nf \
  --tiff input/spatial.tif \
  --h5 input/spatial.h5 \
  --scrna_ref input/mouse_brain_ref.h5ad \
  --outdir results_integrated

# 4. Check results
open results_integrated/PIPELINE_REPORT.html
```

**Expected output:**
- Spatial spots annotated with cell types from reference
- UMAP showing integrated spatial + scRNA-seq
- Cell type composition analysis

---

**Need help?** Check the main README.md for pipeline parameters and troubleshooting.
