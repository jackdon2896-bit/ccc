# Cell-Cell Communication Spatial Transcriptomics Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)](https://www.docker.com/)
[![Wave](https://img.shields.io/badge/wave-enabled-orange.svg)](https://seqera.io/wave/)

Advanced bioinformatics pipeline for spatial transcriptomics analysis with integrated cell-cell communication inference and machine learning tissue classification.

## 🚀 Features

- **📊 Spatial Transcriptomics Analysis**
  - Automated image preprocessing and quality control
  - Cellpose-based cell segmentation
  - ROI extraction and coordinate mapping

- **🔬 Cell-Cell Communication**
  - LIANA integration for ligand-receptor analysis
  - CellPhoneDB database support
  - Network visualization and chord diagrams

- **🤖 Machine Learning Classification**
  - Random Forest tissue classification
  - Feature importance analysis
  - Cross-validation with confidence metrics

- **🧬 Multi-modal Integration**
  - scRNA-seq reference integration
  - Harmony batch correction
  - Cell type annotation transfer

- **☁️ Cloud-Ready**
  - AWS Batch deployment configured
  - Seqera Platform compatible
  - Wave container provisioning
  - Fusion file system support

---

## 📋 Requirements

- **Nextflow:** >= 25.04.0
- **Container Engine:** Docker or Singularity
- **Resources:** 16+ CPUs, 64GB RAM recommended

---

## 🛠️ Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/jackdon2896-bit/ccc.git
cd ccc
```

### 2. Download Test Dataset
```bash
# Mouse brain dataset (RECOMMENDED)
wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_image.tif -O input/mouse_brain.tif

wget https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5 -O input/mouse_brain.h5
```

### 3. Run Pipeline
```bash
nextflow run main.nf \
  --tiff input/mouse_brain.tif \
  --h5 input/mouse_brain.h5 \
  --outdir results
```

### 4. View Results
Open `results/PIPELINE_REPORT.html` in your browser!

---

## 💻 Usage

### Basic Command
```bash
nextflow run main.nf \
  --tiff <path_to_tiff_image> \
  --h5 <path_to_h5_matrix> \
  --outdir <output_directory>
```

### With scRNA-seq Integration
```bash
nextflow run main.nf \
  --tiff input/spatial.tif \
  --h5 input/spatial.h5 \
  --scrna_ref input/reference.h5ad \
  --outdir results
```

### With Custom Parameters
```bash
nextflow run main.nf \
  --tiff input/spatial.tif \
  --h5 input/spatial.h5 \
  --min_genes 500 \
  --max_mito_pct 15 \
  --resolution 0.6 \
  --outdir results
```

---

## 🎛️ Parameters

### Required Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `--tiff` | path | Spatial TIFF image file |
| `--h5` or `--sra_ids` | path/string | H5 matrix file or SRA accessions |

### Optional Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--scrna_ref` | path | null | scRNA-seq reference for integration |
| `--outdir` | path | results | Output directory |
| `--min_genes` | int | 200 | Minimum genes per cell |
| `--min_cells` | int | 3 | Minimum cells per gene |
| `--max_mito_pct` | int | 20 | Maximum mitochondrial % |
| `--cellpose_model` | string | cyto2 | Cellpose model type |
| `--cell_diameter` | int | 30 | Expected cell diameter |
| `--n_neighbors` | int | 15 | Neighbors for graph construction |
| `--resolution` | float | 0.8 | Leiden clustering resolution |
| `--comm_method` | string | liana | Communication method |
| `--ml_model` | string | random_forest | ML classifier type |
| `--cv_folds` | int | 5 | Cross-validation folds |

---

## 📁 Output Structure

```
results/
├── PIPELINE_REPORT.html          # Main interactive report
├── summary_statistics.txt         # Text summary
├── preprocessing/
│   ├── preprocessed.tif
│   └── preprocessing_qc.png
├── segmentation/
│   ├── cell_masks.npy
│   ├── segmentation_overlay.png
│   └── segmentation_stats.txt
├── roi/
│   ├── roi_coordinates.csv
│   └── roi_visualization.png
├── qc/
│   ├── spatial_adata.h5ad
│   ├── qc_metrics.png
│   └── qc_summary.txt
├── clustering/
│   ├── clustered_adata.h5ad
│   ├── umap_clusters.png
│   └── cluster_stats.txt
├── communication/
│   ├── communication_results.h5ad
│   ├── ligand_receptor_network.csv
│   └── communication_summary.txt
├── ml_classification/
│   ├── ml_predictions.h5ad
│   ├── feature_importance.csv
│   ├── classification_report.txt
│   └── confusion_matrix.png
├── integration/              # If --scrna_ref provided
│   ├── integrated_adata.h5ad
│   ├── integration_plot.png
│   └── integration_summary.txt
└── visualizations/
    ├── network_plot.png
    ├── dotplot_interactions.png
    ├── chord_diagram.png
    ├── feature_importance_plot.png
    ├── prediction_confidence.png
    └── cluster_comparison.png
```

---

## 🐳 Container Information

**Pre-built Wave Container:**
```
community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126
```

**Included Tools:**
- Python 3.11 (compact)
- Scanpy, Squidpy, Cellpose
- LIANA, CellPhoneDB
- scvi-tools, Celltypist
- scikit-learn, PyTorch
- sra-tools, FastQC, MultiQC

See [CONTAINER.md](CONTAINER.md) for full details.

---

## ☁️ Cloud Deployment

### AWS Batch
```bash
nextflow run main.nf \
  -profile aws \
  --tiff s3://my-bucket/input/spatial.tif \
  --h5 s3://my-bucket/input/spatial.h5 \
  --outdir s3://my-bucket/results
```

### Seqera Platform
```bash
nextflow run main.nf \
  -profile seqera \
  --tiff <data_link_path> \
  --h5 <data_link_path> \
  --outdir <output_path>
```

Configure AWS details in `nextflow.config` before running.

---

## 📊 Test Datasets

See [input/DATASETS.md](input/DATASETS.md) for curated test datasets:

1. **Mouse Brain Visium** (RECOMMENDED)
   - Always works, well-characterized
   - ~500 MB, 15+ cell types
   
2. **Mouse Kidney Spatial**
   - Disease model with SRA data
   - ~800 MB, 20+ cell types
   
3. **Human Breast Cancer**
   - Clinical application example
   - ~1 GB, 25+ cell types

---

## 🔧 Configuration Profiles

### Local Execution (Docker)
```bash
nextflow run main.nf -profile docker
```

### Local Execution (Singularity)
```bash
nextflow run main.nf -profile singularity
```

### AWS Batch
```bash
nextflow run main.nf -profile aws
```

### Seqera Platform
```bash
nextflow run main.nf -profile seqera
```

### Test Profile
```bash
nextflow run main.nf -profile test
```

---

## 📈 Performance

**Typical Runtime (Mouse Brain Dataset):**
- Local (16 CPUs, 64GB RAM): ~30-45 minutes
- AWS Batch (m5.4xlarge): ~20-30 minutes
- Seqera Platform (optimized): ~15-25 minutes

**Resource Requirements by Module:**
- Preprocessing: 2 CPUs, 4GB RAM
- Segmentation: 4 CPUs, 16GB RAM
- Communication: 8 CPUs, 32GB RAM
- ML Classification: 8 CPUs, 32GB RAM

---

## 🐛 Troubleshooting

### Issue: Container pull fails
**Solution:**
```bash
# Pre-pull the container
docker pull community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126
```

### Issue: Out of memory
**Solution:** Increase Docker memory or use cluster profile:
```bash
nextflow run main.nf -profile aws  # Use cloud resources
```

### Issue: No significant interactions found
**Solution:** Adjust communication parameters:
```bash
nextflow run main.nf --comm_method liana --min_spots 3
```

### Issue: Permission denied (Docker)
**Solution:** Add user mapping in `nextflow.config`:
```groovy
docker.runOptions = '-u $(id -u):$(id -g)'
```

---

## 📚 Citation

If you use this pipeline, please cite:

- **Nextflow:** Di Tommaso et al., Nature Biotechnology 2017
- **Scanpy:** Wolf et al., Genome Biology 2018
- **LIANA:** Dimitrov et al., Nature Methods 2022
- **Cellpose:** Stringer et al., Nature Methods 2021
- **Wave:** Seqera Labs (https://seqera.io/wave/)

---

## 📝 License

MIT License - see LICENSE file for details

---

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---

## 💬 Support

- **Issues:** https://github.com/jackdon2896-bit/ccc/issues
- **Discussions:** https://github.com/jackdon2896-bit/ccc/discussions
- **Email:** seqera-ai@seqera.io

---

## 🙏 Acknowledgments

- Seqera Platform team for Wave container service
- 10x Genomics for public test datasets
- Scanpy, LIANA, and Cellpose development teams

---

**Built with ❤️ by Seqera AI**
