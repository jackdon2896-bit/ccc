# 🧬 Spatial Transcriptomics + scRNA-seq Integration Pipeline

**Production-ready Nextflow pipeline for integrating spatial transcriptomics with scRNA-seq reference data**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

---

## 🎯 Overview

This pipeline integrates **spatial transcriptomics** data with **scRNA-seq reference datasets** to enable:

- ✅ **Cell type deconvolution** in spatial data
- ✅ **Unified UMAP embedding** of both modalities
- ✅ **Cluster identification** across spatial and scRNA-seq
- ✅ **Quality control** and filtering
- ✅ **Publication-quality visualizations**

**Perfect for:** Academic researchers, industry bioinformaticians, spatial biology projects

---

## 🚀 Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/jackdon2896-bit/ccc.git
cd ccc
```

### 2. Edit Parameters

Edit `params_mode_a.json` with your S3 paths:

```json
{
  "scrna_data": "s3://your-bucket/scrna_data.h5ad",
  "spatial_data": "s3://your-bucket/spatial_data.h5ad",
  "output_dir": "s3://your-bucket/results/run_001"
}
```

### 3. Run on Seqera Cloud

**Option A: Seqera Cloud (Recommended - FREE 100 hours!)**

1. Go to https://cloud.seqera.io
2. Sign up (GitHub/Google)
3. Add this pipeline to Launchpad
4. Upload `params_mode_a.json`
5. Launch!

**See detailed guide:** [SEQERA_DEPLOYMENT_GUIDE.md](SEQERA_DEPLOYMENT_GUIDE.md)

**Option B: Local Testing**

```bash
nextflow run main.nf \
  -params-file params_mode_a.json \
  -profile docker \
  -resume
```

---

## 📊 Pipeline Workflow

```
Input Data
    ├── scRNA-seq (AnnData .h5ad)
    └── Spatial (AnnData .h5ad with spatial coords)
         ↓
    INTEGRATE_SCRNA_DATA
    ├── Load and validate data
    ├── Quality control filtering
    ├── Normalization (log1p)
    ├── Feature selection (top 2000 genes)
    ├── PCA (50 components)
    ├── Integration (concatenate)
    ├── Batch correction (harmony/scanorama)
    ├── UMAP embedding
    └── Leiden clustering
         ↓
    SPATIAL_ANALYSIS
    ├── Spatial neighbor graph
    ├── Spatial statistics
    ├── Cell type deconvolution
    ├── Spatial clustering
    └── Visualization (spatial + UMAP)
         ↓
Output
    ├── integrated_image.h5ad (main output)
    ├── integration_qc.png
    ├── integration_metrics.txt
    ├── spatial_analysis_results.h5ad
    ├── spatial_plots.png
    └── deconvolution_results.csv
```

---

## 📁 Input Formats

### scRNA-seq Data

**Supported formats:**
- `.h5ad` (AnnData - recommended)
- `.h5` (10X Genomics filtered_feature_bc_matrix.h5)
- MTX format (matrix.mtx + barcodes.tsv + features.tsv)

**Required structure:**
```python
AnnData object with:
- X: Expression matrix (cells × genes)
- obs: Cell metadata (cell_type, clusters, etc.)
- var: Gene metadata (gene_names, gene_ids)
```

### Spatial Data

**Supported formats:**
- `.h5ad` (AnnData with spatial coordinates)
- `.h5` (Visium filtered_feature_bc_matrix.h5 + spatial folder)
- Visium output directory

**Required structure:**
```python
AnnData object with:
- X: Expression matrix (spots × genes)
- obs: Spot metadata
- var: Gene metadata
- obsm['spatial']: Spatial coordinates (N × 2)
- uns['spatial']: Spatial metadata (images, scalefactors)
```

---

## ⚙️ Parameters

### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scrna_data` | path | required | Path to scRNA-seq data (S3 or local) |
| `spatial_data` | path | required | Path to spatial data (S3 or local) |
| `output_dir` | path | required | Output directory (S3 or local) |
| `mode` | string | `'full'` | Analysis mode: 'full' or 'preprocessed' |

### Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_top_genes` | int | 2000 | Number of highly variable genes |
| `n_pcs` | int | 50 | Number of principal components |
| `resolution` | float | 0.5 | Leiden clustering resolution |
| `spot_diameter` | int | 100 | Spatial spot diameter (µm) |
| `n_neighbors` | int | 15 | Number of neighbors for graphs |

### Resource Limits

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_cpus` | int | 8 | Maximum CPUs per process |
| `max_memory` | string | '64.GB' | Maximum memory per process |
| `max_time` | string | '6.h' | Maximum time per process |

---

## 📦 Outputs

### Integration Results

```
output_dir/
└── integration/
    ├── integrated_image.h5ad          ← Main output (use for Mode B!)
    ├── integration_qc.png             ← QC plots (before/after filtering)
    ├── integration_metrics.txt        ← Metrics (cells, genes, clusters)
    └── versions.txt                   ← Software versions
```

### Spatial Analysis Results

```
output_dir/
└── spatial/
    ├── spatial_analysis_results.h5ad  ← Spatial results with deconvolution
    ├── spatial_plots.png              ← Spatial visualizations
    ├── spatial_metrics.txt            ← Spatial statistics
    └── deconvolution_results.csv      ← Cell type predictions per spot
```

### Execution Reports

```
output_dir/
└── reports/
    ├── timeline.html                  ← Execution timeline
    ├── report.html                    ← Resource usage report
    ├── trace.txt                      ← Detailed trace
    └── dag.html                       ← Pipeline DAG
```

---

## 🐳 Container

Pre-built container with all dependencies:

```
quay.io/biocontainers/scanpy:1.9.3--pyhdfd78af_0
```

**Includes:**
- Python 3.10
- scanpy 1.9.3
- anndata 0.8.0
- numpy, pandas, matplotlib, seaborn
- scikit-learn
- harmonypy (batch correction)
- scanorama (alternative batch correction)

**No manual installation needed!** ✨

---

## 🔧 Configuration Profiles

### Local Testing

```bash
# Docker (default)
nextflow run main.nf -profile docker -params-file params_mode_a.json

# Singularity (HPC)
nextflow run main.nf -profile singularity -params-file params_mode_a.json
```

### Seqera Cloud

```bash
# Community Showcase (FREE 100 hours!)
nextflow run main.nf -profile seqera_showcase -params-file params_mode_a.json

# Your AWS account
nextflow run main.nf -profile awsbatch -params-file params_mode_a.json
```

---

## 💰 Cost Estimate

### Community Showcase (FREE)
- **Platform:** $0
- **Compute:** $0 (100 free CPU hours)
- **Storage:** $0 (included)
- **Total:** **$0** ✨

### Your AWS Account
- **Integration:** ~$4.80 (8 CPUs × 2h)
- **Spatial:** ~$1.20 (4 CPUs × 1h)
- **Storage:** ~$1.25 (50 GB)
- **Total:** **~$7.45 per run**

**Cost optimization:**
- Use `-resume` to avoid recomputing
- Use spot instances (70% cheaper)
- Reduce resources if data < 5000 cells

---

## 🧪 Testing

### Test with Public Data

```bash
# Download test data (10X Visium + PBMC scRNA-seq)
nextflow run main.nf \
  --scrna_data https://datasets.cellxgene.cziscience.com/... \
  --spatial_data https://cf.10xgenomics.com/samples/... \
  --output_dir test_results \
  -profile docker
```

### Validate Output

```python
import scanpy as sc

# Load integrated data
adata = sc.read_h5ad('test_results/integration/integrated_image.h5ad')

# Check structure
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"Clusters: {len(adata.obs['leiden'].unique())}")

# Visualize
sc.pl.umap(adata, color=['leiden', 'batch', 'n_genes'])
```

---

## 🐛 Troubleshooting

### Common Issues

**1. Out of Memory (exit 137)**

Increase memory in `nextflow.config`:
```groovy
process {
    withName: 'INTEGRATE_SCRNA_DATA' {
        memory = '128.GB'  // Double from 64 GB
    }
}
```

**2. S3 Access Denied**

Verify S3 permissions:
```bash
aws s3 ls s3://your-bucket/data/
```

**3. File Not Found**

Use **full S3 paths** (not directories):
```json
"scrna_data": "s3://bucket/path/file.h5ad"  ✅
"scrna_data": "s3://bucket/path/"           ❌
```

**See full guide:** [SEQERA_DEPLOYMENT_GUIDE.md](SEQERA_DEPLOYMENT_GUIDE.md)

---

## 📚 Documentation

- **[SEQERA_DEPLOYMENT_GUIDE.md](SEQERA_DEPLOYMENT_GUIDE.md)** - Complete deployment guide (814 lines!)
  - Quick start (5 minutes)
  - Academic FREE Cloud Pro access
  - Detailed setup instructions
  - Monitoring and getting results
  - Troubleshooting
  - Cost breakdown

---

## 🎓 Academic Use

**Get Seqera Cloud Pro for FREE!**

If you're in academia:
1. Go to https://seqera.io/pricing/
2. Apply for academic access
3. Get unlimited runs + professional support

**Perfect for:**
- Publications and research
- Job applications (shows cloud + bioinformatics skills)
- Building your portfolio

---

## 📖 Citation

If you use this pipeline in your research, please cite:

```
Spatial Transcriptomics + scRNA-seq Integration Pipeline
GitHub: https://github.com/jackdon2896-bit/ccc
```

**Dependencies to cite:**
- **Nextflow:** Di Tommaso et al. (2017) Nature Biotechnology
- **Scanpy:** Wolf et al. (2018) Genome Biology
- **Seqera Platform:** https://seqera.io

---

## 🤝 Contributing

Contributions welcome! Please:
1. Fork repository
2. Create feature branch
3. Test changes thoroughly
4. Submit pull request

---

## 📝 License

MIT License - see LICENSE file

---

## 💪 Support

**Issues or questions?**
- GitHub Issues: https://github.com/jackdon2896-bit/ccc/issues
- Seqera Community: https://community.seqera.io
- Nextflow Slack: Join #seqera-platform

---

## 🚀 What's Next?

After running this pipeline:

1. **Analyze results** - Load `integrated_image.h5ad` in Python/R
2. **Iterate with Mode B** - Use integrated data for faster reruns
3. **Publish findings** - Generate publication-quality figures
4. **Share pipeline** - Add to CV/resume and portfolio

**This pipeline demonstrates enterprise-level bioinformatics skills perfect for job applications!** 💼

---

**Version:** 1.0.0  
**Last Updated:** 2026-04-28  
**Maintained by:** jackdon2896-bit  
**Built with:** Nextflow DSL2, Scanpy, Seqera Platform
