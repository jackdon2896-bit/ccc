# Container Image Documentation

## Container Details

**Image URI:**
```
community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126
```

**Build ID:** `bd-f40b275367e7e126_1`

**Build Status:** Successfully built via Wave by Seqera

**Container Format:** Docker (Docker Hub compatible)

---

## Included Packages

### Core Python Environment
- **Python:** 3.11 (compact installation)

### Spatial Analysis Tools
- **scanpy:** Single-cell analysis framework
- **squidpy:** Spatial molecular data analysis
- **anndata:** Annotated data structures for computational biology
- **cellpose:** Deep learning-based cell segmentation

### Cell Communication Analysis
- **liana-py:** Ligand-receptor interaction analysis
- **cellphonedb:** Cell-cell communication database and methods
- **CellChat integration:** Via liana-py interface

### Machine Learning & Deep Learning
- **scvi-tools:** Single-cell variational inference
- **celltypist:** Automated cell type annotation
- **harmonypy:** Batch effect correction
- **scikit-learn:** Machine learning algorithms
- **pytorch:** Deep learning framework

### Sequencing Tools
- **sra-tools:** SRA data download utilities
- **fastqc:** Quality control for sequencing data
- **multiqc:** Aggregate QC reports

### Data Processing
- **pandas:** Data manipulation
- **numpy:** Numerical computing
- **scipy:** Scientific computing
- **imageio:** Image reading/writing
- **pillow:** Image processing
- **matplotlib:** Plotting
- **seaborn:** Statistical visualization

---

## Usage

### Pull the Container

**With Docker:**
```bash
docker pull community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126
```

**With Singularity:**
```bash
singularity pull docker://community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126
```

### Run Interactive Shell

**Docker:**
```bash
docker run -it --rm \
  -v $(pwd):/workspace \
  -w /workspace \
  community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  /bin/bash
```

**Singularity:**
```bash
singularity shell \
  --bind $(pwd):/workspace \
  --pwd /workspace \
  docker://community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126
```

### Execute Python Scripts

**Docker:**
```bash
docker run --rm \
  -v $(pwd):/workspace \
  -w /workspace \
  community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  python your_script.py
```

---

## Nextflow Configuration

This container is automatically configured in `nextflow.config`:

```groovy
process {
    container = 'community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126'
}
```

The pipeline will automatically use this container for all processes.

---

## Container Size & Performance

- **Optimized:** Python 3.11 compact version for smaller image size
- **Multi-architecture:** Compatible with x86_64 systems
- **Cached:** Pre-built and cached by Wave for instant availability
- **Network:** Hosted on Seqera's Wave CDN for fast downloads

---

## Verification

### Check Python Version
```bash
docker run --rm community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 python --version
```

### Verify Package Installation
```bash
docker run --rm community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  python -c "import scanpy, liana, cellpose; print('✅ All packages loaded')"
```

### List All Installed Packages
```bash
docker run --rm community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  pip list
```

---

## Troubleshooting

### Permission Issues (Docker)
If you encounter permission issues with Docker, add user mapping:
```bash
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/workspace \
  -w /workspace \
  community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  python your_script.py
```

### GPU Access (for Cellpose)
To enable GPU support for Cellpose segmentation:
```bash
docker run --rm --gpus all \
  -v $(pwd):/workspace \
  -w /workspace \
  community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  python your_script.py
```

### Memory Limits
For large datasets, increase Docker memory:
```bash
docker run --rm -m 16g \
  -v $(pwd):/workspace \
  -w /workspace \
  community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126 \
  python your_script.py
```

---

## Building Custom Containers

If you need additional packages, you can extend this container:

```dockerfile
FROM community.wave.seqera.io/library/anndata_cellphonedb_cellpose_celltypist_pruned:f40b275367e7e126

# Add your custom packages
RUN pip install your-package-name

# Or conda packages
RUN conda install -c conda-forge your-conda-package
```

---

## Support

- **Wave Documentation:** https://www.nextflow.io/docs/latest/wave.html
- **Seqera Platform:** https://seqera.io
- **Container Issues:** Check build logs at Wave dashboard or contact Seqera support

---

**Last Updated:** 2025
**Maintained By:** Seqera AI
**License:** All software packages follow their respective licenses
