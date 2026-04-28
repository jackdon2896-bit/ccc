# 🧬 Golden Standard Spatial + scRNA-seq + ML + CCC Pipeline

> **Version 2.0.0** — Production-ready for Seqera Cloud

A fully-integrated Nextflow DSL2 pipeline combining:

| Layer | Tool / Method |
|---|---|
| Spatial transcriptomics | Scanpy · Squidpy · Cellpose |
| scRNA-seq (SRA) | SRA-tools · Scanpy · Harmony |
| Machine learning | Random Forest · GCN (torch-geometric) |
| Trajectory inference | PAGA · DPT (Scanpy) |
| Cell-cell communication | CellChat-style · LIANA-style · COMMOT-style |
| Multi-modal integration | Correlation transfer · Optimal Transport |

---

## 🗂️ Repository Structure

```
ccc/
├── main.nf                        # Entry point
├── nextflow.config                # Full configuration
├── params.json                    # Default parameters (your S3 paths)
├── workflows/
│   └── complete_pipeline.nf      # Master workflow
├── modules/local/
│   ├── sra_download.nf
│   ├── preprocess_spatial.nf
│   ├── preprocess_scrna.nf
│   ├── image_processing.nf        # Cellpose
│   ├── spatial_analysis.nf
│   ├── scrna_analysis.nf
│   ├── ml_analysis.nf             # Random Forest
│   ├── gnn_analysis.nf            # Graph Neural Network
│   ├── trajectory_inference.nf    # PAGA + DPT
│   ├── cell_communication.nf      # CellChat · LIANA · COMMOT
│   ├── final_integration.nf
│   └── generate_report.nf
├── bin/
│   ├── preprocess_spatial.py
│   ├── preprocess_scrna.py
│   ├── cellpose_segment.py
│   ├── spatial_analysis.py
│   ├── scrna_analysis.py
│   ├── ml_analysis.py
│   ├── gnn_analysis.py
│   ├── trajectory_inference.py
│   ├── ccc_cellchat.py
│   ├── ccc_liana.py
│   ├── ccc_commot.py
│   ├── final_integration.py
│   └── generate_report.py
├── conf/
│   ├── test.config
│   └── aws.config
├── environment.yml
└── requirements.txt
```

---

## 🚀 Quick Start — Seqera Cloud

```bash
nextflow run jackdon2896-bit/ccc \
  -params-file params.json \
  -profile seqera \
  -work-dir s3://dinesh-rdr-omics-2026/work
```

## 🔬 Inputs

| Parameter | Description | Default |
|---|---|---|
| `spatial_h5` | Spatial H5/H5AD file | `s3://.../mouse_brain.h5` |
| `spatial_image` | TIFF image | `s3://.../mouse_brain.tif` |
| `sra_accessions` | Comma-separated SRA IDs | `SRR6470906,...` |
| `outdir` | Output S3 path | `s3://.../golden_standard` |

### Brain Region Mapping (SRA samples)

| SRA IDs | Brain Region |
|---|---|
| SRR6470906–908 | Primary motor cortex |
| SRR6470910–912 | Olfactory bulb |
| SRR6470915–917 | Hippocampus |
| SRR6470923–925 | Striatum |

## 📊 Outputs

```
results/
├── preprocessed/          # QC'd spatial + per-sample scRNA H5AD
├── image_processing/      # Cellpose mask + overlay
├── spatial_analysis/      # Clustered spatial H5AD + Squidpy graph
├── scrna_analysis/        # Harmony-integrated scRNA H5AD + markers
├── machine_learning/
│   ├── rf/                # Random Forest predictions
│   └── gnn/               # GNN embeddings
├── machine_learning/trajectory/  # PAGA + pseudotime
├── cell_communication/
│   ├── cellchat/          # LR networks + communication matrix
│   ├── liana/             # Consensus interactions
│   └── commot/            # Spatial communication flows
├── final_integration/     # integrated_data.h5ad (all layers)
├── reports/               # Interactive HTML report
└── pipeline_info/         # Timeline · Trace · DAG
```

## 🔧 Container

Pipeline uses Wave container provisioning:
```
community.wave.seqera.io/library/scanpy_squidpy_cellpose_torch:latest
```

## 📄 Citation

If you use this pipeline, please cite the underlying tools:
Wolf et al. (2018) *Genome Biology* — Scanpy  
Palla et al. (2022) *Nature Methods* — Squidpy  
Stringer et al. (2021) *Nature Methods* — Cellpose  
Jin et al. (2021) *Nature Communications* — CellChat  
