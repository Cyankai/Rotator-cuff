# Rotator-cuff
Spatial transcriptomic profiling of human rotator cuff tear tissue at single-cell scale using 10x Genomics Visium HD. Contains end-to-end R/Seurat workflows for 8um binning, clustering, cell-type mapping, and spatial niche identification.

# Spatial Transcriptomic Architecture of Human Rotator Cuff Tear (Visium HD)

This repository contains the analysis scripts and computational workflows for high-resolution spatial transcriptomic profiling of human rotator cuff tear (RCT) and enthesis specimens generated via the **10x Genomics Visium HD** platform (2 µm and 8 µm binning).

---

## Overview

Rotator cuff tears present complex spatial heterogeneity across the tendon-to-bone insertion, fibrotic transition zones, and inflammatory infiltrates. This pipeline covers:
- **Data Ingestion & Binning:** Parsing Visium HD Square Bin outputs (2 µm / 8 µm resolution) into Seurat v5 spatial objects.
- **QC & Filtering:** Low-count bin filtering, mitochondrial thresholding, and artifact removal across tissue boundaries.
- **Normalization & Dimensionality Reduction:** SCTransform / LogNormalize workflows optimized for high-density spatial grids.
- **Unsupervised Spatial Clustering:** Graph-based clustering and spatial domain identification (fibrocartilage, tendon proper, vascularized repair zones).
- **Marker Analysis & Deconvolution:** Spatial expression mapping of key tendon/fibrosis markers (e.g., *COL1A1*, *COL3A1*, *CHI3L2*, *FN1*, *MMP14*) and integration with single-cell reference datasets.

---

## Directory Structure

```text
├── scripts/
│   ├── 01_load_and_qc.R     # Space Ranger HD import, bin filtering (2um/8um)
│   ├── 02_normalization.R   # SCTransform / LogNormalize and PCA
│   ├── 03_clustering.R      # UMAP, spatial neighborhood graph, clustering
│   ├── 04_marker_analysis.R # Spatial DEGs, zone-specific markers
│   └── 05_visualization.R   # SpatialDimPlot, SpatialFeaturePlot, and export\
├── environment.yml          # Optional conda/renv environment lock
└── README.md
