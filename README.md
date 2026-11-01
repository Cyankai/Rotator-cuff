# Rotator-cuff
High-resolution spatial transcriptomic analysis (10x Genomics Visium HD, 8 µm binning) of human rotator cuff tear tissue, dissecting the spatial coupling between iron overload, hypoxia relief, and myofibroblast activation.

# Spatial Transcriptomic Architecture of Human Rotator Cuff Tear (Visium HD)

This repository contains the end-to-end R analysis pipeline and visualization scripts for high-resolution spatial transcriptomic profiling of human rotator cuff tear (RCT) specimens using the **10x Genomics Visium HD** platform (`square_008um` binned outputs). 

Our analysis resolves **11 distinct cellular subpopulations** across the tendon-to-fibrotic transition zone, reconstructs the fibroblast lineage trajectory across **three evolutionary phases**.

---

## Biological Highlights & Analytical Modules

1. **High-Resolution Spatial Atlas (8 µm Resolution)**:
   * Identification and spatial mapping of 11 niche-specific cell populations:
     * *Tendon & Quiescent Baseline*: `Tendon Core Fibroblasts`, `Quiescent Fibroblasts-1`, `Quiescent Fibroblasts-2`
     * *Remodeling & Transitional States*: `ECM-Remodeling Fibroblasts`, `Fibro-Immune Transitional`, `Pro-fibrotic Fibroblasts`
     * *Pathological Fibrotic Core*: `Active Fibroblasts`, `Myofibroblasts`
     * *Vascular & Myeloid Niches*: `Endothelial Cells`, `Ferritin-High Macrophages`, `TFRC+ Macrophages`
2. **Fibroblast Lineage Trajectory & Transcriptional Cascades**:
   * **Slingshot Pseudotime & Ridgeline Dynamics**: Continuous lineage progression across `Phase I: Homeostasis`, `Phase II: Plasticity Window`, and `Phase III: Fibrogenesis`.
   * **4 Dynamic Gene Modules (LOESS Smoothing)**:
     * *Module 1: Physiological Homeostasis* (*SCX, MKX, TNMD, DCN, FMOD, COL1A1, COL1A2, COMP*)
     * *Module 2: Iron Accumulation & Redox Stress* (*FTH1, FTL, TFRC, SLC39A14, SLC40A1, HMOX1, ACSL4, PTGS2*)
     * *Module 3: Plasticity & Crosslinking (Target)* (*LOX, LOXL1, LOXL2, POSTN, MMP2, CTHRC1, FN1, CCN2/CTGF*)
     * *Module 4: Terminal Fibrogenesis & Scarring* (*ACTA2, TAGLN, COL3A1, COL5A1, TGFB1, SERPINE1, MYL9*)
3. **Spatial Microenvironment & Intercellular Communication**:
   * **Regional DEG & Bimodal GSEA**: Differential expression between healthy tendon and fibrotic regions (`GraphPad-style Volcano Plot` & bidirectional GO BP/CC/MF enrichment).
   * **CellChat Signaling Networks**: Global intercellular communication and macrophage/vascular-driven pathways (`SPP1`, `VEGF`, `ANGPTL`).
   * **Vascular Distance & "Brake-Release" Model**: KD-tree spatial Euclidean distance analysis (`RANN`) to nearest endothelial spots, evaluating the inverse coupling between `HIF-1α Regulon Activity` and `Smad3 / Myofibroblast Activity`.

---
System Requirements & Dependencies
All scripts were developed and tested under R (version >= 4.3.3) on macOS (Apple Silicon) and Ubuntu 22.04 LTS.
Required R Packages

# Core CRAN packages
install.packages(c(
  "Seurat", "tidyverse", "dplyr", "ggplot2", "patchwork", "cowplot",
  "ggridges", "ggrepel", "ggnewscale", "scales", "RColorBrewer", "RANN", "zoo"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment", "slingshot", "ComplexHeatmap", "circlize",
  "clusterProfiler", "org.Hs.eg.db"
))

# GitHub packages (CellChat)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("jinworks/CellChat")
