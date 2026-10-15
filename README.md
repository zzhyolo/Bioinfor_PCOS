# Integrated Network Toxicology and Multi-Dimensional Transcriptomics in DEHP-Associated PCOS

This repository provides the complete analysis scripts and computational workflow for the manuscript:

> **Integrated network toxicology and multi-dimensional transcriptomics reveal DEHP-associated pathogenic mechanisms in Polycystic Ovary Syndrome**

---

## 📂 Repository Structure

The entire analytical workflow is consolidated into two streamlined, self-contained R scripts:

### 1. `01_bulk_transcriptomics_analysis.R`
This script covers the complete bulk transcriptomic, network toxicology, and machine-learning prioritization pipeline:
- **Part 1. Multi-Cohort Integration & Batch Effect Correction**: Combines bulk PCOS and control granulosa-cell transcriptomes (GSE138518, GSE193123, GSE168404) and removes platform batch effects using `ComBat` (`sva` package).
- **Part 2. Differential Expression Analysis**: Identifies differentially expressed genes (DEGs) via `limma` empirical Bayes linear modeling (`|log2FC| > 0.2`, `P < 0.05`) and generates the volcano plot (Figure 2A).
- **Part 3. Candidate Target Intersection & Enrichment**: Intersects DEGs with predicted DEHP targets from public toxicological databases (CTD, STITCH, SwissTargetPrediction, ChEMBL) and performs GO and KEGG pathway enrichment analysis via `clusterProfiler`.
- **Part 4. Core Target Prioritization via Machine Learning**: Implements LASSO logistic regression (10-fold cross-validation, 1-SE rule) and Random Forest feature importance ranking (ntree = 100, mtry = 4, MeanDecreaseGini) to prioritize 22 core candidate genes.
- **Part 5. Consensus Clustering & GSEA**: Conducts consensus clustering ($k = 2$ to $6$, 500 iterations) using `ConsensusClusterPlus`, generates hierarchical clustering heatmaps (`pheatmap`), and performs GSEA across MSigDB Hallmark pathways.

---

### 2. `02_single_cell_spatial_analysis.R`
This script implements single-cell resolution, perturbation simulation, and spatial co-localization analyses:
- **Part 1. Single-Cell RNA Sequencing & Annotation**: Processes normal human ovarian single-cell data (GSE255690), removes inter-sample batch effects using `Harmony` (dims = 1:20), and annotates 8 major ovarian lineages based on canonical lineage markers.
- **Part 2. Single-Cell Virtual Knockout**: Simulates in silico gene deletion of *APOE* and *DHCR7* within the single-cell gene regulatory network using `scTenifoldKnk`.
- **Part 3. Intercellular Communication Modeling**: Reconstructs cell–cell signaling networks with `CellChat` (truncatedMean method, trim = 0.1, min.cells = 1), evaluating the APOE–TREM2 ligand–receptor pair and the DHCR7-dependent desmosterol–NR1H2 and cholesterol–RORA metabolite–nuclear receptor axes.
- **Part 4. Ovarian Spatial Transcriptomics & Co-Localization**: Normalizes human ovarian Stereo-seq data (GSE267315) via `SCTransform`, performs spatial clustering, and maps high-resolution spatial co-localization of target–receptor pairs using `semla`.

---

## 💻 System & Software Requirements

- **R version**: $\ge 4.4.0$
- **Key Bioinformatic Packages**:
  - `limma`, `sva`, `Biobase`
  - `clusterProfiler`, `org.Hs.eg.db`, `msigdbr`
  - `glmnet`, `randomForest`, `caret`, `pROC`
  - `ConsensusClusterPlus`, `pheatmap`, `ggplot2`, `tidyverse`
  - `Seurat`, `harmony`, `CellChat`, `scTenifoldNet` (for scTenifoldKnk), `semla`

---

## 📄 License & Citation

The scripts in this repository are provided for academic research and computational reproducibility. Please cite the original manuscript when referencing these workflows.
