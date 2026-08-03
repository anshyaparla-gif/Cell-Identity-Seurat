# Ovarian Cancer Spatial Transcriptomics Analysis

This repository contains a bioinformatics workflow for analyzing ovarian cancer spatial transcriptomics data using the 10X Genomics Xenium platform and Seurat in R.

## Project Goal

The objective of this project is to investigate the ovarian cancer tumor microenvironment by:

- Identifying major cell populations
- Examining spatial relationships between tumor and immune cells
- Evaluating potential CAR-T therapeutic targets
- Performing downstream spatial and cell-cell communication analyses


# Step 1: Environment Software Setup 

## Software Installed

- R (latest stable version)
- RStudio Desktop
---
## Download Public Data
- https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer

## Required R Packages

The following packages were installed:

```r
install.packages(c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork"
))
```

Additional packages used later include:

```r
install.packages(c(
  "data.table",
  "Matrix",
  "future"
))
```

---

# Dataset

The project uses the public 10X Genomics Xenium dataset:

**FFPE Human Ovarian Cancer with 5K Human Pan Tissue**

Dataset contents include:

- cell_feature_matrix.h5
- cells.parquet
- transcripts.parquet
- experiment.xenium

---

# Loading Data

The expression matrix was loaded into Seurat using:

```r
library(Seurat)

counts <- Read10X_h5("cell_feature_matrix.h5")

xenium <- CreateSeuratObject(
    counts = counts,
    project = "Ovarian_Xenium"
)
```

---

# Initial Quality Check

Basic statistics were examined:

```r
xenium
```

Typical information includes:

- Number of genes
- Number of cells
- Count matrix dimensions

---

# First Analysis

Normalization:

```r
xenium <- NormalizeData(xenium)
```

Identification of variable features:

```r
xenium <- FindVariableFeatures(xenium)
```

Scaling:

```r
xenium <- ScaleData(xenium)
```

Principal Component Analysis:

```r
xenium <- RunPCA(xenium)
```

Neighbor graph:

```r
xenium <- FindNeighbors(xenium, dims = 1:20)
```

Cell clustering:

```r
xenium <- FindClusters(xenium)
```

UMAP visualization:

```r
xenium <- RunUMAP(xenium, dims = 1:20)

DimPlot(xenium)
```

This produced the first visualization of cellular clusters within the Xenium dataset.

---

# Repository Structure

```
.
├── data/
│   └── Xenium dataset (not committed)
├── scripts/
│   └── *.R
├── figures/
├── results/
└── README.md
```

---

# Notes

The raw Xenium dataset is several gigabytes in size and is **not included** in this repository. Download the dataset directly from the 10X Genomics website before running the analysis and link is included above.

---

# Future Work

Planned analyses include:

- Cell type annotation
- Tumor marker visualization
- Spatial distance analysis
- CAR-T target identification
- CellChat ligand-receptor communication
- Statistical analysis of immune infiltration
- Publication-quality figures