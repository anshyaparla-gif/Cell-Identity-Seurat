############################################################
# 03_marker_visualization.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   1. Load the preprocessed Seurat object
#   2. Visualize marker genes for tumor, immune, stromal and endothelial populations
#   3. Save figures to the figures folder
#
############################################################

###############################
# Install / load packages
###############################

required_packages <- c("Seurat", "SeuratObject", "ggplot2", "patchwork")
installed <- rownames(installed.packages())

for(pkg in required_packages){
  if(!(pkg %in% installed)){
    install.packages(pkg)
  }
}

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

###############################
# Project directories
###############################

project_dir <- "C:/Ovarian_CAR_T_Project"
objects_dir <- file.path(project_dir, "objects")
figure_dir <- file.path(project_dir, "figures")

for(d in c(objects_dir, figure_dir)){
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}

setwd(project_dir)

###############################
# Load object
###############################

object_file <- file.path(objects_dir, "02_xenium_preprocessed.rds")

if(!file.exists(object_file)){
  object_file <- file.path(project_dir, "xenium_processed.rds")
}

if(!file.exists(object_file)){
  stop("Preprocessed object not found. Run 02_seurat_preprocessing.R first.")
}

xenium <- readRDS(object_file)

###############################
# Marker genes
###############################

marker_genes <- c(
  "MSLN",
  "MUC16",
  "EPCAM",
  "FOLR1",
  "CD3D",
  "TRAC",
  "CD8A",
  "CD4",
  "LYZ",
  "COL1A1",
  "COL3A1",
  "PECAM1",
  "VWF",
  "VEGFA",
  "CXCL12",
  "CXCR4"
)

marker_genes <- marker_genes[marker_genes %in% rownames(xenium)]

###############################
# Visualizations
###############################

if(length(marker_genes) > 0){
  p1 <- FeaturePlot(xenium, features = marker_genes[1:min(4, length(marker_genes))], combine = FALSE)
  p1 <- wrap_plots(p1, ncol = 2)
  ggsave(file.path(figure_dir, "03_marker_featureplot.pdf"), p1, width = 12, height = 8)

  p2 <- DotPlot(xenium, features = marker_genes, group.by = "seurat_clusters") +
    RotatedAxis() +
    ggtitle("Marker Expression Across Clusters")
  ggsave(file.path(figure_dir, "03_marker_dotplot.pdf"), p2, width = 14, height = 7)
}

cat("Marker visualization figures saved.\n")
