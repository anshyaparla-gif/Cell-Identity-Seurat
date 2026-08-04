############################################################
# 02_seurat_preprocessing.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   1. Load the raw Xenium Seurat object
#   2. Normalize and scale gene expression
#   3. Run dimensional reduction and clustering
#   4. Save the preprocessed object
#
############################################################

###############################
# Install / load packages
###############################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "dplyr",
  "patchwork",
  "future"
)

installed <- rownames(installed.packages())

for(pkg in required_packages){
  if(!(pkg %in% installed)){
    install.packages(pkg)
  }
}

library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(future)

###############################
# Project directories
###############################

project_dir <- "C:/Ovarian_CAR_T_Project"
objects_dir <- file.path(project_dir, "objects")
results_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")

for(d in c(objects_dir, results_dir, figure_dir)){
  if(!dir.exists(d)) dir.create(d, recursive = TRUE)
}

setwd(project_dir)

###############################
# Load raw Seurat object
###############################

raw_object <- file.path(objects_dir, "01_xenium_raw.rds")

if(file.exists(raw_object)){
  xenium <- readRDS(raw_object)
} else if(file.exists(file.path(project_dir, "xenium_processed.rds"))){
  xenium <- readRDS(file.path(project_dir, "xenium_processed.rds"))
} else {
  stop("Raw Seurat object not found. Run 01_setup.R first.")
}

###############################
# Preprocessing workflow
###############################

cat("Normalizing data...\n")
xenium <- NormalizeData(xenium, normalization.method = "LogNormalize", scale.factor = 10000)

cat("Finding variable features...\n")
xenium <- FindVariableFeatures(xenium, selection.method = "vst", nfeatures = 2000)

cat("Scaling data...\n")
xenium <- ScaleData(xenium, features = VariableFeatures(xenium))

cat("Running PCA...\n")
xenium <- RunPCA(xenium, npcs = 50)

cat("Finding neighbors...\n")
xenium <- FindNeighbors(xenium, dims = 1:20)

cat("Clustering cells...\n")
xenium <- FindClusters(xenium, resolution = 0.6)

cat("Running UMAP...\n")
xenium <- RunUMAP(xenium, dims = 1:20)

###############################
# Save processed object
###############################

saveRDS(xenium, file.path(objects_dir, "02_xenium_preprocessed.rds"))
saveRDS(xenium, file.path(project_dir, "xenium_processed.rds"))

cat("\nPreprocessed Seurat object saved.\n")
