############################################################
# 01_setup.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Platform: 10X Genomics Xenium
#
# Purpose:
#   1. Install/load required packages
#   2. Set project directories
#   3. Load Xenium gene expression data
#   4. Create a Seurat object
#   5. Perform initial quality checks
#   6. Save the raw Seurat object
#
# Author: Ansh
############################################################

###############################
# Install Packages (Run Once)
###############################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "dplyr",
  "patchwork",
  "future",
  "data.table"
)

installed <- rownames(installed.packages())

for(pkg in required_packages){

  if(!(pkg %in% installed)){
    install.packages(pkg)
  }

}

###############################
# Load Libraries
###############################

library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(future)
library(data.table)

###############################
# Project Directories
###############################

project_dir <- "C:/Ovarian_CAR_T_Project"

data_dir <- file.path(project_dir, "Xenium_data")

output_dir <- file.path(project_dir, "objects")

figure_dir <- file.path(project_dir, "figures")

results_dir <- file.path(project_dir, "results")

###############################
# Create Directories
###############################

dirs <- c(output_dir,
          figure_dir,
          results_dir)

for(d in dirs){

  if(!dir.exists(d))
    dir.create(d, recursive = TRUE)

}

###############################
# Set Working Directory
###############################

setwd(project_dir)

cat("Working directory:\n")
getwd()

###############################
# Increase Memory (Windows)
###############################

options(future.globals.maxSize = 16 * 1024^3)

###############################
# Locate Xenium Matrix
###############################

matrix_file <- file.path(
  data_dir,
  "cell_feature_matrix.h5"
)

if(!file.exists(matrix_file))
  stop("cell_feature_matrix.h5 not found.")

###############################
# Read Xenium Matrix
###############################

cat("Reading expression matrix...\n")

counts <- Read10X_h5(matrix_file)

cat("Finished loading matrix.\n")

###############################
# Create Seurat Object
###############################

xenium <- CreateSeuratObject(
  counts = counts,
  project = "Ovarian_Xenium",
  min.cells = 3,
  min.features = 200
)

###############################
# Initial QC Statistics
###############################

cat("\n============================\n")
cat("Initial Seurat Summary\n")
cat("============================\n")

print(xenium)

cat("\nNumber of Cells:\n")
print(ncol(xenium))

cat("\nNumber of Genes:\n")
print(nrow(xenium))

###############################
# Metadata Inspection
###############################

head(xenium@meta.data)

summary(xenium$nFeature_RNA)

summary(xenium$nCount_RNA)

###############################
# Save Raw Object
###############################

saveRDS(
  xenium,
  file.path(
    output_dir,
    "01_xenium_raw.rds"
  )
)

cat("\nRaw Seurat object saved.\n")

###############################
# Session Information
###############################

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir,
            "sessionInfo_setup.txt")
)

cat("\nSetup complete.\n")