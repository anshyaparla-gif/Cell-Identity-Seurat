############################################################
# 06_car_t_analysis.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   1. Load the annotated object
#   2. Summarize expression of CAR-T related and tumor-associated genes
#   3. Save a results table for downstream interpretation
#
############################################################

###############################
# Install / load packages
###############################

required_packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2")
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
# Load object
###############################

object_file <- file.path(objects_dir, "04_xenium_annotated.rds")

if(!file.exists(object_file)){
  object_file <- file.path(project_dir, "xenium_processed.rds")
}

if(!file.exists(object_file)){
  stop("Annotated object not found. Run 04_cell_annotation.R first.")
}

xenium <- readRDS(object_file)

###############################
# Genes of interest
###############################

genes_of_interest <- c(
  "CXCL12",
  "CXCR4",
  "VEGFA",
  "KDR",
  "CD3D",
  "CD4",
  "CD8A",
  "NKG7",
  "EPCAM",
  "MSLN",
  "MUC16",
  "FOLR1",
  "PECAM1",
  "VWF"
)

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(xenium)]

###############################
# Average expression per cell type
###############################

if(length(genes_of_interest) > 0){
  avg_exp <- AverageExpression(xenium, features = genes_of_interest, group.by = "cell_type", slot = "data")
  avg_exp <- as.data.frame(avg_exp$RNA)
  write.table(avg_exp, file.path(results_dir, "06_car_t_target_expression.tsv"), sep = "\t", quote = FALSE)
}

cat("CAR-T related expression summary saved.\n")
