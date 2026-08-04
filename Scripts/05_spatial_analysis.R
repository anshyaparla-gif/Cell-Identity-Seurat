############################################################
# 05_spatial_analysis.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   1. Load the annotated Seurat object
#   2. Analyze spatial relationships between tumor and immune cells
#   3. Save distance summaries
#
############################################################

###############################
# Install / load packages
###############################

required_packages <- c("Seurat", "SeuratObject", "ggplot2", "dplyr")
installed <- rownames(installed.packages())

for(pkg in required_packages){
  if(!(pkg %in% installed)){
    install.packages(pkg)
  }
}

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)

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
# Spatial coordinates
###############################

spatial_coords <- NULL

if("x" %in% colnames(xenium@meta.data) && "y" %in% colnames(xenium@meta.data)){
  spatial_coords <- xenium@meta.data[, c("x", "y")]
} else if(methods::is(xenium, "Seurat")){
  tryCatch({
    spatial_coords <- GetTissueCoordinates(xenium)
  }, error = function(e) {
    NULL
  })
}

if(is.null(spatial_coords)){
  warning("No spatial coordinates found. Spatial analysis skipped.")
} else {
  spatial_coords$cell_id <- rownames(spatial_coords)
  colnames(spatial_coords)[1:2] <- c("x", "y")

  tumor_cells <- rownames(xenium@meta.data)[xenium$cell_type %in% c("Tumor", "Epithelial")]
  immune_cells <- rownames(xenium@meta.data)[xenium$cell_type %in% c("T-cell", "Macrophage", "IFN-immune")]

  if(length(tumor_cells) > 0 && length(immune_cells) > 0){
    tumor_coords <- spatial_coords[spatial_coords$cell_id %in% tumor_cells, ]
    immune_coords <- spatial_coords[spatial_coords$cell_id %in% immune_cells, ]

    dist_mat <- as.matrix(dist(rbind(tumor_coords[, c("x", "y")], immune_coords[, c("x", "y")])))[1:nrow(tumor_coords), (nrow(tumor_coords) + 1):nrow(spatial_coords)]
    dist_df <- as.data.frame(dist_mat)
    dist_df$cell_id <- tumor_coords$cell_id
    dist_df <- tidyr::pivot_longer(dist_df, cols = -cell_id, names_to = "neighbor_cell", values_to = "distance_um")

    write.table(
      dist_df,
      file.path(results_dir, "05_spatial_distances.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    summary_stats <- dist_df %>%
      group_by(cell_id) %>%
      summarise(min_distance_um = min(distance_um, na.rm = TRUE))

    write.table(
      summary_stats,
      file.path(results_dir, "05_spatial_distance_summary.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}

cat("Spatial analysis complete.\n")
