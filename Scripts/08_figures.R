############################################################
# 08_figures.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   Recreate the figure-generation workflow from the history file
#   and save publication-style outputs.
############################################################

###############################
# Install / load packages
###############################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "dplyr"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed)) {
    message("Installing ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(dplyr)

###############################
# Project directories
###############################

project_dir <- "C:/Ovarian_CAR_T_Project"
objects_dir <- file.path(project_dir, "objects")
results_dir <- file.path(project_dir, "results")
fig_dir <- file.path(project_dir, "figures")

for (dir_path in c(objects_dir, results_dir, fig_dir)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

setwd(project_dir)

###############################
# Load object
###############################

object_file <- file.path(objects_dir, "04_xenium_annotated.rds")

if (!file.exists(object_file)) {
  object_file <- file.path(project_dir, "xenium_processed.rds")
}

if (!file.exists(object_file)) {
  stop("Annotated object not found. Run 04_cell_annotation.R first.")
}

xenium <- readRDS(object_file)

###############################
# Figure 1: UMAP / cluster overview
###############################

if ("seurat_clusters" %in% names(xenium@meta.data)) {
  umap_plot <- DimPlot(xenium, group.by = "seurat_clusters", label = TRUE) +
    ggtitle("UMAP of Ovarian Xenium Cells")
  ggsave(file.path(fig_dir, "08_umap_clusters.pdf"), umap_plot, width = 8, height = 6)
}

###############################
# Figure 2: Marker expression summary
###############################

marker_genes <- c("MSLN", "MUC16", "EPCAM", "CD3D", "TRAC", "LYZ", "COL1A1", "PECAM1")
marker_genes <- marker_genes[marker_genes %in% rownames(xenium)]

if (length(marker_genes) > 0) {
  marker_plot <- DotPlot(xenium, features = marker_genes, group.by = "cell_type") +
    RotatedAxis() +
    ggtitle("Marker Expression by Cell Type")
  ggsave(file.path(fig_dir, "08_marker_expression_dotplot.pdf"), marker_plot, width = 12, height = 7)
}

###############################
# Figure 3: Cell type composition
###############################

if ("cell_type" %in% colnames(xenium@meta.data)) {
  comp_df <- as.data.frame(table(xenium$cell_type))
  colnames(comp_df) <- c("cell_type", "count")

  comp_plot <- ggplot(comp_df, aes(x = reorder(cell_type, -count), y = count)) +
    geom_col(fill = "steelblue") +
    theme_classic() +
    labs(title = "Cell Type Composition", x = "Cell Type", y = "Number of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(fig_dir, "08_cell_type_composition.pdf"), comp_plot, width = 8, height = 6)
}

###############################
# Figure 4: CellChat summary plots
###############################

cellchat_file <- file.path(results_dir, "cellchat_from_history.rds")
if (file.exists(cellchat_file)) {
  cellchat <- readRDS(cellchat_file)

  group_size <- as.numeric(table(cellchat@idents))
  circle_plot <- netVisual_circle(
    cellchat@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE
  )

  heatmap_plot <- netVisual_heatmap(
    cellchat,
    title.name = "Predicted Cell–Cell Communication Strength"
  )

  pdf(file.path(fig_dir, "08_cellchat_summary_plots.pdf"), width = 10, height = 8)
  print(circle_plot)
  print(heatmap_plot)
  dev.off()
}

message("Figure generation complete. Outputs saved to: ", fig_dir)
