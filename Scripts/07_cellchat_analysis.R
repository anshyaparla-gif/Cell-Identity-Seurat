############################################################
# 07_cellchat_analysis.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   Reproduce the CellChat analysis steps captured in the
#   history file and save the results.
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
  "CellChat",
  "pheatmap"
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
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(CellChat)
library(pheatmap)

###############################
# Project directories
###############################

project_dir <- "C:/Ovarian_CAR_T_Project"
objects_dir <- file.path(project_dir, "objects")
results_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")

for (dir_path in c(objects_dir, results_dir, figure_dir)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

setwd(project_dir)

###############################
# Load Seurat object
###############################

seurat_obj <- NULL

if (exists("xenium_cellchat") && inherits(get("xenium_cellchat"), "Seurat")) {
  seurat_obj <- get("xenium_cellchat")
  message("Using xenium_cellchat from the R session.")
} else if (exists("xenium_full_umap") && inherits(get("xenium_full_umap"), "Seurat")) {
  seurat_obj <- get("xenium_full_umap")
  message("Using xenium_full_umap from the R session.")
} else if (exists("xenium_processed") && inherits(get("xenium_processed"), "Seurat")) {
  seurat_obj <- get("xenium_processed")
  message("Using xenium_processed from the R session.")
} else {
  candidate_files <- c(
    file.path(project_dir, "xenium_processed.rds"),
    file.path(project_dir, "xenium_full_umap.rds"),
    file.path(objects_dir, "04_xenium_annotated.rds")
  )

  for (f in candidate_files) {
    if (file.exists(f)) {
      seurat_obj <- readRDS(f)
      break
    }
  }
}

if (is.null(seurat_obj)) {
  stop("No suitable Seurat object was found. Load xenium_cellchat, xenium_full_umap, or xenium_processed first.")
}

###############################
# Prepare grouping metadata
###############################

if ("cellchat_type_refined" %in% colnames(seurat_obj@meta.data)) {
  group_column <- "cellchat_type_refined"
} else if ("cellchat_type" %in% colnames(seurat_obj@meta.data)) {
  group_column <- "cellchat_type"
} else if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
  group_column <- "cell_type"
} else {
  stop("No suitable cell type metadata column was found for CellChat grouping.")
}

seurat_obj$cellchat_group <- seurat_obj@meta.data[[group_column]]

###############################
# Create CellChat object
###############################

cellchat <- createCellChat(object = seurat_obj, group.by = "cellchat_group")
cellchat@DB <- CellChatDB.human

###############################
# Run CellChat pipeline from history file
###############################

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

###############################
# Reorder identities for plotting
###############################

desired_levels <- c(
  "Tumor",
  "T-cell",
  "Macrophage",
  "IFN-immune",
  "Fibroblast",
  "Endothelial",
  "Epithelial",
  "Muscle"
)

existing_levels <- unique(as.character(seurat_obj$cellchat_group))

if (all(desired_levels %in% existing_levels)) {
  cellchat@idents <- factor(cellchat@idents, levels = desired_levels)
} else {
  message("Some requested levels were not present; using existing levels instead.")
}

###############################
# Create communication plots
###############################

group_size <- as.numeric(table(cellchat@idents))

circle_plot <- netVisual_circle(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE
)

heatmap_plot <- netVisual_heatmap(
  cellchat,
  title.name = "Predicted Cell–Cell Communication Strength Among Human Ovarian Tumor Microenvironment Populations"
)

bubble_plot <- netVisual_bubble(
  cellchat,
  sources.use = "Tumor",
  targets.use = c("T-cell", "Macrophage", "Fibroblast"),
  remove.isolate = FALSE
)

###############################
# Signaling role analyses
###############################

outgoing_roles <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
incoming_roles <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

###############################
# Pathway subsets from history file
###############################

jam_cldn <- subsetCommunication(cellchat, signaling = c("JAM", "CLDN"))
vegf_paths <- subsetCommunication(cellchat, signaling = "VEGF")
cxcl_paths <- subsetCommunication(cellchat, signaling = "CXCL")

###############################
# Marker expression heatmap
###############################

meta <- cellchat@meta
if (!is.null(meta) && nrow(meta) > 0) {
  meta$cell_type <- meta$cellchat_group
  data_input <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))

  avg_expression <- sapply(
    split(colnames(data_input), meta$cell_type),
    function(cells) {
      Matrix::rowMeans(data_input[, cells, drop = FALSE])
    }
  )

  avg_expression <- as.matrix(avg_expression)
  colnames(avg_expression) <- gsub("_", "-", colnames(avg_expression))

  heatmap_genes <- c(
    "MSLN", "MUC16", "FOLR1", "EPCAM", "KRT8", "KRT18", "KRT19",
    "CD3D", "TRAC", "CD4", "CD8A", "GZMB", "GZMA", "NKG7",
    "LYZ", "CD14", "FCGR3A", "FCER1G", "CTSS",
    "IRF7", "IRF8",
    "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2",
    "PECAM1", "VWF", "CDH5", "EMCN",
    "ACTA2", "TAGLN",
    "KRT7", "KRT17"
  )

  heatmap_genes <- intersect(heatmap_genes, rownames(avg_expression))
  heatmap_matrix <- avg_expression[heatmap_genes, ]
  heatmap_matrix <- t(scale(t(heatmap_matrix)))

  heat_colors <- colorRampPalette(c("skyblue", "white", "red"))(100)

  pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = heat_colors,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = NA,
    filename = file.path(figure_dir, "cell_type_average_expression_heatmap.pdf")
  )
}

###############################
# Save outputs
###############################

saveRDS(cellchat, file.path(results_dir, "cellchat_from_history.rds"))

pdf(file.path(figure_dir, "cellchat_history_summary.pdf"), width = 10, height = 8)
print(circle_plot)
print(heatmap_plot)
print(bubble_plot)
print(outgoing_roles)
print(incoming_roles)
dev.off()

message("CellChat analysis complete.")
message("Results saved to: ", results_dir)
message("Figures saved to: ", figure_dir)
