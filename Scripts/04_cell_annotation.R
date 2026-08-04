############################################################
# 04_cell_annotation.R
#
# Project: Ovarian Cancer Spatial Transcriptomics Analysis
# Purpose:
#   1. Load the preprocessed Seurat object
#   2. Annotate cells into broad biological populations
#   3. Save the annotated object
#
############################################################

###############################
# Install / load packages
###############################

required_packages <- c("Seurat", "SeuratObject", "ggplot2")
installed <- rownames(installed.packages())

for(pkg in required_packages){
  if(!(pkg %in% installed)){
    install.packages(pkg)
  }
}

library(Seurat)
library(SeuratObject)
library(ggplot2)

###############################
# Project directories
###############################

project_dir <- "C:/Ovarian_CAR_T_Project"
objects_dir <- file.path(project_dir, "objects")
results_dir <- file.path(project_dir, "results")

for(d in c(objects_dir, results_dir)){
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
# Annotation workflow
###############################

if("cellchat_type_refined" %in% colnames(xenium@meta.data)){
  xenium$cell_type <- xenium$cellchat_type_refined
} else if("cellchat_type" %in% colnames(xenium@meta.data)){
  xenium$cell_type <- xenium$cellchat_type
} else {
  marker_map <- list(
    Tumor = c("MSLN", "MUC16", "EPCAM", "FOLR1"),
    T_cell = c("CD3D", "CD3E", "TRAC", "CD8A", "CD4"),
    Macrophage = c("LYZ", "CD14", "FCGR3A", "CTSS"),
    Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "COL6A1"),
    Endothelial = c("PECAM1", "VWF", "CDH5", "EMCN"),
    Epithelial = c("KRT8", "KRT18", "KRT19", "KRT7"),
    IFN_immune = c("IRF7", "IRF8", "NKG7")
  )

  scores <- lapply(marker_map, function(genes) {
    genes <- genes[genes %in% rownames(xenium)]
    if(length(genes) == 0) return(NULL)
    AddModuleScore(xenium, features = list(genes), name = "score")[["score1"]]
  })

  score_df <- do.call(cbind, scores)
  colnames(score_df) <- names(marker_map)
  xenium$cell_type <- colnames(score_df)[max.col(score_df)]
}

xenium$cell_type <- gsub("_", "-", xenium$cell_type)

###############################
# Save annotated object
###############################

saveRDS(xenium, file.path(objects_dir, "04_xenium_annotated.rds"))
saveRDS(xenium, file.path(project_dir, "xenium_processed.rds"))

write.table(
  as.data.frame(table(xenium$cell_type)),
  file.path(results_dir, "04_cell_type_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Annotation complete and saved.\n")
