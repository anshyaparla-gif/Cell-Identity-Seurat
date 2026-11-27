# ============================================================
# Figure 11: CellChat circle plot
# Final cleaned version for the paper.
# Uses major cell-type groups, excludes unresolved labels,
# and downsamples per group so the workflow remains runnable.
# ============================================================

library(Seurat)
library(CellChat)
library(dplyr)

project_dir <- "C:/Ovarian_CAR_T_Project"
object_path <- file.path(project_dir, "objects", "04_xenium_annotated.rds")
figure_dir <- file.path(project_dir, "figures")
results_dir <- file.path(project_dir, "results")

if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

if (!file.exists(object_path)) {
  stop("Annotated Seurat object not found: ", object_path)
}

obj <- readRDS(object_path)
stopifnot(inherits(obj, "Seurat"))

if (!"Xenium" %in% Assays(obj)) {
  stop("Xenium assay was not found in the annotated Seurat object.")
}

metadata <- obj@meta.data

group_col <- if ("cellchat_type_refined" %in% colnames(metadata)) {
  "cellchat_type_refined"
} else if ("cell_type_refined" %in% colnames(metadata)) {
  "cell_type_refined"
} else {
  stop("Neither 'cellchat_type_refined' nor 'cell_type_refined' was found.")
}

normalize_group <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("_", " ", x)
  x <- gsub("[[:space:]]+", " ", x)

  x[grepl("^(t ?cells?|t-cell(s)?)$", x, ignore.case = TRUE)] <- "T cells"
  x[grepl("macrophage", x, ignore.case = TRUE)] <- "Macrophages"
  x[grepl("muscle|mural|pericyte", x, ignore.case = TRUE)] <- "Muscle/Mural"
  x[grepl("fibroblast|stromal", x, ignore.case = TRUE)] <- "Fibroblast"
  x[grepl("endothelial|pecam1|cdh5|vwf", x, ignore.case = TRUE)] <- "Endothelial"
  x[grepl("epithelial.*tumor|tumor.*epithelial|epithelial tumor|tumor", x, ignore.case = TRUE)] <- "Tumor"
  x[grepl("epithelial", x, ignore.case = TRUE)] <- "Epithelial"
  x[grepl("ifn", x, ignore.case = TRUE)] <- "IFN-immune"
  x[grepl("misc|unresolved|low[ ._-]?feature|low[ ._-]?confidence|unknown", x, ignore.case = TRUE)] <- NA_character_

  x
}

cellchat_groups <- normalize_group(metadata[[group_col]])
allowed_groups <- c(
  "Tumor",
  "T cells",
  "Macrophages",
  "Fibroblast",
  "Endothelial",
  "Epithelial",
  "Muscle/Mural"
)
cellchat_groups[!cellchat_groups %in% allowed_groups] <- NA_character_

obj$cellchat_group_clean <- factor(cellchat_groups, levels = allowed_groups)
keep_cells <- !is.na(obj$cellchat_group_clean)
obj <- subset(obj, cells = colnames(obj)[keep_cells])
obj$cellchat_group_clean <- droplevels(obj$cellchat_group_clean)

# Downsample each retained group so CellChat can finish in a practical time.
set.seed(7381)
max_cells_per_group <- 50
cell_ids_by_group <- split(colnames(obj), obj$cellchat_group_clean)
sampled_cells <- unlist(lapply(cell_ids_by_group, function(ids) {
  if (length(ids) > max_cells_per_group) {
    sample(ids, max_cells_per_group)
  } else {
    ids
  }
}), use.names = FALSE)

obj <- subset(obj, cells = sampled_cells)
obj$cellchat_group_clean <- droplevels(obj$cellchat_group_clean)

DefaultAssay(obj) <- "Xenium"

data.input <- GetAssayData(obj, assay = "Xenium", slot = "data")
meta <- data.frame(cell_type = obj$cellchat_group_clean, row.names = colnames(obj))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

main_png <- file.path(figure_dir, "Figure11_cellchat_circle_major_groups.png")
main_pdf <- file.path(figure_dir, "Figure11_cellchat_circle_major_groups.pdf")

png(main_png, width = 8.5, height = 8, units = "in", res = 300, bg = "white")
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Predicted cell-cell communication among major cell types"
)
dev.off()

pdf(main_pdf, width = 8.5, height = 8, bg = "white")
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Predicted cell-cell communication among major cell types"
)
dev.off()

group_counts <- as.data.frame(table(cellchat@idents))
colnames(group_counts) <- c("cell_type", "n_cells")
write.csv(group_counts, file.path(results_dir, "Figure11_cellchat_group_counts.csv"), row.names = FALSE)

for (f in c(main_png, main_pdf, file.path(results_dir, "Figure11_cellchat_group_counts.csv"))) {
  if (!file.exists(f)) {
    stop("Output file was not created: ", f)
  }
}

cat("Saved:\n")
cat(" - ", main_png, "\n", sep = "")
cat(" - ", main_pdf, "\n", sep = "")
cat(" - ", file.path(results_dir, "Figure11_cellchat_group_counts.csv"), "\n", sep = "")
