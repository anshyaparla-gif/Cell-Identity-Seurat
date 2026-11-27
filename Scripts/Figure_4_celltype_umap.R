
library(Seurat)
library(ggplot2)

figure_dir <- "C:/Ovarian_CAR_T_Project/figures"

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# ============================================================
# Ovarian Cancer Xenium: UMAP with display categories
# Reuses the existing UMAP coordinates in recovered_obj.
# Does not recluster or overwrite original cell-type labels.
# ============================================================

# 1. Map existing clusters to display categories.
#    Clusters 3 and 4 are split provisionally.
cluster_to_category <- c(
  "0"  = "Tumor",
  "1"  = "Muscle/Mural",
  "2"  = "Tumor",
  "3"  = "Fibroblast",
  "4"  = "Mixed Stromal",
  "5"  = "Macrophages",
  "6"  = "Endothelial",
  "7"  = "Muscle/Mural",
  "8"  = "T cells",
  "9"  = "Miscellaneous",
  "10" = "Tumor",
  "11" = "Epithelial",
  "12" = "Low-Feature",
  "13" = "Low-Feature"
)

cluster_ids <- as.character(recovered_obj$seurat_clusters)

recovered_obj$cell_category_display <- unname(
  cluster_to_category[cluster_ids]
)

# 2. Check that all cells received a category.
print(table(
  recovered_obj$cell_category_display,
  useNA = "ifany"
))

if (anyNA(recovered_obj$cell_category_display)) {
  stop("Some cells have no display category. Check cluster IDs.")
}

# 3. Optional: combine clusters 3 and 4 into one category.
#    Change to TRUE if you prefer the original broader grouping.
combine_stromal <- FALSE

if (combine_stromal) {
  recovered_obj$cell_category_display[
    recovered_obj$cell_category_display == "Mixed Stromal"
  ] <- "Fibroblast/Stromal"
  recovered_obj$cell_category_display[
    recovered_obj$cell_category_display == "Fibroblast"
  ] <- "Fibroblast/Stromal"
}

# 4. Plot using the EXISTING UMAP coordinates.
#    No text labels on the embedding; retain the side legend.
p_categories <- DimPlot(
  recovered_obj,
  reduction = "umap",
  group.by = "cell_category_display",
  label = FALSE,
  raster = TRUE
) +
  ggtitle("UMAP of Cell and Tissue Types") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

print(p_categories)

# 5. Save at the dimensions you used previously.
ggsave(
  filename = file.path(figure_dir, "Figure4_celltype_umap.png"),
  plot = p_categories,
  width = 10,
  height = 8,
  dpi = 300
)

# 6. Confirm that the file was created.
print(file.exists(
  "C:/Ovarian_CAR_T_Project/Figure4_celltype_umap.png"
))



# ============================================================
# Final UMAP category mapping
# Combines Miscellaneous and Low-Feature clusters
# ============================================================

cluster_to_category <- c(
  "0"  = "Tumor",
  "1"  = "Muscle / Mural",
  "2"  = "Tumor",
  "3"  = "Fibroblast",
  "4"  = "Mixed Stromal",
  "5"  = "Macrophages",
  "6"  = "Endothelial",
  "7"  = "Muscle / Mural",
  "8"  = "T cells",
  "9"  = "Miscellaneous / Unresolved",
  "10" = "Tumor",
  "11" = "Ciliated Epithelial",
  "12" = "Miscellaneous / Unresolved",
  "13" = "Miscellaneous / Unresolved"
)

cluster_ids <- as.character(recovered_obj$seurat_clusters)

recovered_obj$cell_category_display <- unname(
  cluster_to_category[cluster_ids]
)

# Check the category assignments
print(table(
  recovered_obj$cell_category_display,
  useNA = "ifany"
))

# Confirm no cells were left without a category
stopifnot(!anyNA(recovered_obj$cell_category_display))

# Generate UMAP using the existing coordinates
p_categories <- DimPlot(
  recovered_obj,
  reduction = "umap",
  group.by = "cell_category_display",
  label = FALSE,
  raster = TRUE
) +
  ggtitle("UMAP of Cell and Tissue Types") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

print(p_categories)

# Save final figure
ggsave(
  filename = file.path(figure_dir, "Figure4_celltype_umap.png"),
  plot = p_categories,
  width = 10,
  height = 8,
  dpi = 300
)