
# ============================================================
# Figure_4.R
# Figure 4: Marker-expression UMAPs
# Uses existing UMAP coordinates and Xenium expression data.
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

project_dir <- "C:/Ovarian_CAR_T_Project"
object_path <- file.path(project_dir, "xenium_processed.rds")
figure_dir <- file.path(project_dir, "figures")

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# Load the saved object if it is not already in memory
if (!exists("recovered_obj")) {
  recovered_obj <- readRDS(object_path)
}

DefaultAssay(recovered_obj) <- "Xenium"

# Ten marker genes
marker_genes <- c(
  "MSLN", "MUC16", "FOLR1",  # Tumor-associated
  "DCN", "LUM",               # Fibroblast / mixed stromal
  "C1QC",                     # Macrophages
  "PECAM1",                   # Endothelial
  "TRAC",                     # T cells
  "MYH11",                    # Muscle / mural
  "FOXJ1"                     # Ciliated epithelial
)

# Validate markers
missing_genes <- setdiff(
  marker_genes,
  rownames(recovered_obj[["Xenium"]])
)

if (length(missing_genes) > 0) {
  stop(
    "Missing markers: ",
    paste(missing_genes, collapse = ", ")
  )
}

if (!"umap" %in% names(recovered_obj@reductions)) {
  stop("Existing UMAP coordinates were not found.")
}

cat("All ten markers and existing UMAP coordinates validated.\n")

# Generate the expression panels with much smaller points
p_marker_grid <- FeaturePlot(
  object = recovered_obj,
  features = marker_genes,
  reduction = "umap",
  cols = c("grey90", "#08519C"),
  min.cutoff = 0,
  max.cutoff = "q99",
  keep.scale = "feature",
  order = TRUE,
  raster = TRUE,
  pt.size = 0.01,
  ncol = 5
) &
  theme_classic(base_size = 10) &
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )

# Save under NEW filenames to preserve the previous version
png_path <- file.path(
  figure_dir,
  "Figure4_marker_validation_smallpoints.png"
)

pdf_path <- file.path(
  figure_dir,
  "Figure4_marker_validation_smallpoints.pdf"
)

cat("Saving revised PNG...\n")

ggsave(
  filename = png_path,
  plot = p_marker_grid,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

cat("Saving revised PDF...\n")

ggsave(
  filename = pdf_path,
  plot = p_marker_grid,
  width = 16,
  height = 8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

# Confirm that both output files exist
for (f in c(png_path, pdf_path)) {
  if (!file.exists(f)) {
    stop("Output file was not created: ", f)
  }
  
  cat(
    "Saved: ", f,
    " (", round(file.info(f)$size / 1024), " KB)\n",
    sep = ""
  )
}

cat("Revised Figure 4 generation complete.\n")
