# ============================================================
# Figure 7 immune-marker violin plots
# Final paper version: six-gene violin plot analysis used to
# support the T-cell marker interpretation.
#
# This script keeps the figure aligned with the final paper
# numbering and saves outputs as Figure7_immune_marker_violins.*
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

project_dir <- "C:/Ovarian_CAR_T_Project"
figure_dir <- file.path(project_dir, "figures")
output_dir <- file.path(project_dir, "results")

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (exists("recovered_obj", inherits = TRUE)) {
  obj <- recovered_obj
} else {
  object_path <- file.path(project_dir, "xenium_processed.rds")
  if (!file.exists(object_path)) {
    stop("Saved Seurat object not found: ", object_path)
  }
  obj <- readRDS(object_path)
}

stopifnot(inherits(obj, "Seurat"))

DefaultAssay(obj) <- "Xenium"

# Six-gene panel used for the final Figure 7 layout.
figure7_genes <- c(
  "CD3D",
  "TRAC",
  "CD4",
  "CD8A",
  "CD8B",
  "FOXP3"
)

required_metadata <- "cell_category_display"
if (!required_metadata %in% colnames(obj[[]])) {
  stop("Metadata column '", required_metadata, "' was not found in the Seurat object.")
}

missing_genes <- setdiff(figure7_genes, rownames(obj[["Xenium"]]))
if (length(missing_genes) > 0) {
  stop("Missing Figure 7 genes in Xenium assay: ", paste(missing_genes, collapse = ", "))
}

figure7_data <- FetchData(
  obj,
  vars = c(required_metadata, figure7_genes),
  layer = "data"
) |>
  filter(cell_category_display %in% c("T cells", "Macrophages", "Tumor")) |>
  mutate(
    cell_category_display = factor(
      cell_category_display,
      levels = c("T cells", "Macrophages", "Tumor")
    )
  )

figure7_long <- figure7_data |>
  pivot_longer(
    cols = all_of(figure7_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) |>
  mutate(
    Gene = factor(Gene, levels = figure7_genes)
  )

figure7_summary <- figure7_long |>
  group_by(Gene, cell_category_display) |>
  summarise(
    Total_Cells = n(),
    Positive_Cells = sum(Expression > 0),
    Percent_Detected = mean(Expression > 0) * 100,
    Mean_All_Cells = mean(Expression),
    Mean_Positive_Cells = ifelse(
      sum(Expression > 0) > 0,
      mean(Expression[Expression > 0]),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  arrange(Gene, cell_category_display)

write.csv(
  figure7_summary,
  file.path(output_dir, "Figure7_expression_summary.csv"),
  row.names = FALSE
)

figure7_plot <- ggplot(
  figure7_long,
  aes(
    x = Gene,
    y = Expression,
    fill = cell_category_display
  )
) +
  geom_violin(
    position = position_dodge(width = 0.85),
    scale = "width",
    trim = TRUE,
    alpha = 0.85,
    color = "black",
    linewidth = 0.25
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.85),
    width = 0.12,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "T cells" = "#7B3294",
      "Macrophages" = "#008837",
      "Tumor" = "#E66101"
    )
  ) +
  labs(
    title = "Figure 7. Immune Marker Expression Across Cell Types",
    x = NULL,
    y = "Normalized expression",
    fill = "Cell category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

png_path <- file.path(figure_dir, "Figure7_immune_marker_violins.png")
pdf_path <- file.path(figure_dir, "Figure7_immune_marker_violins.pdf")

ggsave(
  filename = png_path,
  plot = figure7_plot,
  width = 11,
  height = 6,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = pdf_path,
  plot = figure7_plot,
  width = 11,
  height = 6,
  bg = "white",
  limitsize = FALSE
)

for (f in c(png_path, pdf_path, file.path(output_dir, "Figure7_expression_summary.csv"))) {
  if (!file.exists(f)) {
    stop("Output file was not created: ", f)
  }
}

cat("Saved:\n")
cat(" - ", png_path, "\n", sep = "")
cat(" - ", pdf_path, "\n", sep = "")
cat(" - ", file.path(output_dir, "Figure7_expression_summary.csv"), "\n", sep = "")

