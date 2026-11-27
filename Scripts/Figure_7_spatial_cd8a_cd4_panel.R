# ============================================================
# Figure 7 replacement line figure
# Spatial maps for tumor-antigen score, CD8A, and CD4.
# The figure uses the Xenium FOV coordinates from recovered_obj
# and highlights expression-positive cells on the same spatial
# frame of reference.
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

if (exists("recovered_obj", inherits = TRUE)) {
  obj <- recovered_obj
} else {
  if (!file.exists(object_path)) {
    stop("Saved Seurat object not found: ", object_path)
  }
  obj <- readRDS(object_path)
}

stopifnot(inherits(obj, "Seurat"))

DefaultAssay(obj) <- "Xenium"

if (!"fov" %in% SeuratObject::Images(obj)) {
  stop("Spatial image 'fov' was not found in the Seurat object.")
}

tumor_markers <- c("MSLN", "MUC16", "EPCAM", "FOLR1")
panel_features <- c(tumor_markers, "CD8A", "CD4")

missing_features <- setdiff(panel_features, rownames(obj[["Xenium"]]))
if (length(missing_features) > 0) {
  stop(
    "Missing features in Xenium assay: ",
    paste(missing_features, collapse = ", ")
  )
}

coords <- Seurat::GetTissueCoordinates(obj, image = "fov")

required_coord_cols <- c("x", "y", "cell")
if (!all(required_coord_cols %in% colnames(coords))) {
  stop(
    "Expected coordinate columns are missing. Found: ",
    paste(colnames(coords), collapse = ", ")
  )
}

expr <- FetchData(obj, vars = panel_features, layer = "data")
expr$cell <- rownames(expr)

coord_index <- match(as.character(coords$cell), expr$cell)

plot_data <- data.frame(
  cell = as.character(coords$cell),
  x = coords$x,
  y = coords$y,
  tumor_score = rowMeans(
    expr[coord_index, tumor_markers, drop = FALSE],
    na.rm = TRUE
  ),
  cd8a = expr[coord_index, "CD8A"],
  cd4 = expr[coord_index, "CD4"],
  stringsAsFactors = FALSE
)

plot_data <- plot_data[complete.cases(plot_data), ]

make_panel <- function(data, value_col, title, legend_title) {
  positive_data <- data[data[[value_col]] > 0, , drop = FALSE]

  high_color <- switch(
    value_col,
    tumor_score = "#D73027",
    cd8a = "#2C7BB6",
    cd4 = "#31A354",
    "#08519C"
  )

  ggplot(data, aes(x = x, y = y)) +
    geom_point(color = "grey88", size = 0.03, alpha = 0.18) +
    geom_point(
      data = positive_data,
      aes(color = .data[[value_col]]),
      size = 0.06,
      alpha = 0.95
    ) +
    scale_color_gradient(
      low = "#D9EAF7",
      high = high_color,
      name = legend_title
    ) +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 9),
      legend.position = "right",
      plot.margin = margin(8, 8, 8, 8)
    )
}

p_tumor <- make_panel(
  plot_data,
  value_col = "tumor_score",
  title = "Tumor antigen composite",
  legend_title = "Mean log-expression\n(MSLN, MUC16, EPCAM, FOLR1)"
)

p_cd8a <- make_panel(
  plot_data,
  value_col = "cd8a",
  title = "CD8A",
  legend_title = "CD8A log-expression"
)

p_cd4 <- make_panel(
  plot_data,
  value_col = "cd4",
  title = "CD4",
  legend_title = "CD4 log-expression"
)

panel_plot <- p_tumor + p_cd8a + p_cd4 +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

png_path <- file.path(figure_dir, "Figure7_spatial_cd8a_cd4_panel.png")
pdf_path <- file.path(figure_dir, "Figure7_spatial_cd8a_cd4_panel.pdf")

ggsave(
  filename = png_path,
  plot = panel_plot,
  width = 15,
  height = 5.5,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = pdf_path,
  plot = panel_plot,
  width = 15,
  height = 5.5,
  bg = "white",
  limitsize = FALSE
)

for (f in c(png_path, pdf_path)) {
  if (!file.exists(f)) {
    stop("Output file was not created: ", f)
  }
}

cat("Saved:", png_path, "\n")
cat("Saved:", pdf_path, "\n")

