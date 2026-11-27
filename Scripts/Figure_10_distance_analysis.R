# ============================================================
# Figure 10 distance analysis
# Compares tumor antigen upregulation with distance to T cells.
# Produces two outputs:
#   1) Continuous correlation scatterplot
#   2) High-vs-low tumor antigen distance comparison
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

project_dir <- "C:/Ovarian_CAR_T_Project"
object_path <- file.path(project_dir, "xenium_processed.rds")
figure_dir <- file.path(project_dir, "figures")
results_dir <- file.path(project_dir, "results")

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
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

required_features <- c("MSLN", "MUC16", "EPCAM", "FOLR1")
missing_features <- setdiff(required_features, rownames(obj[["Xenium"]]))
if (length(missing_features) > 0) {
  stop(
    "Missing features in Xenium assay: ",
    paste(missing_features, collapse = ", ")
  )
}

if (!"cell_category_display" %in% colnames(obj@meta.data)) {
  stop("Metadata column 'cell_category_display' is missing.")
}

coords <- Seurat::GetTissueCoordinates(obj, image = "fov")
required_coord_cols <- c("x", "y", "cell")
if (!all(required_coord_cols %in% colnames(coords))) {
  stop(
    "Expected coordinate columns are missing. Found: ",
    paste(colnames(coords), collapse = ", ")
  )
}

coord_df <- data.frame(
  cell = as.character(coords$cell),
  x = coords$x,
  y = coords$y,
  category = as.character(obj@meta.data[as.character(coords$cell), "cell_category_display"]),
  stringsAsFactors = FALSE
)

coord_df <- coord_df[!is.na(coord_df$category), ]

tumor_cells <- coord_df |>
  filter(category == "Tumor")

t_cells <- coord_df |>
  filter(category == "T cells")

if (nrow(tumor_cells) == 0) {
  stop("No tumor cells were found using cell_category_display == 'Tumor'.")
}

if (nrow(t_cells) == 0) {
  stop("No T cells were found using cell_category_display == 'T cells'.")
}

tumor_expr <- FetchData(obj, vars = required_features, layer = "data")
tumor_expr$cell <- rownames(tumor_expr)

tumor_df <- tumor_cells |>
  inner_join(tumor_expr, by = "cell") |>
  mutate(
    antigen_score = rowMeans(across(all_of(required_features))),
    antigen_group = if_else(
      antigen_score >= median(antigen_score, na.rm = TRUE),
      "High antigen",
      "Low antigen"
    )
  )

nearest_distances <- function(query_xy, ref_xy) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for nearest-neighbor distance calculation.")
  }
  nn <- RANN::nn2(data = ref_xy, query = query_xy, k = 1)
  as.numeric(nn$nn.dists[, 1])
}

t_ref <- as.matrix(t_cells[, c("x", "y")])
tumor_query <- as.matrix(tumor_df[, c("x", "y")])
tumor_df$distance_to_nearest_T <- nearest_distances(tumor_query, t_ref)

spearman_test <- cor.test(
  tumor_df$antigen_score,
  tumor_df$distance_to_nearest_T,
  method = "spearman",
  exact = FALSE
)

scatter_plot <- ggplot(tumor_df, aes(x = antigen_score, y = distance_to_nearest_T)) +
  geom_point(alpha = 0.25, size = 0.5, color = "#D73027") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  labs(
    title = "Tumor antigen score vs. distance to nearest T cell",
    subtitle = paste0(
      "Spearman rho = ", round(unname(spearman_test$estimate), 3),
      ", p = ", format.pval(spearman_test$p.value, digits = 3, eps = 0.001)
    ),
    x = "Tumor antigen score (mean log-expression of MSLN, MUC16, EPCAM, FOLR1)",
    y = "Distance to nearest T cell"
  ) +
  theme_classic(base_size = 11)

group_plot <- ggplot(tumor_df, aes(x = antigen_group, y = distance_to_nearest_T, fill = antigen_group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.55, alpha = 0.8) +
  geom_jitter(width = 0.14, alpha = 0.2, size = 0.4) +
  scale_fill_manual(values = c("Low antigen" = "#9ECAE1", "High antigen" = "#FC9272")) +
  labs(
    title = "Distance to nearest T cell by tumor antigen group",
    subtitle = "Tumor cells split at the median tumor antigen score",
    x = NULL,
    y = "Distance to nearest T cell"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

scatter_png <- file.path(figure_dir, "Figure10_antigen_distance_scatter.png")
scatter_pdf <- file.path(figure_dir, "Figure10_antigen_distance_scatter.pdf")
group_png <- file.path(figure_dir, "Figure10_antigen_distance_high_low.png")
group_pdf <- file.path(figure_dir, "Figure10_antigen_distance_high_low.pdf")

ggsave(scatter_png, scatter_plot, width = 8.5, height = 6, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(scatter_pdf, scatter_plot, width = 8.5, height = 6, bg = "white", limitsize = FALSE)
ggsave(group_png, group_plot, width = 7.2, height = 5.8, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(group_pdf, group_plot, width = 7.2, height = 5.8, bg = "white", limitsize = FALSE)

summary_path <- file.path(results_dir, "Figure10_distance_summary.csv")
write.csv(
  data.frame(
    spearman_rho = unname(spearman_test$estimate),
    spearman_p_value = spearman_test$p.value,
    n_tumor_cells = nrow(tumor_df),
    n_t_cells = nrow(t_cells)
  ),
  summary_path,
  row.names = FALSE
)

for (f in c(scatter_png, scatter_pdf, group_png, group_pdf, summary_path)) {
  if (!file.exists(f)) {
    stop("Output file was not created: ", f)
  }
}

cat("Saved:\n")
cat(" - ", scatter_png, "\n", sep = "")
cat(" - ", scatter_pdf, "\n", sep = "")
cat(" - ", group_png, "\n", sep = "")
cat(" - ", group_pdf, "\n", sep = "")
cat(" - ", summary_path, "\n", sep = "")

