
# ============================================================
# Figure 3: Spatial Distribution of Cell Categories
# Project: Human Ovarian Cancer Xenium Spatial Transcriptomics
#
# Features:
#   - Timestamped progress logs
#   - Elapsed-time reporting
#   - Cell and coordinate validation
#   - PNG saved before PDF
#   - Raster-based PDF to avoid hundreds of thousands of
#     individual vector points
#   - Does not modify the saved Seurat object
# ============================================================

options(warn = 1)

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

png_path <- file.path(
  figure_dir, "Figure3_spatial_distribution.png"
)
pdf_path <- file.path(
  figure_dir, "Figure3_spatial_distribution.pdf"
)
log_path <- file.path(
  results_dir, "Figure3_progress.log"
)

# ------------------------------------------------------------
# Progress logging
# ------------------------------------------------------------

start_time <- Sys.time()

log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", timestamp, "] ", message)
  
  cat(line, "\n")
  flush.console()
  
  cat(line, "\n", file = log_path, append = TRUE)
}

run_step <- function(step_name, code) {
  log_message(paste0("START: ", step_name))
  step_start <- Sys.time()
  
  result <- tryCatch(
    force(code),
    error = function(e) {
      log_message(
        paste0("FAILED: ", step_name, " | ", conditionMessage(e))
      )
      stop(e)
    }
  )
  
  elapsed <- as.numeric(
    difftime(Sys.time(), step_start, units = "secs")
  )
  
  log_message(
    sprintf("COMPLETE: %s | Elapsed: %.1f seconds",
            step_name, elapsed)
  )
  
  result
}

# Start a fresh log for this run.
writeLines(
  paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] Figure 3 run started."
  ),
  con = log_path
)

log_message("Figure 3 script initialized.")

# ------------------------------------------------------------
# Step 1: Validate packages and files
# ------------------------------------------------------------

run_step("Validate required packages and saved object", {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed.")
  }
  
  if (!file.exists(object_path)) {
    stop("Saved Seurat object not found: ", object_path)
  }
  
  log_message(paste("Object file size:",
                    round(file.info(object_path)$size / 1e9, 2),
                    "GB"))
})

# ------------------------------------------------------------
# Step 2: Load the object
# ------------------------------------------------------------

if (exists("recovered_obj", envir = .GlobalEnv)) {
  log_message(
    "Using existing recovered_obj in memory; skipping RDS load."
  )
} else {
  recovered_obj <- run_step("Load saved Seurat object", {
    readRDS(object_path)
  })
}

log_message(
  paste("Object dimensions:",
        nrow(recovered_obj), "features x",
        ncol(recovered_obj), "cells")
)

# ------------------------------------------------------------
# Step 3: Validate spatial image and metadata
# ------------------------------------------------------------

run_step("Validate spatial image and metadata", {
  if (!"fov" %in% SeuratObject::Images(recovered_obj)) {
    stop("Spatial image 'fov' was not found in the object.")
  }
  
  if (!"cell_category_display" %in%
      colnames(recovered_obj@meta.data)) {
    stop("Metadata column 'cell_category_display' is missing.")
  }
  
  log_message("Spatial image 'fov' and category metadata found.")
})

# ------------------------------------------------------------
# Step 4: Extract tissue coordinates
# ------------------------------------------------------------

coords <- run_step("Extract Xenium tissue coordinates", {
  Seurat::GetTissueCoordinates(
    recovered_obj,
    image = "fov"
  )
})

log_message(paste("Coordinate rows:", nrow(coords)))
log_message(
  paste("Coordinate columns:", paste(colnames(coords), collapse = ", "))
)

required_coord_cols <- c("x", "y", "cell")

if (!all(required_coord_cols %in% colnames(coords))) {
  stop(
    "Expected coordinate columns are missing. Found: ",
    paste(colnames(coords), collapse = ", ")
  )
}

# ------------------------------------------------------------
# Step 5: Match coordinates to metadata by cell ID
# ------------------------------------------------------------

coords_df <- run_step("Match coordinates to cell metadata", {
  data.frame(
    cell = as.character(coords$cell),
    x = coords$x,
    y = coords$y,
    stringsAsFactors = FALSE
  )
})

metadata <- recovered_obj@meta.data

metadata_df <- data.frame(
  cell = rownames(metadata),
  category = as.character(metadata$cell_category_display),
  stringsAsFactors = FALSE
)

match_index <- run_step("Match cell identifiers", {
  match(coords_df$cell, metadata_df$cell)
})

coords_df$category <- metadata_df$category[match_index]

matched_count <- sum(!is.na(match_index))
unmatched_count <- sum(is.na(match_index))

log_message(paste("Matched cells:", matched_count))
log_message(paste("Unmatched coordinate rows:", unmatched_count))

if (matched_count == 0) {
  stop("No coordinate cell IDs matched the metadata.")
}

if (unmatched_count > 0) {
  log_message(
    "WARNING: Some coordinate rows lack matching metadata."
  )
}

# ------------------------------------------------------------
# Step 6: Standardize categories and validate colors
# ------------------------------------------------------------

run_step("Standardize cell categories", {
  unresolved_pattern <- grepl(
    "miscellaneous|unresolved|low.feature|low.confidence",
    coords_df$category,
    ignore.case = TRUE
  )
  
  coords_df$category[
    is.na(coords_df$category) |
      coords_df$category == "" |
      unresolved_pattern
  ] <- "Unresolved / low-confidence"
  
  expected_categories <- c(
    "Tumor",
    "Fibroblast",
    "Mixed Stromal",
    "Macrophages",
    "Endothelial",
    "T cells",
    "Muscle / Mural",
    "Ciliated Epithelial",
    "Unresolved / low-confidence"
  )
  
  category_colors <- c(
    "Tumor" = "#D73027",
    "Fibroblast" = "#1A9850",
    "Mixed Stromal" = "#91CF60",
    "Macrophages" = "#4575B4",
    "Endothelial" = "#FDAE61",
    "T cells" = "#984EA3",
    "Muscle / Mural" = "#A6761D",
    "Ciliated Epithelial" = "#66C2A5",
    "Unresolved / low-confidence" = "#BDBDBD"
  )
  
  observed_categories <- unique(coords_df$category)
  unexpected <- setdiff(observed_categories, expected_categories)
  
  if (length(unexpected) > 0) {
    stop(
      "Unexpected category labels: ",
      paste(unexpected, collapse = ", ")
    )
  }
  
  coords_df$category <- factor(
    coords_df$category,
    levels = expected_categories
  )
  
  log_message("Category counts:")
  print(table(coords_df$category, useNA = "ifany"))
  flush.console()
})

# ------------------------------------------------------------
# Step 7: Build spatial plot
# ------------------------------------------------------------

p_spatial <- run_step("Build spatial plot", {
  ggplot2::ggplot(
    coords_df,
    ggplot2::aes(x = x, y = y, color = category)
  ) +
    ggplot2::geom_point(
      size = 0.01,
      alpha = 0.85,
      shape = 16
    ) +
    ggplot2::scale_color_manual(
      values = category_colors,
      breaks = expected_categories,
      drop = TRUE
    ) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = "Spatial Distribution of Cell Categories",
      subtitle = "Human ovarian cancer Xenium dataset",
      color = "Cell category",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 15, hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10, hjust = 0.5
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
})

log_message(
  "Plot object created. PNG rendering/export begins next."
)

# ------------------------------------------------------------
# Step 8: Save PNG first
# ------------------------------------------------------------

run_step("Save 300-DPI PNG", {
  ggplot2::ggsave(
    filename = png_path,
    plot = p_spatial,
    width = 12,
    height = 9,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = TRUE
  )
  
  if (!file.exists(png_path)) {
    stop("PNG file was not created.")
  }
  
  log_message(
    paste("PNG file size:",
          round(file.info(png_path)$size / 1e6, 1),
          "MB")
  )
})

# ------------------------------------------------------------
# Step 9: Save a raster-based PDF
# ------------------------------------------------------------

if (requireNamespace("png", quietly = TRUE)) {
  run_step("Embed PNG in PDF", {
    image_array <- png::readPNG(png_path)
    
    grDevices::pdf(
      file = pdf_path,
      width = 12,
      height = 9,
      useDingbats = FALSE
    )
    
    tryCatch(
      {
        grid::grid.newpage()
        grid::grid.raster(
          image_array,
          width = grid::unit(1, "npc"),
          height = grid::unit(1, "npc"),
          interpolate = FALSE
        )
      },
      finally = grDevices::dev.off()
    )
    
    rm(image_array)
    
    if (!file.exists(pdf_path)) {
      stop("PDF file was not created.")
    }
    
    log_message(
      paste("PDF file size:",
            round(file.info(pdf_path)$size / 1e6, 1),
            "MB")
    )
  })
} else {
  log_message(
    paste(
      "WARNING: Package 'png' is not installed.",
      "PNG was saved successfully; PDF export was skipped",
      "to avoid a potentially expensive vector export."
    )
  )
}

# ------------------------------------------------------------
# Step 10: Finish and report outputs
# ------------------------------------------------------------

total_elapsed <- as.numeric(
  difftime(Sys.time(), start_time, units = "mins")
)

log_message(
  sprintf("RUN FINISHED | Total elapsed time: %.2f minutes",
          total_elapsed)
)

log_message(paste("PNG output:", png_path))

if (file.exists(pdf_path)) {
  log_message(paste("PDF output:", pdf_path))
}

log_message(paste("Progress log:", log_path))

cat("\nFigure 3 script finished.\n")
cat("PNG:", png_path, "\n")

if (file.exists(pdf_path)) {
  cat("PDF:", pdf_path, "\n")
}

cat("LOG:", log_path, "\n")