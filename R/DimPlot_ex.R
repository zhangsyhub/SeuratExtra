#' @title Journal-Quality Dimension Reduction Plot
#' @description A customized DimPlot with corner-arrow axes and optimized legends, suitable for high-impact publications.
#' @param object A Seurat object.
#' @param reduction Character. The reduction to use (e.g., "umap", "tsne"). Default is "umap".
#' @param group.by Character. The meta.data column to group cells by. Default is "celltype".
#' @param axis_title Character. Prefix for axis labels. Default is "UMAP".
#' @param pt.size Numeric. Size of the points. Default is 0.1.
#' @param raster Logical. Whether to rasterize the points (requires ggrastr). Recommended for >50k cells. Default is FALSE.
#' @param user_colors Vector. Custom color palette. If NULL, an optimized high-contrast palette is used.
#' @param legend_ncol Integer. Number of columns for the legend. If NULL, calculated automatically for balance.
#' @param font_size_base Numeric. Base font size for the plot. Default is 7.
#' @param font_size_axis Numeric. Font size for axis labels. Default is 6.
#' @param font_size_legend Numeric. Font size for legend text. Default is 6.
#' @param arrow_len Numeric. Relative length of the axis arrows. Default is 0.15.
#' @param order Logical. If TRUE, points are ordered by group; if FALSE (default), points are shuffled to prevent density bias.
#'
#' @import ggplot2
#' @import Seurat
#' @importFrom grid unit
#' @importFrom grDevices colorRampPalette
#' @importFrom ggrastr rasterise
#'
#' @return A ggplot object.
#' @export
DimPlot_ex <- function(object,
                       reduction = "umap",
                       group.by = "celltype",
                       axis_title = "UMAP",
                       pt.size = 0.1,
                       raster = FALSE,
                       user_colors = NULL,
                       legend_ncol = NULL,
                       font_size_base = 7,
                       font_size_axis = 6,
                       font_size_legend = 6,
                       arrow_len = 0.15,
                       order = FALSE) {

  # 1. Dependency & Input Check
  if (raster && !requireNamespace("ggrastr", quietly = TRUE)) {
    warning("Package 'ggrastr' is required for rasterization. Falling back to vector plot.")
    raster <- FALSE
  }

  if (!reduction %in% names(object@reductions)) {
    stop(paste("Reduction", reduction, "not found in object."))
  }

  if (!group.by %in% colnames(object@meta.data)) {
    stop(paste("Group", group.by, "not found in meta.data."))
  }

  # 2. Extract and Prepare Data
  embed <- Seurat::Embeddings(object, reduction)
  df <- as.data.frame(embed[, 1:2])
  colnames(df) <- c("DIM_1", "DIM_2")
  df$group <- as.factor(object@meta.data[[group.by]])

  # Shuffling/Ordering Logic
  if (isTRUE(order)) {
    df <- df[order(df$group), ]
  } else {
    set.seed(42) # Ensure reproducibility
    df <- df[sample(nrow(df)), ]
  }

  # 3. Palette Logic
  n_groups <- length(levels(df$group))
  if (is.null(user_colors)) {
    classic_cols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                      "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
    if(n_groups <= length(classic_cols)) {
      cols <- classic_cols[1:n_groups]
    } else {
      cols <- grDevices::colorRampPalette(classic_cols)(n_groups)
    }
  } else {
    cols <- user_colors
    if(length(cols) < n_groups) warning("Provided colors are fewer than groups!")
  }

  # 4. Legend Layout Logic
  if (!is.null(legend_ncol)) {
    ncol_final <- legend_ncol
  } else {
    ncol_final <- ceiling(n_groups / 10) # Default to max 10 rows per column
  }
  nrow_final <- ceiling(n_groups / ncol_final)

  # 5. Coordinate & Arrow Calculations
  x_range <- range(df$DIM_1)
  y_range <- range(df$DIM_2)

  # Add a small buffer (5%) to prevent arrows from overlapping with edge points
  x_start <- x_range[1] - diff(x_range) * 0.05
  y_start <- y_range[1] - diff(y_range) * 0.05

  x_end <- x_start + diff(x_range) * arrow_len
  y_end <- y_start + diff(y_range) * arrow_len

  text_offset_x <- diff(x_range) * 0.02
  text_offset_y <- diff(y_range) * 0.02

  # 6. Plotting
  p <- ggplot2::ggplot(df, ggplot2::aes(x = DIM_1, y = DIM_2, color = group))

  if (raster) {
    p <- p + ggrastr::rasterise(ggplot2::geom_point(size = pt.size, stroke = 0), dpi = 600)
  } else {
    p <- p + ggplot2::geom_point(size = pt.size, stroke = 0)
  }

  p <- p +
    ggplot2::scale_color_manual(values = cols) +
    # Corner Arrows
    ggplot2::annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_start,
                      arrow = ggplot2::arrow(type = "closed", length = grid::unit(0.05, "inches"), angle = 20),
                      linewidth = 0.3, color = "black") +
    ggplot2::annotate("segment", x = x_start, xend = x_start, y = y_start, yend = y_end,
                      arrow = ggplot2::arrow(type = "closed", length = grid::unit(0.05, "inches"), angle = 20),
                      linewidth = 0.3, color = "black") +
    # Axis Text Labels
    ggplot2::annotate("text", x = x_start + (x_end - x_start)/2, y = y_start - text_offset_y,
                      label = paste0(axis_title, " 1"), size = font_size_axis / 2.83,
                      vjust = 1, hjust = 0.5) +
    ggplot2::annotate("text", x = x_start - text_offset_x, y = y_start + (y_end - y_start)/2,
                      label = paste0(axis_title, " 2"), size = font_size_axis / 2.83,
                      angle = 90, vjust = 0, hjust = 0.5) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      text = ggplot2::element_text(size = font_size_base, family = "sans"),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 20, l = 20, unit = "pt"),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = font_size_legend),
      legend.key.size = grid::unit(0.12, "inches"),
      legend.spacing.x = grid::unit(0.05, "inches")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      nrow = nrow_final,
      byrow = FALSE,
      override.aes = list(size = 2.5, alpha = 1)
    ))

  return(p)
}
