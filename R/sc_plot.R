#' Enhanced Visualization (DotPlot, Heatmap, VlnPlot, DimPlot, FeaturePlot) for Seurat
#'
#' @description A unified visualization function that supports five major plot types
#' with consistent aesthetic logic, professional corner-arrow axes, and clustering support.
#'
#' @param object A Seurat object.
#' @param features Vector or named list of features. Required for all types except 'dimplot'.
#' @param type Plot type: "featureplot", "dimplot", "dotplot", "heatmap", or "vlnplot".
#' @param group.by Metadata column to group cells by. Default is "celltype".
#' @param split.by Metadata column to split the plot by (faceting).
#' @param show_axis (Dim/Feature only) If TRUE, shows standard axes. If FALSE (default), shows clean corner arrows.
#' @param cols (DimPlot/VlnPlot) Vector of discrete colors. If NULL, uses a classic 20-color palette.
#' @param palette (Dot/Heat/Feature) Continuous palette name. Default "Spectral" for high contrast.
#' @param reduction (Dim/Feature only) Reduction to use (e.g., "umap", "harmony_umap"). Default "umap".
#' @param axis_title (Dim/Feature only) Title for axes labels. Default "UMAP".
#' @param pt.size Point size for DimPlot and FeaturePlot. Default 0.1.
#' @param raster (Dim/Feature only) If TRUE, rasterizes points using ggrastr to avoid large file sizes.
#' @param order (Dim/Feature only) If TRUE, brings high expression or specific groups to the front.
#' @param legend_ncol (DimPlot only) Number of columns for the legend.
#' @param arrow_len (Dim/Feature only) Length of corner arrows relative to plot range. Default 0.15.
#' @param cluster_rows (Heatmap only) If TRUE, clusters groups (Y-axis) based on average expression.
#' @param cluster_cols (Heatmap only) If TRUE, clusters features (X-axis).
#' @param dot.size.range (DotPlot only) Size range for dots. Default c(0, 6).
#' @param base_size Base font size for the entire plot. Default 7.
#' @param font_size_axis (Dim/Feature only) Font size for axis labels. Default 6.
#' @param font_size_legend Font size for legend text. Default 6.
#' @param x.text.angle Angle for x-axis labels (Dot/Vln/Heat). Default 90.
#' @param lw_param Line width for panel borders.
#' @param border_col Color for panel borders.
#' @param fill_col Background color for facet strips (labels).
#' @param color_name Legend title for color scale.
#' @param size_name (DotPlot only) Legend title for size scale.
#' @param bar.width Width of the color bar legend.
#' @param bar.height Height of the color bar legend.
#' @param rotate_vln (VlnPlot only) If TRUE, rotates violins to horizontal.
#'
#' @return A ggplot2 object.
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis magma plasma
#' @importFrom grDevices colorRampPalette
#' @importFrom stats dist hclust
#' @importFrom grid unit
#' @export
sc_plot <- function(
    object,
    features = NULL,
    type = c("featureplot", "dimplot", "dotplot", "heatmap", "vlnplot"),
    group.by = "celltype",
    split.by = NULL,
    show_axis = FALSE,
    cols = NULL,
    palette = "Spectral",
    reduction = "umap",
    axis_title = "UMAP",
    pt.size = 0.1,
    raster = FALSE,
    order = TRUE,
    legend_ncol = NULL,
    arrow_len = 0.15,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    dot.size.range = c(0, 6),
    base_size = 7,
    font_size_axis = 6,
    font_size_legend = 6,
    x.text.angle = 90,
    lw_param = 0.1,
    border_col = "black",
    fill_col = "grey92",
    color_name = "Expression",
    size_name = "% Exp",
    bar.width = 0.8,
    bar.height = 3,
    rotate_vln = FALSE
) {

  type <- match.arg(type)

  # =========================================================================
  # 1. Feature Validation & Gene Group Mapping
  # =========================================================================
  genes_plot <- NULL
  gene_group_map <- NULL

  if (type != "dimplot") {
    if (is.null(features)) stop("Error: 'features' argument is required for this plot type.")

    all_features <- if (is.list(features)) unique(unlist(features)) else features
    valid_features <- intersect(all_features, rownames(object))

    if (length(valid_features) == 0) stop("Error: None of the features provided were found in the object matrix.")

    if (is.list(features)) {
      features_list <- lapply(features, function(x) intersect(x, valid_features))
      features_list <- features_list[sapply(features_list, length) > 0]
      genes_plot <- unique(unlist(features_list))

      # Prepare the mapping table: Feature -> Gene Group
      gene_group_map <- stack(features_list)
      colnames(gene_group_map) <- c("features.plot", "gene_group")
      gene_group_map$features.plot <- as.character(gene_group_map$features.plot)
      gene_group_map$gene_group <- factor(gene_group_map$gene_group, levels = names(features_list))
    } else {
      genes_plot <- valid_features
    }
  }

  # =========================================================================
  # 2. Color Palette Setup
  # =========================================================================
  color_pal <- switch(palette,
                      "Spectral" = rev(RColorBrewer::brewer.pal(11, "Spectral")),
                      "RdYlBu"   = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                      "Viridis"  = viridis::viridis(100),
                      "Magma"    = viridis::magma(100),
                      "Plasma"   = viridis::plasma(100),
                      "Seurat"   = c("lightgrey", "blue"),
                      rev(RColorBrewer::brewer.pal(11, "Spectral"))
  )

  # =========================================================================
  # 3. Plotting Logic
  # =========================================================================

  p <- ggplot()

  # Standard Theme for Non-coordinate plots (Dot/Heat/Vln)
  base_theme <- theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = x.text.angle, hjust = 1, vjust = 0.5, face = "italic", color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title = element_blank(),
      panel.border = element_rect(color = border_col, fill = NA, linewidth = lw_param),
      strip.background = element_rect(color = border_col, fill = fill_col, linewidth = lw_param),
      strip.text = element_text(face = "bold", color = "black"),
      panel.spacing = unit(1, "pt"),
      panel.grid = element_blank()
    )

  # -------------------------------------------------------------------------
  # TYPE: FEATUREPLOT
  # -------------------------------------------------------------------------
  if (type == "featureplot") {

    # Raster Check
    if (raster && !requireNamespace("ggrastr", quietly = TRUE)) {
      warning("Package 'ggrastr' required for rasterization. Using vector points.")
      raster <- FALSE
    }

    if (!reduction %in% names(object@reductions)) stop(paste("Reduction", reduction, "not found."))

    # Data Prep
    embed <- Seurat::Embeddings(object, reduction)
    df_coords <- as.data.frame(embed[, 1:2])
    colnames(df_coords) <- c("DIM_1", "DIM_2")

    data_exp <- Seurat::FetchData(object, vars = genes_plot)
    df <- cbind(df_coords, data_exp) %>%
      tidyr::pivot_longer(cols = all_of(genes_plot), names_to = "features.plot", values_to = "expression")

    # Ordering
    if (isTRUE(order)) {
      df <- df %>% dplyr::arrange(expression)
    } else {
      set.seed(42); df <- df[sample(nrow(df)), ]
    }
    df$features.plot <- factor(df$features.plot, levels = genes_plot)

    # Plot Construction
    p <- ggplot(df, aes(x = DIM_1, y = DIM_2, color = expression))

    if (raster) {
      p <- p + ggrastr::rasterise(geom_point(size = pt.size, stroke = 0), dpi = 300)
    } else {
      p <- p + geom_point(size = pt.size, stroke = 0)
    }

    p <- p + scale_color_gradientn(colours = color_pal, name = color_name)

    # Axis Logic
    if (show_axis) {
      p <- p + theme_bw(base_size = base_size) +
        labs(x = paste0(axis_title, " 1"), y = paste0(axis_title, " 2")) +
        theme(aspect.ratio = 1, panel.grid = element_blank(),
              axis.text = element_text(color = "black", size = font_size_axis),
              axis.title = element_text(color = "black", size = font_size_axis),
              strip.text = element_text(face = "bold", size = base_size))
    } else {
      x_r <- range(df$DIM_1); y_r <- range(df$DIM_2)
      x_s <- x_r[1] - diff(x_r)*0.05; x_e <- x_s + diff(x_r)*arrow_len
      y_s <- y_r[1] - diff(y_r)*0.05; y_e <- y_s + diff(y_r)*arrow_len

      p <- p +
        annotate("segment", x=x_s, xend=x_e, y=y_s, yend=y_s, arrow=arrow(type="closed", length=unit(0.05,"in")), linewidth=0.3, color="black") +
        annotate("segment", x=x_s, xend=x_s, y=y_s, yend=y_e, arrow=arrow(type="closed", length=unit(0.05,"in")), linewidth=0.3, color="black") +
        annotate("text", x=x_s+diff(c(x_s,x_e))/2, y=y_s-diff(y_r)*0.02, label=paste(axis_title,"1"), size=font_size_axis/2.83, vjust=1) +
        annotate("text", x=x_s-diff(x_r)*0.02, y=y_s+diff(c(y_s,y_e))/2, label=paste(axis_title,"2"), size=font_size_axis/2.83, angle=90, vjust=0) +
        theme_void() + coord_fixed() +
        theme(strip.text = element_text(face="bold", size=base_size))
    }
    p <- p + facet_wrap(~features.plot) +
      guides(color = guide_colorbar(barwidth = bar.width, barheight = bar.height)) +
      theme(
        legend.title = element_text(size = base_size, face = "bold"),
        legend.text = element_text(size = base_size)
      )
  }

  # -------------------------------------------------------------------------
  # TYPE: DIMPLOT
  # -------------------------------------------------------------------------
  else if (type == "dimplot") {

    if (raster && !requireNamespace("ggrastr", quietly = TRUE)) raster <- FALSE
    if (!reduction %in% names(object@reductions)) stop(paste("Reduction", reduction, "not found."))

    embed <- Seurat::Embeddings(object, reduction)
    df <- as.data.frame(embed[, 1:2])
    colnames(df) <- c("DIM_1", "DIM_2")
    df$group <- as.factor(object@meta.data[[group.by]])

    if (isTRUE(order)) df <- df[order(df$group), ] else { set.seed(42); df <- df[sample(nrow(df)), ] }

    n_groups <- length(levels(df$group))
    if (is.null(cols)) {
      classic_cols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3")
      final_cols <- if(n_groups <= 20) classic_cols[1:n_groups] else grDevices::colorRampPalette(classic_cols)(n_groups)
    } else { final_cols <- cols }
    if (!is.null(legend_ncol)) {
      ncol_final <- legend_ncol
    } else {
      ncol_final <- ceiling(n_groups / 10) # Default to max 10 rows per column
    }
    nrow_final <- ceiling(n_groups / ncol_final)
    p <- ggplot(df, aes(x = DIM_1, y = DIM_2, color = group))
    if (raster) p <- p + ggrastr::rasterise(geom_point(size = pt.size, stroke = 0), dpi = 300) else p <- p + geom_point(size = pt.size, stroke = 0)

    p <- p + scale_color_manual(values = final_cols)

    if (show_axis) {
      p <- p + theme_bw(base_size = base_size) + labs(x = paste(axis_title,"1"), y = paste(axis_title,"2")) +
        theme(aspect.ratio = 1, panel.grid = element_blank())
    } else {
      x_r <- range(df$DIM_1); y_r <- range(df$DIM_2)
      x_s <- x_r[1]-diff(x_r)*0.05; x_e <- x_s+diff(x_r)*arrow_len
      y_s <- y_r[1]-diff(y_r)*0.05; y_e <- y_s+diff(y_r)*arrow_len

      p <- p +
        annotate("segment", x=x_s, xend=x_e, y=y_s, yend=y_s, arrow=arrow(type="closed", length=unit(0.05,"in")), linewidth=0.3, color="black") +
        annotate("segment", x=x_s, xend=x_s, y=y_s, yend=y_e, arrow=arrow(type="closed", length=unit(0.05,"in")), linewidth=0.3, color="black") +
        annotate("text", x=x_s+diff(c(x_s,x_e))/2, y=y_s-diff(y_r)*0.02, label=paste(axis_title,"1"), size=font_size_axis/2.83, vjust=1) +
        annotate("text", x=x_s-diff(x_r)*0.02, y=y_s+diff(c(y_s,y_e))/2, label=paste(axis_title,"2"), size=font_size_axis/2.83, angle=90, vjust=0) +
        theme_void() + coord_fixed()
    }

    p <- p + theme(legend.text = element_text(size = font_size_legend), legend.key.size = unit(0.12, "in")) +
      guides(color = guide_legend(title = NULL, nrow = nrow_final,
                                  byrow = FALSE,
                                  override.aes = list(size = 2.5, alpha = 1)))
  }

  # -------------------------------------------------------------------------
  # TYPE: DOTPLOT
  # -------------------------------------------------------------------------
  else if (type == "dotplot") {

    run_group <- group.by
    lookup <- NULL

    # 1. Handle Split Logic
    if (!is.null(split.by)) {
      # Create safe merged ID
      object$temp_id <- paste(object@meta.data[[group.by]], object@meta.data[[split.by]], sep = "_SEP_")
      run_group <- "temp_id"
      lookup <- data.frame(id = object$temp_id, g = object@meta.data[[group.by]], s = object@meta.data[[split.by]]) %>% unique()
    }

    # 2. Run Seurat DotPlot
    p_base <- Seurat::DotPlot(object, features = genes_plot, group.by = run_group)
    df_plot <- p_base$data

    # 3. Restore Metadata
    if (!is.null(split.by)) {
      df_plot <- dplyr::left_join(df_plot, lookup, by = "id")
      # Fix Factor Levels (Critical for Order)
      levs <- levels(object@meta.data[[group.by]])
      if(is.null(levs)) levs <- unique(as.character(object@meta.data[[group.by]]))
      df_plot$id_for_plot <- factor(df_plot$g, levels = rev(levs))
      df_plot$split_for_plot <- df_plot$s
    } else {
      df_plot$id_for_plot <- df_plot$id
    }

    # 4. Join Gene Group Map (The Fix for 'missing gene_group' error)
    if (is.list(features) && !is.null(gene_group_map)) {
      df_plot <- dplyr::left_join(df_plot, gene_group_map, by = "features.plot")
    }

    # 5. Build Plot
    p <- ggplot(df_plot, aes(x = features.plot, y = id_for_plot)) +
      geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
      scale_color_gradientn(colours = color_pal, name = color_name) +
      scale_size(range = dot.size.range, name = size_name) +
      guides(color = guide_colorbar(barwidth = bar.width, barheight = bar.height, order = 1),
             size = guide_legend(title = size_name, order = 2)) +
      base_theme

    # 6. Apply Faceting
    if (!is.null(split.by) && is.list(features)) {
      p <- p + facet_grid(. ~ split_for_plot + gene_group, scales = "free_x", space = "free_x")
    } else if (!is.null(split.by)) {
      p <- p + facet_grid(. ~ split_for_plot, scales = "fixed", space = "fixed")
    } else if (is.list(features)) {
      p <- p + facet_grid(. ~ gene_group, scales = "free_x", space = "free_x")
    }
  }

  # -------------------------------------------------------------------------
  # TYPE: HEATMAP
  # -------------------------------------------------------------------------
  else if (type == "heatmap") {

    groups_to_avg <- group.by
    if (!is.null(split.by)) groups_to_avg <- c(group.by, split.by)

    # 1. Calc Avg Exp
    avg_exp <- Seurat::AverageExpression(object, features = genes_plot, group.by = groups_to_avg, slots = "data")[[1]]
    mat_scaled <- t(scale(t(avg_exp)))

    # 2. Clustering
    if (cluster_cols) genes_plot <- rownames(mat_scaled)[stats::hclust(stats::dist(mat_scaled))$order]
    if (cluster_rows) {
      groups_order <- colnames(mat_scaled)[stats::hclust(stats::dist(t(mat_scaled)))$order]
    }

    # 3. Reshape
    df_plot <- as.data.frame(mat_scaled) %>%
      dplyr::mutate(features.plot = rownames(.)) %>%
      tidyr::pivot_longer(cols = -features.plot, names_to = "full", values_to = "avg.exp.scaled")

    # 4. Handle Split
    if (!is.null(split.by)) {
      # Robust splitting of the combined group name
      meta_comb <- unique(object@meta.data[, c(group.by, split.by)])
      colnames(meta_comb) <- c("real_group", "real_split")
      # Simulate Seurat's Paste
      meta_comb$full_sim <- paste(meta_comb$real_group, meta_comb$real_split, sep = "_")

      df_plot <- dplyr::left_join(df_plot, meta_comb, by = c("full" = "full_sim"))

      # Fallback if join failed
      if(any(is.na(df_plot$real_group))) {
        df_plot <- df_plot %>% tidyr::separate(full, into = c("real_group", "real_split"), sep = "_", extra = "merge", fill = "right")
      }
      df_plot$id_for_plot <- df_plot$real_group
      df_plot$split_for_plot <- df_plot$real_split
    } else {
      df_plot$id_for_plot <- df_plot$full
      if(cluster_rows) df_plot$id_for_plot <- factor(df_plot$id_for_plot, levels = groups_order)
    }

    df_plot$features.plot <- factor(df_plot$features.plot, levels = genes_plot)

    # 5. Join Gene Group Map (The Fix)
    if (is.list(features) && !is.null(gene_group_map)) {
      df_plot <- dplyr::left_join(df_plot, gene_group_map, by = "features.plot")
    }

    # 6. Plot
    p <- ggplot(df_plot, aes(x = features.plot, y = id_for_plot, fill = avg.exp.scaled)) +
      geom_tile() +
      scale_fill_gradientn(colours = color_pal, name = color_name) +
      guides(fill = guide_colorbar(barwidth = bar.width, barheight = bar.height)) +
      base_theme + theme(panel.grid = element_blank())

    # 7. Apply Faceting
    # Note: If clustering is ON, we generally don't facet by gene group as it breaks the tree visual
    use_gene_facets <- is.list(features) && !cluster_cols

    if (!is.null(split.by) && use_gene_facets) {
      p <- p + facet_grid(. ~ split_for_plot + gene_group, scales = "free_x", space = "free_x")
    } else if (!is.null(split.by)) {
      p <- p + facet_grid(. ~ split_for_plot, scales = "fixed", space = "fixed")
    } else if (use_gene_facets) {
      p <- p + facet_grid(. ~ gene_group, scales = "free_x", space = "free_x")
    }
  }

  # -------------------------------------------------------------------------
  # TYPE: VLNPLOT
  # -------------------------------------------------------------------------
  else if (type == "vlnplot") {

    cols_to_fetch <- c(genes_plot, group.by)
    if (!is.null(split.by)) cols_to_fetch <- c(cols_to_fetch, split.by)

    df_plot <- Seurat::FetchData(object, vars = cols_to_fetch) %>%
      tidyr::pivot_longer(cols = all_of(genes_plot), names_to = "features.plot", values_to = "expression")

    df_plot$id_for_plot <- df_plot[[group.by]]
    if (!is.null(split.by)) df_plot$split_for_plot <- df_plot[[split.by]]

    # Color Logic
    n_groups <- if(!is.null(split.by)) length(unique(df_plot[[split.by]])) else length(unique(df_plot$id_for_plot))

    if (is.null(cols)) {
      classic_cols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
      vln_cols <- if(n_groups <= 10) classic_cols[1:n_groups] else grDevices::colorRampPalette(classic_cols)(n_groups)
    } else {
      vln_cols <- cols
    }

    df_plot$features.plot <- factor(df_plot$features.plot, levels = genes_plot)

    # Join Gene Map (For robustness, though Vln usually facets by Feature anyway)
    if (is.list(features) && !is.null(gene_group_map)) {
      df_plot <- dplyr::left_join(df_plot, gene_group_map, by = "features.plot")
    }

    p <- ggplot(df_plot, aes(x = id_for_plot, y = expression, fill = if(!is.null(split.by)) .data[[split.by]] else id_for_plot)) +
      geom_violin(scale = "width", trim = TRUE) +
      scale_fill_manual(values = vln_cols) +
      base_theme +
      theme(legend.position = "none") + ylab("Exp")

    if (rotate_vln) p <- p + coord_flip()

    # Faceting
    p <- p + facet_wrap(~features.plot, scales = "free_y")
  }

  return(p)
}
