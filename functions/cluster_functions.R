
plotAbundancesCustom <-
  function (x,
            k,
            by = c("sample_id", "cluster_id"),
            group_by = "condition",
            shape_by = NULL,
            col_clust = TRUE,
            distance = c("euclidean",
                         "maximum",
                         "manhattan",
                         "canberra",
                         "binary",
                         "minkowski"),
            linkage = c(
              "average",
              "ward.D",
              "single",
              "complete",
              "mcquitty",
              "median",
              "centroid",
              "ward.D2"
            ),
            k_pal = CATALYST:::.cluster_cols)
  {
    library(ggplot2)
    by <- match.arg(by)
    CATALYST:::.check_sce(x, TRUE)
    k <- CATALYST:::.check_k(x, k)
    CATALYST:::.check_cd_factor(x, group_by)
    CATALYST:::.check_cd_factor(x, shape_by)
    CATALYST:::.check_pal(k_pal)
    linkage <- match.arg(linkage)
    distance <- match.arg(distance)
    stopifnot(is.logical(col_clust), length(col_clust) == 1)
    shapes <- CATALYST:::.get_shapes(x, shape_by)
    if (is.null(shapes))
      shape_by <- NULL
    if (by == "sample_id") {
      nk <- nlevels(CATALYST::cluster_ids(x, k))
      if (length(k_pal) < nk)
        k_pal <- colorRampPalette(k_pal)(nk)
    }
    ns <-
      table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    m <- match(df$sample_id, x$sample_id)
    for (i in c(shape_by, group_by))
      df[[i]] <- x[[i]][m]
    if (by == "sample_id" &&
        col_clust && length(unique(df$sample_id)) >
        1) {
      d <- dist(t(fq), distance)
      h <- hclust(d, linkage)
      o <- colnames(fq)[h$order]
      df$sample_id <- factor(df$sample_id, o)
    }
    p <- ggplot2::ggplot(df, ggplot2::aes_string(y = "Freq")) + ggplot2::labs(x = NULL,
                                                   y = "Proportion [%]") + ggplot2::theme_bw() + ggplot2::theme(
                                                     panel.grid = ggplot2::element_blank(),
                                                     strip.text = ggplot2::element_text(face = "bold"),
                                                     strip.background = ggplot2::element_rect(fill = NA,
                                                                                     color = NA),
                                                     axis.text = ggplot2::element_text(color = "black"),
                                                     axis.text.x = ggplot2::element_text(
                                                       angle = 45,
                                                       hjust = 1,
                                                       vjust = 1
                                                     ),
                                                     legend.key.height = grid::unit(0.8, "lines")
                                                   )
    switch(
      by,
      sample_id = p + (if (!is.null(group_by))
        ggplot2::facet_wrap(group_by,
                   scales = "free_x")) + ggplot2::geom_bar(
                     ggplot2::aes_string(x = "sample_id",
                                fill = "cluster_id"),
                     position = "fill",
                     stat = "identity"
                   ) +
        ggplot2::scale_fill_manual("cluster_id", values = k_pal) + ggplot2::scale_x_discrete(expand = c(0,
                                                                                      0)) + ggplot2::scale_y_continuous(expand = c(0, 0), labels = seq(0,
                                                                                                                                              100, 25)) + ggplot2::theme(
                                                                                                                                                panel.border = ggplot2::element_blank(),
                                                                                                                                                panel.spacing.x = grid::unit(1,
                                                                                                                                                                       "lines")
                                                                                                                                              ),
      cluster_id = {
        p <-
          p + ggplot2::scale_shape_manual(values = shapes) + ggplot2::guides(
            col = ggplot2::guide_legend(order = 1,
                               override.aes = list(size = 3)),
            shape = ggplot2::guide_legend(override.aes = list(size = 3))
          )
        if (is.null(group_by)) {
          p + ggplot2::geom_boxplot(
            ggplot2::aes_string(x = "cluster_id"),
            alpha = 0.2,
            position = ggplot2::position_dodge(),
            outlier.color = NA
          ) +
            ggplot2::geom_point(ggplot2::aes_string("cluster_id", shape = shape_by),
                       position = ggplot2::position_jitter(width = 0.2))
        } else {
          p + ggplot2::facet_wrap("cluster_id", scales = "free_y",
                         ncol = 4) + ggplot2::geom_boxplot(
                           ggplot2::aes_string(
                             x = group_by,
                             color = group_by,
                             fill = group_by
                           ),
                           position = ggplot2::position_dodge(),
                           alpha = 0.2,
                           outlier.color = NA,
                           show.legend = FALSE
                         ) +
            ggplot2::geom_point(ggplot2::aes_string(
              x = group_by,
              col = group_by,
              shape = shape_by
            ),
            position = ggplot2::position_jitter(width = 0.2))
        }
      }
    )
  }

calcClusterFreqBySample <- function (x, k){
    CATALYST:::.check_sce(x, TRUE)
    k <- CATALYST:::.check_k(x, k)
    ns <-
      table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
}
    
plotClusterExprsCustom <-
  function (x,
            k,
            features = SummarizedExperiment::rowData(x)$marker_name[SummarizedExperiment::rowData(x)$used_for_clustering],
            assay = "exprs")
  {
    CATALYST:::.check_sce(x, TRUE)
    k <- CATALYST:::.check_k(x, k)
    x$cluster_id <- CATALYST::cluster_ids(x, k)
    features <- CATALYST:::.get_features(x, features)
    ms <-
      t(CATALYST:::.agg(x[features,], "cluster_id", "median", assay = assay))
    d <- dist(ms, method = "euclidean")
    o <- hclust(d, method = "average")$order
    cd <- SummarizedExperiment::colData(x)
    es <- assay(x[features,], assay)
    df <-
      data.table::as.data.table(data.frame(t(es), cd, check.names = FALSE))
    df <-
      data.table::melt(
        df,
        id.vars = names(cd),
        variable.name = "antigen",
        value.name = "expression"
      )
    df$avg <- "no"
    avg <- df
    avg$cluster_id <- "avg"
    avg$avg <- "yes"
    df <- rbind(df, avg)
    fq <- tabulate(x$cluster_id) / ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    df$cluster_id <- factor(df$cluster_id,
                            levels = rev(c("avg",
                                           levels(x$cluster_id)[o])),
                            labels = rev(c(
                              "average",
                              paste0(names(fq), " (", fq, "%)")[o]
                            )))
    ggplot2::ggplot(df,
                    ggplot2::aes_string(
             x = "expression",
             y = "cluster_id",
             col = "avg",
             fill = "avg"
           )) + ggplot2::facet_wrap( ~ antigen, scales = "free_x",
                            nrow = 2) + ggridges::geom_density_ridges(alpha = 0.2) + ggridges::theme_ridges() +
      ggplot2::theme(
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold")
      )
  }

plotFreqHeatmapCustom <- function (x,
                                   k,
                                   m = NULL,
                                   normalize = TRUE,
                                   row_anno = TRUE,
                                   col_anno = TRUE,
                                   row_clust = TRUE,
                                   col_clust = TRUE,
                                   row_dend = TRUE,
                                   col_dend = TRUE,
                                   bars = TRUE,
                                   perc = FALSE,
                                   hm_pal = rev(RColorBrewer::brewer.pal(11,
                                                                         "RdBu")),
                                   k_pal = CATALYST:::.cluster_cols,
                                   m_pal = k_pal)
{
  library(ComplexHeatmap)
  args <- as.list(environment())
  CATALYST:::.check_args_plotFreqHeatmap(args)
  x$cluster_id <- CATALYST::cluster_ids(x, k)
  ns <- table(x$cluster_id, x$sample_id)
  fq <- prop.table(ns, 2)
  y <- as.matrix(unclass(fq))
  if (normalize)
    y <- CATALYST:::.z_normalize(asin(sqrt(y)))
  if (!isFALSE(row_anno)) {
    left_anno <- CATALYST:::.anno_clusters(x, k, m, k_pal, m_pal)
  }
  else
    left_anno <- NULL
  if (!isFALSE(col_anno)) {
    top_anno <-
      CATALYST:::.anno_factors(x, levels(x$sample_id), col_anno,
                               "colum")
  }
  else
    top_anno <- NULL
  if (bars) {
    right_anno <- CATALYST:::.anno_counts(x$cluster_id, perc)
  }
  else
    right_anno <- NULL
  ComplexHeatmap::Heatmap(
    matrix = y,
    name = paste0("normalized\n"[normalize],
                  "frequency"),
    col = hm_pal,
    na_col = "lightgrey",
    rect_gp = gpar(col = "white"),
    column_title = "sample_id",
    column_title_side = "bottom",
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    show_row_dend = row_dend,
    show_column_dend = col_dend,
    show_row_names = is.null(left_anno),
    row_names_side = "left",
    top_annotation = top_anno,
    left_annotation = left_anno,
    right_annotation = right_anno
  )
}

library(SingleCellExperiment)
library(FlowSOM)
library(S4Vectors)
library(ggplot2)
library(ggpubr)
library(grid) 
library(gridExtra) 
library(cowplot) 
library(patchwork)

# SOM added to metadata
clusterSCE <-
  function (x,
            assayType = "exprs",
            features = "type",
            xdim = 10,
            ydim = 10,
            maxK = 20,
            verbose = TRUE,
            seed = 1){
    # set seed
    if (!is.null(seed)) {
      set.seed(seed)
    }
    # check input
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(
      is.logical(verbose),
      length(verbose) == 1,
      vapply(list(xdim,
                  ydim, maxK, seed), function(arg)
                    is.numeric(arg) &&
               length(arg) == 1, logical(1))
    )
    # extract markers to cluster
    features <- CATALYST:::.get_features(x, features)
    if (is.null(CATALYST::marker_classes(x))) {
      SummarizedExperiment::rowData(x)$marker_class <-
        factor(c("state", "type")[as.numeric(rownames(x) %in%
                                               features) + 1], levels = c("type", "state", "none"))
    }
    SummarizedExperiment::rowData(x)$used_for_clustering <- rownames(x) %in% features
    if (verbose)
      message("o running FlowSOM clustering...")
    fsom <-
      FlowSOM::ReadInput(flowCore::flowFrame(t(SummarizedExperiment::assay(x, assayType))))
    som <-
      FlowSOM::BuildSOM(
        fsom,
        colsToUse = features,
        silent = FALSE,
        xdim = xdim,
        ydim = ydim
      )
    som <- FlowSOM::BuildMST(som)
    if (verbose)
      message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <-
      suppressWarnings(
        ConsensusClusterPlus::ConsensusClusterPlus(
          t(som$map$codes),
          maxK = maxK,
          reps = 100,
          distance = "euclidean",
          seed = seed,
          plot = NULL
        )
      )
    dev.off()
    # compute metacluster ecdf values
    Kvec <- 2:maxK
    names(Kvec) <- Kvec
    mc_dt <- rbindlist(lapply(Kvec, function(i){
      consensus_matrix <- mc[[i]]$consensusMatrix
      consensus_values <- consensus_matrix[lower.tri(consensus_matrix)]
      ecdf_data <- ecdf(consensus_values)
      list(ConsensusIndex = consensus_values, CDF = ecdf_data(consensus_values))
    }), idcol = 'k')
    mc_dt[, k:=factor(k, levels=Kvec)]
    
    k <- xdim * ydim
    mcs <- seq_len(maxK)[-1]
    codes <-
      data.frame(seq_len(k), purrr::map(mc[-1], "consensusClass"))
    codes <-
      dplyr::mutate_all(codes, function(u)
        factor(u, levels = sort(unique(u))))
    colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s",
                                                      mcs))
    x$cluster_id <- factor(som$map$mapping[, 1])
    S4Vectors::metadata(x)$cluster_codes <- codes
    S4Vectors::metadata(x)$SOM <- som
    # S4Vectors::metadata(x)$delta_area <- CATALYST:::.plot_delta_area(mc)
    S4Vectors::metadata(x)$mc_dt <- mc_dt
    x <-
      CATALYST::mergeClusters(
        x,
        k = sprintf("meta%s", maxK),
        id = "all",
        table = data.frame(old_cluster = seq_len(maxK), new_cluster = "all")
      )
    return(x)
  }


# plot ecdf
plot_ecdf <- function(sce, interactive = TRUE) {
  require(ggplot2)
  mc_dt <- metadata(sce)$mc_dt
  stopifnot(!is.null(mc_dt))
  ggp <-
    ggplot(mc_dt, aes(x = ConsensusIndex, y = CDF, color = k)) + geom_step() + theme_bw()
  if (interactive)
    ggp <- plotly::ggplotly(ggp)
  return(ggp)
}

add_max_curvature <- function(dt, column, h = 1) {
  # Ensure dt is a data.table
  if (!inherits(dt, "data.table")) {
    stop("dt must be a data.table")
  }
  
  # Calculate lagged and lead values of the target column
  dt[, paste0(column, "_lagged") := data.table::shift(get(column), type = "lag")]
  dt[, paste0(column, "_lead") := data.table::shift(get(column), type = "lead")]
  
  # Calculate the second derivative
  dt[, second_derivative := (get(paste0(column, "_lead")) - 2 * get(column) + get(paste0(column, "_lagged"))) / (h ^
                                                                                                                   2)]
  
  # Identify the index of the maximum curvature point
  max_curv_index <- which.max(abs(dt$second_derivative))
  
  # Add a column to indicate the maximum curvature point
  dt[, `maximum curvature` := 'no']
  dt[max_curv_index, `maximum curvature` := 'maximum curvature']
  
  # Optionally, clean up by removing the intermediate columns
  dt[, c(paste0(column, "_lagged"),
         paste0(column, "_lead"),
         "second_derivative") := NULL]
  
  return(invisible(dt))
}

# plot delta area
plot_delta_area <- function(sce, interactive = TRUE) {
  mc_dt <- metadata(sce)$mc_dt
  stopifnot(!is.null(mc_dt))
  
  mc_dt <- mc_dt[order(k, ConsensusIndex)]
  
  area_under_curve <- mc_dt[, .(AUC = {
    # Calculate the widths (differences in ConsensusIndex) and heights (CDF values)
    widths <- diff(.SD$ConsensusIndex)
    heights <- head(.SD$CDF,-1)
    # Calculate the area of each rectangle and sum them up
    sum(widths * heights)
  }), by = k]
  
  # Compute delta Area
  area_under_curve[, DeltaAUC := c(NA, diff(AUC))]
  area_under_curve[k == '2', DeltaAUC := AUC]
  
  add_max_curvature(area_under_curve, "DeltaAUC")
  
  # View the results
  ggp <-
    ggplot(area_under_curve, aes(x = k, y = DeltaAUC, group = 1)) + geom_line() + geom_point(aes(color =
                                                                                                   `maximum curvature`)) +
    scale_color_manual(values = c('maximum curvature' = 'red')) + theme_bw()
  if (interactive) {
    ggp <- plotly::ggplotly(ggp)
    for (i in seq_along(ggp$x$data))
      if (length(ggp$x$data[[i]]$legendgroup) > 0)
        if (ggp$x$data[[i]]$legendgroup == "no")
          ggp$x$data[[i]]$showlegend <- FALSE
  }
  
  return(ggp)
}

# plot PAC
plot_pac <- function(sce,
                     interactive = TRUE,
                     x1 = .05,
                     x2 = 1 - x1) {
  mc_dt <- metadata(sce)$mc_dt
  stopifnot(!is.null(mc_dt))
  
  pac_dt <-
    merge(mc_dt[ConsensusIndex <= x1, .(y1 = max(CDF)), by = k], mc_dt[ConsensusIndex <= x2, .(y2 =
                                                                                                 max(CDF)), by = k])
  pac_dt[, PAC := y2 - y1]
  pac_dt[, `minimum PAC` := 'no']
  pac_dt[PAC == min(PAC), `minimum PAC` := 'minimum PAC']
  
  add_max_curvature(pac_dt, "PAC")
  
  # View the results
  ggp <-
    ggplot(pac_dt, aes(x = k,
                       y = PAC,
                       group = 1)) + geom_line() + geom_point(aes(color = `maximum curvature`, shape = `minimum PAC`)) +
    scale_color_manual(values = c('maximum curvature' = 'red')) + scale_shape_manual(values = c('minimum PAC' = 17, 'no' = 16), breaks= 'minimum PAC') + theme_bw()
  if (interactive) {
    ggp <- plotly::ggplotly(ggp)
    for (i in seq_along(ggp$x$data))
      if (length(ggp$x$data[[i]]$legendgroup) > 0)
        if (ggp$x$data[[i]]$legendgroup == "no")
          ggp$x$data[[i]]$showlegend <- FALSE
  }
  return(ggp)
}

# function for plotting star chart (FlowSOM)
plotStarsCustom <- function (fsom, markers = fsom$map$colsUsed, overall = TRUE, nodeValues = NULL, nodeColors = NULL, colorPalette = FlowSOM_colors, 
                             list_insteadof_ggarrange = FALSE, ...) 
{
  fsom <- FlowSOM::UpdateFlowSOM(fsom)
  if (length(names(list(...))) > 0 && "backgroundColor" %in% 
      names(list(...))) {
    warning(paste0("\"backgroundColor\" is deprecated, ", 
                   "please use \"backgroundColors\" instead."))
  }
  channels <- GetChannels(fsom, markers)
  p <- PlotFlowSOMCustom(fsom = fsom, nodeValues = nodeValues, nodeColors = nodeColors, ...)
  if (!is.null(names(colorPalette))) {
    names(colorPalette) <- GetChannels(fsom, names(colorPalette))
  }
  if(overall){
    p <- AddStars(p = p, fsom = fsom, markers = channels, colorPalette = colorPalette)
    if (!is.null(names(colorPalette))) {
      names(colorPalette) <- fsom$prettyColnames[GetChannels(fsom, 
                                                             names(colorPalette))]
    }
    l1 <- PlotStarLegend(channels, colorPalette)
    l2 <- ggpubr::get_legend(p, position = "bottom")
    if (list_insteadof_ggarrange) {
      p <- p + ggplot2::theme(legend.position = "none")
      l2 <- ggpubr::as_ggplot(l2)
      return(list(tree = p, starLegend = l1, backgroundLegend = l2))
    }
    else {
      p <- ggpubr::ggarrange(p, ggpubr::ggarrange(l1, l2, 
                                                  ncol = 1), NULL, ncol = 3, widths = c(3, 1, 0.3), 
                             legend = "none")
    }
  }
  return(p)
}

# function for plotting with FlowSOM (FlowSOM)
PlotFlowSOMCustom <- function (fsom, nodeValues = NULL, nodeColors = NULL, view = "MST", nodeSizes = fsom$map$pctgs, maxNodeSize = 1, 
                               refNodeSize = max(nodeSizes), equalNodeSize = FALSE, backgroundValues = NULL, 
                               backgroundColors = NULL, backgroundLim = NULL, title = NULL) 
{
  requireNamespace("ggplot2")
  fsom <- FlowSOM::UpdateFlowSOM(fsom)
  nNodes <- NClusters(fsom)
  isEmpty <- fsom$map$pctgs == 0
  if (length(nodeSizes) != nNodes) {
    stop(paste0("Length of 'nodeSizes' should be equal to number of clusters ", 
                "in FlowSOM object (", nNodes, " clusters and ", 
                length(nodeSizes), " node sizes)."))
  }
  if (length(backgroundValues) != nNodes && !is.null(backgroundValues)) {
    stop(paste0("Length of 'backgroundValues' should be equal to number of ", 
                "clusters in FlowSOM object (", nNodes, " clusters and ", 
                length(backgroundValues), " backgroundValues)."))
  }
  if (deparse(substitute(backgroundColors)) == "backgroundColor") {
    warning(paste0("In the new FlowSOM version \"backgroundColors\"", 
                   " is used instead of \"backgroundColor\""))
  }
  layout <- FlowSOM:::ParseLayout(fsom, view)
  if (is.matrix(view) || is.data.frame(view)) 
    view <- "matrix"
  autoNodeSize <- FlowSOM:::AutoMaxNodeSize(layout = layout, overlap = ifelse(view %in% 
                                                                                c("grid"), -0.3, 1))
  maxNodeSize <- autoNodeSize * maxNodeSize
  if (equalNodeSize) {
    scaledNodeSize <- rep(maxNodeSize, nNodes)
  }
  else {
    scaledNodeSize <- FlowSOM:::ParseNodeSize(nodeSizes, maxNodeSize, 
                                              refNodeSize)
  }
  if (any(isEmpty)) {
    scaledNodeSize[isEmpty] <- min(maxNodeSize, 0.05)
  }
  plot_df <- data.frame(x = layout$x, y = layout$y, size = scaledNodeSize, 
                        bg_size = scaledNodeSize * 1.5)
  p <- ggplot2::ggplot(plot_df)
  if (!is.null(backgroundValues)) {
    p <- FlowSOM:::AddBackground(p, backgroundValues = backgroundValues, 
                                 backgroundColors = backgroundColors, backgroundLim = backgroundLim)
  }
  if (view == "MST") {
    p <- FlowSOM:::AddMST(p, fsom)
  }
  if(is.null(nodeValues)){
    p <- FlowSOM:::AddNodes(p = p, values = as.character(isEmpty), colorPalette = c(`TRUE` = "gray", 
                                                                                    `FALSE` = "white"), showLegend = FALSE)
  } else {
    p <- FlowSOM:::AddNodes(p = p, values = nodeValues, colorPalette = nodeColors, showLegend = FALSE)
  }
  p <- p + guides(fill_new = guide_legend(title = "Clusters")) # diese Zeile muss auch in die original PlotStars Funktion von FlowSOM
  p <- p + ggplot2::coord_fixed() + ggplot2::theme_void()
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  } 
  return(p)
}


# own function for enabling selection of groups and facetting
plotMarkerCustom <- function (sce, marker, facet_by = "", subselection_col = "", subselection=NULL, assayType = "exprs", colorPalette = grDevices::colorRampPalette(c("#00007F","blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                                                                                                                                                      "red", "#7F0000")), backgroundValues = NULL)
{
  if (facet_by == "") {
    # no faceting
    if (subselection_col == "") {
      # no subselection
      color_values <- round(S4Vectors::metadata(sce)$SOM$map$medianValues[, marker], 2)
      colors <- colorPalette(100)[as.numeric(cut(color_values, breaks = 100))]
      # plot star chart for single marker with color_values and colors
      p <- plotStarsCustom(metadata(sce)$SOM, overall = FALSE, marker = marker, nodeValues = color_values, nodeColors = colorPalette, backgroundValues = backgroundValues)
      p <- p + ggtitle(marker) + theme(plot.title = element_text(hjust = 0.5))
      
      
      # TODO add legend for coloring
      #color_legend <- colorPalette(100)
      my_breaks <- seq(min(color_values), max(color_values), (max(color_values)-min(color_values))/5)
      legend_data <- data.frame(value = color_values, constant_y = 1)
      lplot <- ggplot(legend_data, aes(x = value, y = constant_y, fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colors = colorPalette(100), guide = "colourbar", breaks = my_breaks) +
        theme_classic() +
        theme(legend.key.height= unit(2,'cm'),
              legend.title = element_blank())
      # Draw Only Legend without plot 
      legend <- get_legend(lplot)
      
      # Add legend to the plot
      p <- p + theme(legend.position = "left")  # cluster legend
      p <- ggarrange(p, legend, ncol=2, widths = c(6, 1))
    } else {
      # subselection by group
      sce_filtered <- CATALYST::filterSCE(sce, get(subselection_col) == subselection)
      median_cond <- data.table::data.table(t(SummarizedExperiment::assay(sce_filtered, assayType)))
      median_cond[, cluster_id := sce_filtered$cluster_id]
      median_cond <- median_cond[, .(my_marker = median(get(marker))), by = cluster_id]
      missing_clusters <- median_cond[, setdiff(levels(cluster_id), cluster_id)]
      if (length(missing_clusters) != 0)
        median_cond <- rbind(median_cond, data.table::data.table(cluster_id = missing_clusters, my_marker= NA))
      data.table::setkey(median_cond, cluster_id)
      lev <- round(median_cond[, my_marker], 2)
      yval <- seq(min(lev, na.rm = T), max(lev, na.rm = T), by = (max(lev, na.rm = T) - min(lev, na.rm = T))/length(colorPalette(100)))
      
      # plot star chart for subselection with new calculated color_values
      colors <- colorPalette(100)[findInterval(median_cond[, my_marker], yval)]
      p <- plotStarsCustom(metadata(sce)$SOM, overall = FALSE, marker = marker, nodeValues = lev, nodeColors = colorPalette, backgroundValues = backgroundValues)
      p <- p + ggtitle(paste0(marker, ", ", subselection)) + theme(plot.title = element_text(hjust = 0.5))
      
      # TODO add legend for coloring
      my_breaks <- seq(min(lev, na.rm = T), max(lev, na.rm = T),by = (max(lev, na.rm = T) - min(lev, na.rm = T))/5)
      legend_data <- data.frame(value = my_breaks[-1], constant_y = 1)
      lplot <- ggplot(legend_data, aes(x = value, y = constant_y, fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colors = colorPalette(100), guide = "colourbar", breaks = my_breaks) +
        theme_classic() +
        theme(legend.key.height= unit(2, 'cm'),
              legend.title = element_blank())
      # Draw Only Legend without plot 
      legend <- get_legend(lplot)
      
      # Add legend to the plot
      p <- p + theme(legend.position = "left")  # cluster legend
      p <- ggarrange(p, legend, ncol=2, widths = c(6, 1))
    }
    return(p)
  } else {
    # facet
    metadata(sce)$experiment_info <- as.data.frame(metadata(sce)$experiment_info)
    metadata(sce)$experiment_info[[facet_by]] <- as.factor(metadata(sce)$experiment_info[[facet_by]])
    cond_levels <- levels(CATALYST::ei(sce)[[facet_by]])
    both_cond <- data.table::rbindlist(sapply(cond_levels, function(cond){
      if (subselection_col != ""){
        sce_filtered <- CATALYST::filterSCE(CATALYST::filterSCE(sce, marker_name == marker), get(facet_by) == cond, get(subselection_col) == subselection)
      }
      else {
        sce_filtered <- CATALYST::filterSCE(CATALYST::filterSCE(sce, marker_name == marker), get(facet_by) == cond)
      }
      median_cond <- data.table::data.table(t(SummarizedExperiment::assay(sce_filtered, assayType)))
      median_cond[, cluster_id := sce_filtered$cluster_id]
      median_cond <- median_cond[, .(my_marker = median(get(marker))), by = cluster_id]
      missing_clusters <- median_cond[, setdiff(levels(cluster_id), cluster_id)]
      if (length(missing_clusters) != 0)
        median_cond <- rbind(median_cond, data.table::data.table(cluster_id = missing_clusters, my_marker= NA))
      data.table::setkey(median_cond, cluster_id)
    }, simplify = FALSE), idcol = "condition")
    lev <- round(both_cond[, my_marker], 2)
    yval <- seq(min(lev, na.rm = T), max(lev, na.rm = T), by = (max(lev, na.rm = T) - min(lev, na.rm = T))/length(colorPalette(100)))
    plots <- list()
    for (cond in cond_levels) {
      # draw new plot with only selected condition
      color_values <- round(both_cond[condition == cond, my_marker], 2)
      colors <- colorPalette(100)[findInterval(color_values, yval)]
      p <- plotStarsCustom(metadata(sce)$SOM, overall = FALSE, marker = marker, nodeValues = color_values, nodeColors = colorPalette, backgroundValues = backgroundValues)
      if(subselection_col == ""){
        title_plot <- marker
      } else {
        title_plot <- paste0(marker, ", ", subselection)
      }
      p <- p + ggtitle(paste0(title_plot, ", ", cond)) + theme(plot.title = element_text(hjust = 0.5))
      plots[[cond]] <- p
    }
    # TODO add legend for coloring
    my_breaks <- seq(min(lev, na.rm = T), max(lev, na.rm = T),by = (max(lev, na.rm = T) - min(lev, na.rm = T))/5)
    legend_data <- data.frame(value = my_breaks[-1], constant_y = 1)
    lplot <- ggplot(legend_data, aes(x = value, y = constant_y, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colors = colorPalette(100), guide = "colourbar", breaks = my_breaks) +
      theme_classic() +
      theme(legend.key.height= unit(2, 'cm'),
            legend.title = element_blank())
    # Draw Only Legend without plot 
    legend <- get_legend(lplot)
    
    # hide one cluster legend to avoid double plotting
    plots[[1]] <- plots[[1]] + theme(legend.position = "left")
    plots[[2]] <- plots[[2]] + theme(legend.position = "none")
    
    # Combine the plots and legend with patchwork
    p <- wrap_plots(plots, nrow = 1, heights = c(6,1)) + legend
    return(p)
  }
}

addClusterAll <- function(sce){
  sce$cluster_id <- as.factor("all")
  S4Vectors::metadata(sce)$cluster_codes <- data.frame(all = as.factor("all"))
  return(sce)
}

