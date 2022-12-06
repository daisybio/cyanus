
plotStarsCustom <-
  function (sce, 
            backgroundValues = NULL,
            markers = SummarizedExperiment::rowData(sce)$marker_name[SummarizedExperiment::rowData(sce)$used_for_clustering], 
            colorPalette = grDevices::colorRampPalette(
              c(
                "#00007F",
                "blue",
                "#007FFF",
                "cyan",
                "#7FFF7F",
                "yellow",
                "#FF7F00",
                "red",
                "#7F0000"
              )
            ), 
            list_insteadof_ggarrange = FALSE, ...) 
  {
    if (length(names(list(...))) > 0 && "backgroundColor" %in% 
        names(list(...))) {
      warning(paste0("\"backgroundColor\" is deprecated, ", 
                     "please use \"backgroundColors\" instead."))
    }
    p <- PlotFlowSOMCustom(sce, view = "MST", 
                             backgroundValues = backgroundValues,
                             maxNodeSize = 1.5)
    channels <- rowData(sce)[rownames(rowData(sce)) %in% markers, "channel_name"]
    p <- AddStarsCustom(p = p, sce=sce, markers = markers, colorPalette = colorPalette)
    if (!is.null(names(colorPalette))) {
      names(colorPalette) <- markers
    }
    l1 <- FlowSOM::PlotStarLegend(markers, colorPalette)
    l2 <- ggpubr::get_legend(p)
    if (list_insteadof_ggarrange) {
      p <- p + ggplot2::theme(legend.position = "none")
      l2 <- ggpubr::as_ggplot(l2)
      return(list(tree = p, starLegend = l1, backgroundLegend = l2))
    }
    else {
      p <- ggpubr::ggarrange(p, ggpubr::ggarrange(l1, l2, 
                                                  ncol = 1), NULL, ncol = 3, widths = c(3, 1, 0.3), 
                             legend = "none")
      return(p)
    }
  }

PlotFlowSOMCustom <- function (sce, view = "MST", 
                               nodeSizes = S4Vectors::metadata(sce)$SOM_MST$size,
                               maxNodeSize = 1.5, 
                               refNodeSize = max(nodeSizes),
                               equalNodeSize = FALSE, 
                               backgroundValues = NULL,
                               backgroundColors = NULL, 
                               backgroundLim = NULL, 
                               title = NULL) 
{
  requireNamespace("ggplot2")
  nNodes <- nrow(metadata(sce)$SOM_MST$l)
  if (length(backgroundValues) != nNodes && !is.null(backgroundValues)) {
    stop(paste0("Length of 'backgroundValues' should be equal to number of ", 
                "clusters in FlowSOM object (", nNodes, " clusters and ", 
                length(backgroundValues), " backgroundValues)."))
  }
  if (deparse(substitute(backgroundColors)) == "backgroundColor") {
    warning(paste0("In the new FlowSOM version \"backgroundColors\"", 
                   " is used instead of \"backgroundColor\""))
  }
  layout <- as.data.frame(metadata(sce)$SOM_MST$l)
  colnames(layout) <- c("x", "y")
  autoNodeSize <- min(stats::dist(layout[, c(1, 2)]))
  maxNodeSize <- autoNodeSize * maxNodeSize
  if (equalNodeSize) {
    scaledNodeSize <- rep(maxNodeSize, nNodes)
  }
  else {
    scaledNodeSize <- data.frame(size = S4Vectors::metadata(sce)$SOM_MST$size) %>% 
      dplyr::mutate(scaled = (.data$size/refNodeSize) * maxNodeSize^2) %>% 
      dplyr::mutate(sqrt_scaled = sqrt(.data$scaled)) %>% 
      dplyr::pull(.data$sqrt_scaled)
  }
  plot_df <- data.frame(x = layout$x, y = layout$y, size = scaledNodeSize, 
                        bg_size = scaledNodeSize * 1.5)
  p <- ggplot2::ggplot(plot_df)
  if (!is.null(backgroundValues)) {
    p <- AddBackgroundCustom(p, backgroundValues = backgroundValues, 
                               backgroundColors = backgroundColors, backgroundLim = backgroundLim)
  }
  if (view == "MST") {
    edges <- ParseEdgesCustom(sce)
    p <- p + ggplot2::geom_segment(data = edges, ggplot2::aes(x = .data$x, 
                                                              y = .data$y, xend = .data$xend, yend = .data$yend), 
                                   linewidth = 0.2)
  }
  nodeInfo <- ggplot2::ggplot_build(p)$plot$data
  p <- p + ggforce::geom_circle(data = nodeInfo, 
                                ggplot2::aes(x0 = .data$x, y0 = .data$y, r = .data$size), 
                                fill = "white", 
                                size = 0.2)
  p <- p + ggplot2::coord_fixed() + ggplot2::theme_void()
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  return(p)
}


ParseEdgesCustom <- function (sce) 
{
  edgeList <- as.data.frame(igraph::as_edgelist(metadata(sce)$SOM_MST$graph), 
                            stringsAsFactors = FALSE)
  coords <- metadata(sce)$SOM_MST$l
  segmentPlot <- lapply(seq_len(nrow(edgeList)), function(row_id) {
    node_ids <- as.numeric(edgeList[row_id, ])
    row <- c(coords[node_ids[1], 1], coords[node_ids[1], 
                                            2], coords[node_ids[2], 1], coords[node_ids[2], 
                                                                               2])
    return(row)
  })
  segmentPlot <- do.call(rbind, segmentPlot)
  colnames(segmentPlot) <- c("x", "y", "xend", "yend")
  return(as.data.frame(segmentPlot))
}


AddStarsCustom <- function (p, 
                              sce, 
                              markers = rowData(sce)$marker_name, 
                              colorPalette = NULL) 
{
  nNodes <- nrow(metadata(sce)$SOM_MST$l)
  nMarkers <- length(markers)
  nodeInfo <- ggplot2::ggplot_build(p)$plot$data
  data <- S4Vectors::metadata(sce)$SOM_medianValues[, markers, drop = FALSE]
  data[is.na(data)] <- 0
  data <- FlowSOM:::ScaleStarHeights(data, nodeInfo$size)
  markers_tmp <- rep(1, ncol(data))
  names(markers_tmp) <- colnames(data)
  starValues <- lapply(seq_len(nNodes), function(cl) {
    nodeData <- FlowSOM:::ParseArcs(x = nodeInfo$x[cl], y = nodeInfo$y[cl], 
                                    arcValues = markers_tmp, arcHeights = data[cl, ])
    return(nodeData)
  })
  starValues <- do.call(rbind, starValues)
  starValues$Markers <- factor(starValues$Markers, levels = colnames(data))
  p <- FlowSOM:::AddStarsPies(p, starValues, colorPalette, showLegend = FALSE)
  return(p)
}

AddBackgroundCustom <- function (p, backgroundValues, backgroundColors = NULL, backgroundLim = NULL) 
{
  if (is.character(backgroundValues)) {
    backgroundValues <- factor(backgroundValues)
  }
  p <- FlowSOM:::AddScale(p, backgroundValues, backgroundColors, backgroundLim, 
                          labelLegend = "Background")
  p <- p + ggforce::geom_circle(ggplot2::aes(x0 = .data$x, 
                                             y0 = .data$y, r = .data$bg_size, fill = backgroundValues), 
                                col = NA, alpha = 0.4)
  return(p)
}



plotMarkerCustom <- function(sce,
                             marker = SummarizedExperiment::rowData(sce)$marker_name[SummarizedExperiment::rowData(sce)$used_for_clustering][1], 
                             facet_by = "", subselection_col = "",
                             subselection=NULL, assayType = "exprs",
                             refMarkers = SummarizedExperiment::rowData(sce)$marker_name[SummarizedExperiment::rowData(sce)$used_for_clustering],
                             colorPalette = grDevices::colorRampPalette(
                               c(
                                 "#00007F",
                                 "blue",
                                 "#007FFF",
                                 "cyan",
                                 "#7FFF7F",
                                 "yellow",
                                 "#FF7F00",
                                 "red",
                                 "#7F0000"
                               )
                             ),
                             backgroundValues = NULL,
                             lim = NULL, 
                             ...){
  mfis <- GetClusterMFIsCustom(sce)
  channels <- rowData(sce)[rowData(sce)$marker_name %in% marker, "channel_name"]
  marker_dict <- rowData(sce)[rowData(sce)$marker_name %in% marker, "marker_name"]
  names(marker_dict) <- rowData(sce)[rowData(sce)$marker_name %in% marker, "channel_name"]
  if(is.null(subselection) & facet_by == "" & subselection_col == ""){
    if (is.null(lim)) 
      lim <- c(min(mfis[, channels]), max(mfis[, channels]))
    plotList <- lapply(seq_along(channels), function(channelI) {
      p <- PlotVariableCustom(sce, variable = mfis[, channels[channelI]], 
                              variableName = "MFI", colorPalette = colorPalette, 
                              backgroundValues = backgroundValues,
                              lim = lim, ...)
      p <- p + ggplot2::ggtitle(marker_dict[channelI])
    })
  }
  else if(facet_by != ""){
      metadata(sce)$experiment_info[[facet_by]] <- as.factor(metadata(sce)$experiment_info[[facet_by]])
      cond_levels <- levels(CATALYST::ei(sce)[[facet_by]])
      mfis_cond <- data.table::rbindlist(sapply(cond_levels, function(cond){
        if (subselection_col != "") 
          sce_filtered <- CATALYST::filterSCE(CATALYST::filterSCE(sce, marker_name %in% marker), get(facet_by) == cond, get(subselection_col) == subselection)
        else 
          sce_filtered <- CATALYST::filterSCE(CATALYST::filterSCE(sce, marker_name %in% marker), get(facet_by) == cond)
        median_cond <- data.table::data.table(t(SummarizedExperiment::assay(sce_filtered, assayType)))
        median_cond[, cluster_id := sce_filtered$cluster_id]
        tmp_list <- lapply(marker, function(m){
          tmpDF <- median_cond[, .(median(get(m), na.rm = TRUE)), by = cluster_id]
          missing_clusters <- median_cond[, setdiff(levels(cluster_id), cluster_id)]
          if (length(missing_clusters) != 0)
            tmpDF <- rbind(tmpDF, data.table::data.table(cluster_id = missing_clusters, V1= NA))
          setnames(tmpDF, 'V1', m)
        })
        median_cond <- Reduce(merge, tmp_list)
        data.table::setkey(median_cond, cluster_id)
      }, simplify = FALSE), idcol = "condition")
      if (is.null(lim))
        lim <- c(min(mfis_cond[, ..marker], na.rm = TRUE), max(mfis_cond[, ..marker], na.rm = TRUE))
      plotList <- lapply(seq_along(channels), function(channelI) {
        lapply(cond_levels, function(cond){
          p <- PlotVariableCustom(sce, variable = mfis_cond[condition == cond, get(marker_dict[channelI])], 
                                  variableName = "MFI", colorPalette = colorPalette, 
                                  backgroundValues = backgroundValues,
                                  lim = lim, ...)
          if(subselection_col != ""){
            p <- p + ggplot2::ggtitle(paste(marker_dict[channelI], ',', cond, ',', subselection))
          }else{
            p <- p + ggplot2::ggtitle(paste(marker_dict[channelI], ',', cond))
          }
        })
      })
      plotList <- unlist(plotList, recursive = FALSE)
    }
  else{
    sce_filtered <- CATALYST::filterSCE(sce, get(subselection_col) == subselection)
    median_cond <- data.table::data.table(t(SummarizedExperiment::assay(sce_filtered, assayType)))
    median_cond[, cluster_id := sce_filtered$cluster_id]
    tmp_list <- lapply(marker, function(m){
      tmpDF <- median_cond[, .(median(get(m), na.rm = TRUE)), by = cluster_id]
      missing_clusters <- median_cond[, setdiff(levels(cluster_id), cluster_id)]
      if (length(missing_clusters) != 0)
        tmpDF <- rbind(tmpDF, data.table::data.table(cluster_id = missing_clusters, V1= NA))
      setnames(tmpDF, 'V1', m)
    })
    median_cond <- Reduce(merge, tmp_list)
    data.table::setkey(median_cond, cluster_id)
    plotList <- lapply(seq_along(channels), function(channelI) {
      p <- PlotVariableCustom(sce, variable = median_cond[, get(marker_dict[channelI])], 
                              variableName = "MFI", colorPalette = colorPalette, 
                              backgroundValues = backgroundValues,
                              lim = lim, ...)
      p <- p + ggplot2::ggtitle(paste(marker_dict[channelI], ',', subselection))
    })
  }
  p <- ggpubr::ggarrange(plotlist = plotList, common.legend = TRUE, 
                         legend = "right")
  return(p)
  
}


GetClusterMFIsCustom <- function (sce) 
{
  MFIs <-S4Vectors::metadata(sce)$SOM_medianValues[, drop = FALSE]
  rownames(MFIs) <- seq_len(nrow(MFIs))
  colnames(MFIs) <- rowData(sce)[, "channel_name"]
  return(MFIs)
}


PlotVariableCustom <- function (sce, variable, variableName = "", colorPalette = FlowSOM_colors, 
                                lim = NULL, ...) 
{
  if (length(variable) != nrow(metadata(sce)$SOM_MST$l)) {
    stop(paste0("Length of 'variable' should be equal to number of clusters in", 
                " FlowSOM object (", nrow(metadata(sce)$SOM_MST$l), " clusters and ", 
                length(variable), " variables)."))
  }
  p <- PlotFlowSOMCustom(sce = sce, ...)
  p <- FlowSOM::AddNodes(p = p, values = variable, colorPalette = colorPalette, 
                         lim = lim, label = variableName)
  return(p)
}


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
    p <- ggplot2::ggplot(df, ggplot2::aes(y = Freq)) + ggplot2::labs(x = NULL,
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

clusterSCE <-
  function (x,
            assayType = "exprs",
            features = "type",
            xdim = 10,
            ydim = 10,
            maxK = 20,
            verbose = TRUE,
            seed = 1)
  {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(
      is.logical(verbose),
      length(verbose) == 1,
      vapply(list(xdim,
                  ydim, maxK, seed), function(arg)
                    is.numeric(arg) &&
               length(arg) == 1, logical(1))
    )
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
    S4Vectors::metadata(x)$SOM_codes <- som$map$codes
    S4Vectors::metadata(x)$SOM_medianValues <- som$map$medianValues
    S4Vectors::metadata(x)$SOM_MST <- som$MST
    # this was changed because of the new FlowSOM clustering
    S4Vectors::metadata(x)$SOM_MST$size <- som$map$pctgs
    S4Vectors::metadata(x)$delta_area <- CATALYST:::.plot_delta_area(mc)
    x <-
      CATALYST::mergeClusters(
        x,
        k = sprintf("meta%s", maxK),
        id = "all",
        table = data.frame(old_cluster = seq_len(maxK), new_cluster = "all")
      )
    return(x)
  }

addClusterAll <- function(sce){
  sce$cluster_id <- as.factor("all")
  S4Vectors::metadata(sce)$cluster_codes <- data.frame(all = as.factor("all"))
  return(sce)
}

