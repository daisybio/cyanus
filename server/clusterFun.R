plotAbundancesCustom <-
  function (x,
            k = "meta20",
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
      nk <- nlevels(cluster_ids(x, k))
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
    p <- ggplot(df, aes_string(y = "Freq")) + labs(x = NULL,
                                                   y = "Proportion [%]") + theme_bw() + theme(
                                                     panel.grid = element_blank(),
                                                     strip.text = element_text(face = "bold"),
                                                     strip.background = element_rect(fill = NA,
                                                                                     color = NA),
                                                     axis.text = element_text(color = "black"),
                                                     axis.text.x = element_text(
                                                       angle = 45,
                                                       hjust = 1,
                                                       vjust = 1
                                                     ),
                                                     legend.key.height = unit(0.8, "lines")
                                                   )
    switch(
      by,
      sample_id = p + (if (!is.null(group_by))
        facet_wrap(group_by,
                   scales = "free_x")) + geom_bar(
                     aes_string(x = "sample_id",
                                fill = "cluster_id"),
                     position = "fill",
                     stat = "identity"
                   ) +
        scale_fill_manual("cluster_id", values = k_pal) + scale_x_discrete(expand = c(0,
                                                                                      0)) + scale_y_continuous(expand = c(0, 0), labels = seq(0,
                                                                                                                                              100, 25)) + theme(
                                                                                                                                                panel.border = element_blank(),
                                                                                                                                                panel.spacing.x = unit(1,
                                                                                                                                                                       "lines")
                                                                                                                                              ),
      cluster_id = {
        p <-
          p + scale_shape_manual(values = shapes) + guides(
            col = guide_legend(order = 1,
                               override.aes = list(size = 3)),
            shape = guide_legend(override.aes = list(size = 3))
          )
        if (is.null(group_by)) {
          p + geom_boxplot(
            aes_string(x = "cluster_id"),
            alpha = 0.2,
            position = position_dodge(),
            outlier.color = NA
          ) +
            geom_point(aes_string("cluster_id", shape = shape_by),
                       position = position_jitter(width = 0.2))
        } else {
          p + facet_wrap("cluster_id", scales = "free_y",
                         ncol = 4) + geom_boxplot(
                           aes_string(
                             x = group_by,
                             color = group_by,
                             fill = group_by
                           ),
                           position = position_dodge(),
                           alpha = 0.2,
                           outlier.color = NA,
                           show.legend = FALSE
                         ) +
            geom_point(aes_string(
              x = group_by,
              col = group_by,
              shape = shape_by
            ),
            position = position_jitter(width = 0.2))
        }
      }
    )
  }


plotClusterExprsCustom <-
  function (x,
            k = "meta20",
            features = "type",
            assay = "exprs")
  {
    CATALYST:::.check_sce(x, TRUE)
    k <- CATALYST:::.check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    features <- CATALYST:::.get_features(x, features)
    ms <-
      t(CATALYST:::.agg(x[features, ], "cluster_id", "median", assay = assay))
    d <- dist(ms, method = "euclidean")
    o <- hclust(d, method = "average")$order
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
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
    ggplot(df,
           aes_string(
             x = "expression",
             y = "cluster_id",
             col = "avg",
             fill = "avg"
           )) + facet_wrap(~ antigen, scales = "free_x",
                           nrow = 2) + ggridges::geom_density_ridges(alpha = 0.2) + ggridges::theme_ridges() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      )
  }

plotFreqHeatmapCustom <- function (x,
                                   k = "meta20",
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
  x$cluster_id <- cluster_ids(x, k)
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
            assayType,
            features = "type",
            xdim = 10,
            ydim = 10,
            maxK = 20,
            verbose = TRUE,
            seed = 1)
  {
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
    if (is.null(marker_classes(x))) {
      rowData(x)$marker_class <-
        factor(c("state", "type")[as.numeric(rownames(x) %in%
                                               features) + 1], levels = c("type", "state", "none"))
    }
    rowData(x)$used_for_clustering <- rownames(x) %in% features
    if (verbose)
      message("o running FlowSOM clustering...")
    fsom <-
      FlowSOM::ReadInput(flowCore::flowFrame(t(assay(x, assayType))))
    som <-
      FlowSOM::BuildSOM(
        fsom,
        colsToUse = features,
        silent = TRUE,
        xdim = xdim,
        ydim = ydim
      )
    if (verbose)
      message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <-
      suppressWarnings(suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(
          t(som$map$codes),
          maxK = maxK,
          reps = 100,
          distance = "euclidean",
          seed = seed,
          plot = NULL
        )
      ))
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
    metadata(x)$cluster_codes <- codes
    metadata(x)$SOM_codes <- som$map$codes
    metadata(x)$delta_area <- CATALYST:::.plot_delta_area(mc)
    return(x)
  }
