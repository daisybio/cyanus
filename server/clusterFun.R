plotClusterExprsCustom <-
  function (x,
            method = c("flowSOM", "clusterX", "rphenoGraph"),
            k = "meta20",
            features = "type")
  {
    library(data.table)
    library(ggplot2)
    
    method <- match.arg(method)
    CATALYST:::.check_sce(x)
    if (method == "flowSOM"){
      k <- CATALYST:::.check_k(x, k)
      x$cluster_id <- cluster_ids(x, k)
    } else x$cluster_id <- colData(x)[[sprintf("%s_id", method)]]
    features <- CATALYST:::.get_features(x, features)
    ms <- t(CATALYST:::.agg(x[features,], "cluster_id", "median"))
    d <- dist(ms, method = "euclidean")
    o <- hclust(d, method = "average")$order
    cd <- colData(x)
    es <- assay(x[features,], "exprs")
    df <- as.data.table(data.frame(t(es), cd, check.names = FALSE))
    df <- melt(
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
           )) + facet_wrap( ~ antigen, scales = "free_x",
                            nrow = 2) + ggridges::geom_density_ridges(alpha = 0.2) + ggridges::theme_ridges() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      )
  }

plotFreqHeatmapCustom <- function (x, 
                                   method = c("flowSOM", "clusterX", "rphenoGraph"),
                                   k = "meta20", m = NULL, normalize = TRUE, row_anno = TRUE, 
                                   col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, row_dend = TRUE, 
                                   col_dend = TRUE, bars = TRUE, perc = FALSE, hm_pal = rev(RColorBrewer::brewer.pal(11, 
                                                                                                       "RdBu")), k_pal = CATALYST:::.cluster_cols, m_pal = k_pal) 
{
  library(ComplexHeatmap)
  
  method <- match.arg(method)
  args <- as.list(environment())
  #CATALYST:::.check_args_plotFreqHeatmap(args)
  if (method == "flowSOM"){
    x$cluster_id <- cluster_ids(x, k)
  } else x$cluster_id <- colData(x)[[sprintf("%s_id", method)]]
  ns <- table(x$cluster_id, x$sample_id)
  fq <- prop.table(ns, 2)
  y <- as.matrix(unclass(fq))
  if (normalize) 
    y <- CATALYST:::.z_normalize(asin(sqrt(y)))
  if (!isFALSE(row_anno)) {
    left_anno <- CATALYST:::.anno_clusters(x, k, m, k_pal, m_pal)
  }
  else left_anno <- NULL
  if (!isFALSE(col_anno)) {
    top_anno <- CATALYST:::.anno_factors(x, levels(x$sample_id), col_anno, 
                              "colum")
  }
  else top_anno <- NULL
  if (bars) {
    right_anno <- CATALYST:::.anno_counts(x$cluster_id, perc)
  }
  else right_anno <- NULL
  ComplexHeatmap::Heatmap(matrix = y, name = paste0("normalized\n"[normalize], 
                                    "frequency"), col = hm_pal, na_col = "lightgrey", rect_gp = gpar(col = "white"), 
          column_title = "sample_id", column_title_side = "bottom", 
          cluster_rows = row_clust, cluster_columns = col_clust, 
          show_row_dend = row_dend, show_column_dend = col_dend, 
          show_row_names = is.null(left_anno), row_names_side = "left", 
          top_annotation = top_anno, left_annotation = left_anno, 
          right_annotation = right_anno)
}

clusterSCE <-
  function (x,
            method = c("flowSOM", "clusterX", "rphenoGraph"),
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
    method <- match.arg(method)
    features <- CATALYST:::.get_features(x, features)
    if (is.null(marker_classes(x))) {
      rowData(x)$marker_class <-
        factor(c("state", "type")[as.numeric(rownames(x) %in%
                                               features) + 1], levels = c("type", "state", "none"))
    }
    rowData(x)[[sprintf("used_for_clustering_%s", method)]] <-
      rownames(x) %in% features
    
    assays = c("exprs" = "Transformed", "counts" = "Raw")
    if (method == "flowSOM") {
      if (verbose)
        message("o running FlowSOM clustering...")
      fsom <-
        FlowSOM::ReadInput(flowCore::flowFrame(t(assay(x, assayType))))
      som <- FlowSOM::BuildSOM(
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
      metadata(x)$cluster_run$flowSOM <- list(
        features = features,
        assayType = assays[assayType],
        xdim = xdim,
        ydim = ydim,
        maxK = maxK
      )
    } else if (method == "clusterX") {
      library(cytofkit)
      clusterx <-
        cytofkit::ClusterX(t(assay(x, assayType)[features, ]))
      colData(x)$clusterX_id <- clusterx$cluster
      colData(x)$clusterX_id <- as.factor(colData(x)$clusterX_id)
      metadata(x)$cluster_run$clusterX <- list(features = features,
                                               assayType = assays[assayType])
    } else if (method == "rphenoGraph") {
      library(cytofkit)
      rphenograph <-
        cytofkit::Rphenograph(t(assay(x, assayType)[features, ]), maxK)
      colData(x)$rphenoGraph_id[as.numeric(rphenograph$names)] <-
        rphenograph$membership
      colData(x)$rphenoGraph_id <- as.factor(colData(x)$rphenoGraph_id)
      metadata(x)$cluster_run$rphenoGraph <- list(features = features,
                                                  assayType = assays[assayType],
                                                  k = maxK)
    } else
      stop("which method was selected?")
    return(x)
  }