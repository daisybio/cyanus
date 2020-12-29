.get_features <- function (x, fs)
{
  if (is.null(fs)) {
    fs <- rownames(x)
  }
  else {
    stopifnot(is.character(fs))
    foo <- tryCatch(
      error = function(e)
        e,
      match.arg(fs,
                c("type", "state", "none"))
    )
    if (!inherits(foo, "error")) {
      fs <- foo
      stopifnot(!is.null(marker_classes(x)))
      fs <- rownames(x)[marker_classes(x) == fs]
      if (length(fs) == 0)
        stop("No features matched the specified marker class.")
    }
    else
      stopifnot(fs %in% rownames(x))
  }
  return(fs)
}

.triangle <- function (m) 
{
  n <- ncol(m)
  nm <- matrix(0, ncol = n, nrow = n)
  fm <- m
  nm[upper.tri(nm)] <- m[upper.tri(m)]
  fm <- t(nm) + nm
  diag(fm) <- diag(m)
  nm <- fm
  nm[upper.tri(nm)] <- NA
  diag(nm) <- NA
  m[lower.tri(nm)]
}


.plot_delta_area <- function (mc)
{
  library(ggplot2)
  maxK <- length(mc)
  v <- lapply(mc[seq_len(maxK)[-1]], function(x)
    .triangle(x$ml))
  h <- lapply(v, function(x) {
    h <- graphics::hist(x, breaks = seq(0, 1, 0.01), plot = FALSE)
    h$counts <- cumsum(h$counts) / sum(h$counts)
    return(h)
  })
  areaK <- vapply(h, function(x)
    cumsum(x$counts * 0.01)[100],
    numeric(1))
  deltaK <- c(areaK[1], diff(areaK) / areaK[seq_len(maxK - 2)])
  df <- data.frame(k = seq_len(maxK)[-1], y = deltaK)
  y_max <- ceiling(max(df$y) * 2) / 2
  ggplot(df, aes_string(x = "k", y = "y")) + theme_classic() +
    geom_line(color = "steelblue", lty = 2) + geom_point(size = 2.5,
                                                         color = "navy") + coord_fixed(4) + scale_x_continuous(breaks = seq(2,
                                                                                                                            20, 2), expand = c(0, 0.5)) + scale_y_continuous(
                                                                                                                              limits = c(0,
                                                                                                                                         y_max),
                                                                                                                              expand = c(0, 0.125),
                                                                                                                              breaks = function(x)
                                                                                                                                seq(x[1] +
                                                                                                                                      0.125, x[2], 0.5)
                                                                                                                            ) + ylab("Relative change\nin area under CDF curve") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey",
                                      size = 0.2)
    )
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
    features <- .get_features(x, features)
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
      fsom <- FlowSOM::ReadInput(flowCore::flowFrame(t(assay(x, assayType))))
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
      codes <- data.frame(seq_len(k), purrr::map(mc[-1], "consensusClass"))
      codes <-
        dplyr::mutate_all(codes, function(u)
          factor(u, levels = sort(unique(u))))
      colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s",
                                                        mcs))
      x$cluster_id <- factor(som$map$mapping[, 1])
      metadata(x)$cluster_codes <- codes
      metadata(x)$SOM_codes <- som$map$codes
      metadata(x)$delta_area <- .plot_delta_area(mc)
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
        cytofkit::ClusterX(t(assay(x, assayType)[features,]))
      colData(x)$clusterX_id <- clusterx$cluster
      metadata(x)$cluster_run$clusterX <- list(
        features = features,
        assayType = assays[assayType]
      )
    } else if (method == "rphenoGraph") {
      library(cytofkit)
      rphenograph <-
        cytofkit::Rphenograph(t(assay(x, assayType)[features,]), maxK)
      colData(x)$rphenoGraph_id[as.numeric(rphenograph$names)] <- rphenograph$membership
      metadata(x)$cluster_run$rphenoGraph <- list(
        features = features,
        assayType = assays[assayType],
        k = maxK
      )
    } else
      stop("which method was selected?")
    return(x)
  }