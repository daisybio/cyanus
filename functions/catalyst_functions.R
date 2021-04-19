library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)
library(ggvenn)



# function anno_features from CATALYST with different color palette
.anno_factors <- function (x, ids, which, type = c("row", "column")) 
{
  type <- match.arg(type)
  cd <- colData(x)
  df <- data.frame(cd, check.names = FALSE)
  df <- select_if(df, ~!is.numeric(.))
  df <- mutate_all(df, ~droplevels(factor(.x)))
  m <- match(ids, df$sample_id)
  ns <- split(df, df$sample_id) %>% lapply(mutate_all, droplevels) %>% 
    lapply(summarize_all, nlevels) %>% do.call(what = "rbind")
  keep <- names(which(colMeans(ns) == 1))
  keep <- setdiff(keep, c("sample_id", "cluster_id"))
  if (is.character(which)) 
    keep <- intersect(keep, which)
  if (length(keep) == 0) 
    return(NULL)
  df <- df[m, keep, drop = FALSE]
  lvls <- lapply(as.list(df), levels)
  nlvls <- vapply(lvls, length, numeric(1))
  #pal <- brewer.pal(8, "Set2")
  pal <-  c("#999999","#009E73","#E69F00", "#56B4E9", "#F0E442","#D55E00","#0072B2", "#CC79A7")
  names(is) <- is <- colnames(df)
  cols <- lapply(is, function(i) {
    if (nlvls[i] > length(pal)) 
      pal_i <- colorRampPalette(pal)(max(nlvls))
    else pal_i <- pal
    u <- pal_i[seq_len(nlvls[i])]
    names(u) <- lvls[[i]]
    u
  })
  ComplexHeatmap::HeatmapAnnotation(which = type, df = df, col = cols, gp = gpar(col = "white"))
}

plotExprHeatmap <- function (x, features = NULL, by = c("sample_id", "cluster_id", 
                                                        "both"), k = "meta20", m = NULL, assay = "exprs", fun = c("median", 
                                                                                                                  "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, 
                             row_anno = TRUE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
                             row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, 
                             bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
                             k_pal = CATALYST:::.cluster_cols, m_pal = k_pal, distance = c("euclidean", 
                                                                                           "maximum", "manhattan", "canberra", "binary", "minkowski"), 
                             linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
                                         "median", "centroid", "ward.D2")) 
{
  args <- as.list(environment())
  CATALYST:::.check_args_plotExprHeatmap(args)
  distance <- match.arg(distance)
  linkage <- match.arg(linkage)
  scale <- match.arg(scale)
  fun <- match.arg(fun)
  by <- match.arg(by)
  x <- x[unique(CATALYST:::.get_features(x, features)), ]
  if (by != "sample_id") {
    CATALYST:::.check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
  }
  if (by == "both") 
    by <- c("cluster_id", "sample_id")
  .do_agg <- function() {
    z <- CATALYST:::.agg(x, by, fun, assay)
    if (length(by) == 1) 
      return(z)
    set_rownames(do.call("rbind", z), levels(x$cluster_id))
  }
  .do_scale <- function() {
    if (scale == "first") {
      z <- assay(x, assay)
      z <- CATALYST:::.scale_exprs(z, 1, q)
      assay(x, assay, FALSE) <- z
      return(x)
    }
    else CATALYST:::.scale_exprs(z, 1, q)
  }
  z <- switch(scale, first = {
    x <- .do_scale()
    .do_agg()
  }, last = {
    z <- .do_agg()
    .do_scale()
  }, never = {
    .do_agg()
  })
  if (length(by) == 1) 
    z <- t(z)
  if (scale != "never" && !(assay == "counts" && fun == "sum")) {
    qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  }
  else lgd_aes <- list()
  lgd_aes$title_gp <- grid::gpar(fontsize = 10, fontface = "bold", 
                                 lineheight = 0.8)
  if (!isFALSE(row_anno)) {
    left_anno <- switch(by[1], sample_id = .anno_factors(x, 
                                                         levels(x$sample_id), row_anno, "row"), .anno_clusters(x, 
                                                                                                               k, m, k_pal, m_pal))
  }
  else left_anno <- NULL
  if (!isFALSE(col_anno) && length(by) == 2) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, 
                              "colum")
  }
  else top_anno <- NULL
  if (bars) {
    right_anno <- CATALYST:::.anno_counts(x[[by[1]]], perc)
  }
  else right_anno <- NULL
  if (bin_anno) {
    cell_fun <- function(j, i, x, y, ...) grid.text(gp = grid::gpar(fontsize = 8), 
                                                    sprintf("%.2f", z[i, j]), x, y)
  }
  else cell_fun <- NULL
  a <- ifelse(assay == "exprs", "expression", assay)
  f <- switch(fun, median = "med", fun)
  hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
                                            "scaled", a), last = sprintf("%s %s\n%s", "scaled", 
                                                                         fun, a), never = paste(fun, a, sep = "\n"))
  if (length(by) == 2) {
    col_title <- features
  }
  else if (length(features) == 1 && features %in% c("type", 
                                                    "state")) {
    col_title <- paste0(features, "_markers")
  }
  else col_title <- ""
  ComplexHeatmap::Heatmap(matrix = z, name = hm_title, col = circlize::colorRamp2(seq(min(z), 
                                                                                      max(z), l = n <- 100), grDevices::colorRampPalette(hm_pal)(n)), 
                          column_title = col_title, column_title_side = ifelse(length(by) == 
                                                                                 2, "top", "bottom"), cell_fun = cell_fun, cluster_rows = row_clust, 
                          cluster_columns = col_clust, show_row_dend = row_dend, 
                          show_column_dend = col_dend, clustering_distance_rows = distance, 
                          clustering_method_rows = linkage, clustering_distance_columns = distance, 
                          clustering_method_columns = linkage, show_row_names = (is.null(left_anno) || 
                                                                                   isTRUE(by == "sample_id")) && !perc, row_names_side = ifelse(by[1] == 
                                                                                                                                                  "cluster_id" || isFALSE(row_anno) && !row_dend || 
                                                                                                                                                  isFALSE(row_clust), "left", "right"), top_annotation = top_anno, 
                          left_annotation = left_anno, right_annotation = right_anno, 
                          rect_gp = gpar(col = "white"), heatmap_legend_param = lgd_aes)
}

plotDiffHeatmap <- function (x, y, k = NULL, top_n = 20, fdr = 0.05, lfc = 1, all = FALSE, 
                             sort_by = c("padj", "lfc", "none"), y_cols = list(padj = "p_adj", 
                                                                               lfc = "logFC", target = "marker_id"), assay = "exprs", 
                             fun = c("median", "mean", "sum"), normalize = TRUE, col_anno = TRUE, 
                             row_anno = TRUE, hm_pal = NULL, fdr_pal = c("lightgrey", 
                                                                         "lightgreen"), lfc_pal = c("blue3", "white", "red3")) 
{
  fun <- match.arg(fun)
  sort_by <- match.arg(sort_by)
  args <- as.list(environment())
  defs <- as.list(formals("plotDiffHeatmap")$y_cols[-1])
  miss <- !names(defs) %in% names(args$y_cols)
  if (any(miss)) 
    y_cols <- args$y_cols <- c(args$y_cols, defs[miss])[names(defs)]
  .check_args_plotDiffHeatmap(args)
  stopifnot(y_cols[[sort_by]] %in% names(y))
  y_cols <- y_cols[y_cols %in% names(y)]
  if (is.null(k)) {
    kids <- levels(y$cluster_id)
    same <- vapply(cluster_codes(x), function(u) identical(levels(u), 
                                                           kids), logical(1))
    if (!any(same)) 
      stop("Couldn't match any clustering", " in input data 'x' with results in 'y'.")
    k <- names(cluster_codes(x))[same][1]
  }
  else {
    k <- .check_k(x, k)
  }
  x$cluster_id <- cluster_ids(x, k)
  y <- data.frame(y, check.names = FALSE)
  y <- mutate_if(y, is.factor, as.character)
  if (any(rownames(x) %in% unlist(y))) {
    features <- intersect(rownames(x), y[[y_cols$target]])
    if (length(features) == 0) 
      stop("Couldn't match features between", " results 'y' and input data 'x'.")
    i <- y[[y_cols$target]] %in% rownames(x)
    type <- "ds"
  }
  else {
    i <- TRUE
    type <- "da"
  }
  y <- rename(y, target = y_cols$target, padj = y_cols$padj, 
              lfc = y_cols$lfc)
  i <- i & !is.na(y$padj) & y$cluster_id %in% levels(x$cluster_id)
  if (!all) {
    i <- i & y$padj < fdr
    if (!is.null(y$lfc)) 
      i <- i & abs(y$lfc) > lfc
  }
  y <- y[i, , drop = FALSE]
  if (nrow(y) == 0) 
    stop("No results remaining;", " perhaps 'x' or 'y' has been filtered,", 
         " or features couldn't be matched.")
  if (sort_by != "none") {
    o <- order(abs(y[[sort_by]]), decreasing = (sort_by == 
                                                  "lfc"))
    y <- y[o, , drop = FALSE]
  }
  if (top_n > nrow(y)) 
    top_n <- nrow(y)
  top <- y[seq_len(top_n), ]
  if (!isFALSE(col_anno)) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, 
                              "column")
  }
  else top_anno <- NULL
  if (is.null(hm_pal)) 
    hm_pal <- rev(brewer.pal(11, ifelse(type == "ds", "RdYlBu", 
                                        "RdBu")))
  if (row_anno) {
    s <- factor(ifelse(top$padj < fdr, "yes", "no"), levels = c("no", 
                                                                "yes"))
    if (!is.null(top$lfc)) {
      lfc_lims <- range(top$lfc, na.rm = TRUE)
      if (all(lfc_lims > 0)) {
        lfc_brks <- c(0, lfc_lims[2])
        lfc_pal <- lfc_pal[-1]
      }
      else if (all(lfc_lims < 0)) {
        lfc_brks <- c(lfc_lims[1], 0)
        lfc_pal <- lfc_pal[-3]
      }
      else lfc_brks <- c(lfc_lims[1], 0, lfc_lims[2])
      lfc_anno <- top$lfc
      anno_cols <- list(logFC = colorRamp2(lfc_brks, lfc_pal))
    }
    else {
      lfc_anno <- NULL
      anno_cols <- list()
    }
    names(fdr_pal) <- levels(s)
    anno_cols$significant <- fdr_pal
    right_anno <- rowAnnotation(logFC = lfc_anno, significant = s, 
                                foo = row_anno_text(scientific(top$padj, 2), gp = gpar(fontsize = 8)), 
                                col = anno_cols, gp = gpar(col = "white"), show_annotation_name = FALSE, 
                                simple_anno_size = unit(4, "mm"))
  }
  else right_anno <- NULL
  switch(type, da = {
    ns <- table(x$cluster_id, x$sample_id)
    fq <- prop.table(ns, 2)
    fq <- fq[top$cluster_id, ]
    y <- as.matrix(unclass(fq))
    if (normalize) y <- .z_normalize(asin(sqrt(y)))
    Heatmap(matrix = y, name = paste0("normalized\n"[normalize], 
                                      "frequency"), col = hm_pal, na_col = "lightgrey", 
            rect_gp = gpar(col = "white"), cluster_rows = FALSE, 
            cluster_columns = FALSE, row_names_side = "left", 
            top_annotation = top_anno, right_annotation = right_anno)
  }, ds = {
    y <- assay(x, assay)
    cs <- .split_cells(x, c("cluster_id", "sample_id"))
    z <- t(mapply(function(k, g) vapply(cs[[k]], function(cs) {
      if (length(cs) == 0) return(NA)
      get(fun)(y[g, cs, drop = FALSE])
    }, numeric(1)), k = top$cluster_id, g = top$target))
    rownames(z) <- sprintf("%s(%s)", top$target, top$cluster_id)
    if (normalize) z <- .z_normalize(z)
    Heatmap(matrix = z, name = paste0("z-normalized\n"[normalize], 
                                      "expression"), col = hm_pal, cluster_rows = FALSE, 
            cluster_columns = FALSE, top_annotation = top_anno, 
            row_names_side = "left", rect_gp = gpar(col = "white"), 
            right_annotation = right_anno, heatmap_legend_param = list(title_gp = gpar(fontsize = 10, 
                                                                                       fontface = "bold", lineheight = 0.8)))
  })
}