diffcyt_method <- function (d_input, experiment_info = NULL, marker_info = NULL, 
          design = NULL, formula = NULL, contrast, analysis_type = c("DA", 
                                                                     "DS"), method_DA = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", 
                                                                                          "diffcyt-DA-GLMM"), method_DS = c("diffcyt-DS-limma", 
                                                                                                                            "diffcyt-DS-LMM"), markers_to_test = NULL, clustering_to_use = NULL, 
          cols_to_include = NULL, subsampling = FALSE, n_sub = NULL, 
          seed_sub = NULL, transform = TRUE, cofactor = 5, cols_clustering = NULL, 
          xdim = 10, ydim = 10, meta_clustering = FALSE, meta_k = 40, 
          seed_clustering = NULL, min_cells = 3, min_samples = NULL, 
          normalize = FALSE, norm_factors = "TMM", trend_method = "none", 
          block_id = NULL, trend = TRUE, use_weights = TRUE, plot = FALSE, 
          path = ".", verbose = TRUE) 
{
  analysis_type <- match.arg(analysis_type)
  method_DA <- match.arg(method_DA)
  method_DS <- match.arg(method_DS)
  if (!is(d_input, "SingleCellExperiment")) {
    if (is.null(experiment_info) | is.null(marker_info)) {
      stop("'experiment_info' and 'marker_info' must be provided (unless using a SingleCellExperiment ", 
           "object from CATALYST as input)")
    }
    if (verbose) 
      message("preparing data...")
    d_se <- prepareData(d_input, experiment_info, marker_info, 
                        cols_to_include, subsampling, n_sub, seed_sub)
    if (transform) {
      if (verbose) 
        message("transforming data...")
      d_se <- diffcyt::transformData(d_se, cofactor)
    }
    if (verbose) 
      message("generating clusters...")
    d_se <- generateClusters(d_se, cols_clustering, xdim, 
                             ydim, meta_clustering, meta_k, seed_clustering)
  }
  else if (is(d_input, "SingleCellExperiment")) {
    if (verbose) 
      message("using SingleCellExperiment object from CATALYST as input")
    if (is.null(clustering_to_use)) {
      stopifnot("cluster_id" %in% colnames(colData(d_input)))
      if (verbose) 
        message("using cluster IDs stored in column named 'cluster_id' in 'colData' of ", 
                "SingleCellExperiment object from CATALYST")
      clustering_name <- colnames(metadata(d_input)$cluster_codes)[1]
    }
    else if (!is.null(clustering_to_use)) {
      stopifnot(as.character(clustering_to_use) %in% colnames(metadata(d_input)$cluster_codes))
      stopifnot("cluster_id" %in% colnames(colData(d_input)))
      if (verbose) 
        message("using cluster IDs from clustering stored in column '", 
                clustering_to_use, "' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST")
      code_id <- colData(d_input)$cluster_id
      cluster_id <- metadata(d_input)$cluster_codes[, 
                                                    clustering_to_use][code_id]
      stopifnot(length(cluster_id) == nrow(colData(d_input)), 
                length(code_id) == nrow(colData(d_input)))
      colData(d_input)$code_id <- code_id
      colData(d_input)$cluster_id <- cluster_id
      clustering_name <- clustering_to_use
    }
    stopifnot("sample_id" %in% colnames(colData(d_input)))
    stopifnot("cluster_id" %in% colnames(colData(d_input)))
    stopifnot("cluster_codes" %in% names(metadata(d_input)))
    if (!("experiment_info" %in% names(metadata(d_input)))) {
      if (verbose) 
        message("generating 'experiment_info' from input object")
      m <- match(levels(droplevels(factor(d_input$sample_id))), 
                 d_input$sample_id)
      experiment_info <- data.frame(colData(d_input)[m, 
                                                     ], check.names = FALSE, row.names = NULL)
      metadata <- as.list(c(metadata(d_input), experiment_info))
    }
    else {
      experiment_info <- metadata(d_input)$experiment_info
      metadata <- metadata(d_input)
    }
    cs_by_s <- split(seq_len(ncol(d_input)), colData(d_input)$sample_id)
    cs <- unlist(cs_by_s[as.character(experiment_info$sample_id)])
    es <- t(assays(d_input)[["exprs"]])[cs, , drop = FALSE]
    d_se <- SummarizedExperiment(assays = list(exprs = es), 
                                 rowData = colData(d_input)[cs, ], colData = rowData(d_input), 
                                 metadata = metadata)
  }
  if (verbose) 
    message("calculating features...")
  d_counts <- calcCounts(d_se)
  if (analysis_type == "DS") {
    d_medians <- calcMedians(d_se)
    d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
    d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
  }
  if (analysis_type == "DA" && method_DA == "diffcyt-DA-edgeR") {
    if (verbose) 
      message("calculating DA tests using method 'diffcyt-DA-edgeR'...")
    res <- testDA_edgeR(d_counts, design, contrast, trend_method, 
                        min_cells, min_samples, normalize, norm_factors)
  }
  if (analysis_type == "DA" && method_DA == "diffcyt-DA-voom") {
    if (verbose) 
      message("calculating DA tests using method 'diffcyt-DA-voom'...")
    res <- testDA_voom(d_counts, design, contrast, block_id, 
                       min_cells, min_samples, normalize, norm_factors, 
                       plot, path)
  }
  if (analysis_type == "DA" && method_DA == "diffcyt-DA-GLMM") {
    if (verbose) 
      message("calculating DA tests using method 'diffcyt-DA-GLMM'...")
    res <- testDA_GLMM(d_counts, formula, contrast, min_cells, 
                       min_samples, normalize, norm_factors)
  }
  if (analysis_type == "DS" && method_DS == "diffcyt-DS-limma") {
    if (verbose) 
      message("calculating DS tests using method 'diffcyt-DS-limma'...")
    res <- diffcyt::testDS_limma(d_counts, d_medians, design, contrast, 
                        block_id, trend, use_weights, markers_to_test, min_cells, 
                        min_samples, plot, path)
  }
  if (analysis_type == "DS" && method_DS == "diffcyt-DS-LMM") {
    if (verbose) 
      message("calculating DS tests using method 'diffcyt-DS-LMM'...")
    res <- lmm_method(d_counts, d_medians, formula, contrast, 
                      use_weights, markers_to_test, min_cells, min_samples)
  }
  if (!is(d_input, "SingleCellExperiment")) {
    if (analysis_type == "DA") {
      return(list(res = res, d_se = d_se, d_counts = d_counts))
    }
    else if (analysis_type == "DS") {
      return(list(res = res, d_se = d_se, d_counts = d_counts, 
                  d_medians = d_medians, d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
                  d_medians_by_sample_marker = d_medians_by_sample_marker))
    }
  }
  else if (is(d_input, "SingleCellExperiment")) {
    metadata(res) <- as.list(c(metadata(res), clustering_name = clustering_name))
    if (analysis_type == "DA") {
      return(list(res = res, d_counts = d_counts))
    }
    else if (analysis_type == "DS") {
      return(list(res = res, d_counts = d_counts, d_medians = d_medians, 
                  d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
                  d_medians_by_sample_marker = d_medians_by_sample_marker))
    }
  }
}

lmm_method <- function (d_counts, d_medians, formula, contrast, use_weights = TRUE, 
                        markers_to_test = NULL, min_cells = 3, min_samples = NULL) 
{
  if (is.null(min_samples)) {
    min_samples <- ncol(d_counts)/2
  }
  if (!is.null(markers_to_test)) {
    markers_to_test <- markers_to_test
  }
  else {
    markers_to_test <- metadata(d_medians)$id_state_markers
  }
  counts <- assays(d_counts)[["counts"]]
  cluster_id <- rowData(d_counts)$cluster_id
  tf <- counts >= min_cells
  ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
  counts <- counts[ix_keep, , drop = FALSE]
  cluster_id <- cluster_id[ix_keep]
  if (is.logical(use_weights)) {
    if (use_weights) {
      weight <- colSums(counts)
    }
    else {
      weight <- NULL
    }
  }
  else if (is.numeric(weight)) {
    stopifnot(length(weight) == ncol(d_counts))
  }
  if (ncol(contrast) == 1 & nrow(contrast) > 1) {
    contrast <- t(contrast)
  }
  state_names <- names(assays(d_medians))[markers_to_test]
  meds <- do.call("rbind", {
    lapply(as.list(assays(d_medians)[state_names]), function(a) a[as.character(cluster_id), 
                                                                  , drop = FALSE])
  })
  meds_all <- do.call("rbind", as.list(assays(d_medians)[state_names]))
  p_vals <- rep(NA, nrow(meds))
  for (i in seq_len(nrow(meds))) {
    p_vals[i] <- tryCatch({
      y <- meds[i, ]
      data_i <- cbind(cbind(y, weight), formula$data)
      if (formula$random_terms) {
        if (!use_weights){
          fit <- lmer(formula$formula, data = data_i) #weights = NULL
        } else {
          fit <- lmer(formula$formula, data = data_i, weights = weight)
        }
        
      }
      else {
        if (!use_weights){
          fit <- lm(formula$formula, data = data_i) #weights=NULL
        } else {
          fit <- lm(formula$formula, data = data_i, weights = weight)
        }
      }
      test <- multcomp::glht(fit, contrast)
      summary(test)$test$pvalues
    }, error = function(e) NA)
  }
  p_adj <- p.adjust(p_vals, method = "fdr")
  stopifnot(length(p_vals) == length(p_adj))
  out <- data.frame(p_val = p_vals, p_adj = p_adj, stringsAsFactors = FALSE)
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster_id) * 
                                     length(state_names), ncol = ncol(out)))
  colnames(row_data) <- colnames(out)
  cluster_id_nm <- as.numeric(cluster_id)
  s <- seq(0, nlevels(cluster_id) * (length(state_names) - 
                                       1), by = nlevels(cluster_id))
  r1 <- rep(cluster_id_nm, length(state_names))
  r2 <- rep(s, each = length(cluster_id_nm))
  stopifnot(length(s) == length(state_names))
  stopifnot(length(r1) == length(r2))
  rows <- r1 + r2
  row_data[rows, ] <- out
  clus <- factor(rep(levels(cluster_id), length(state_names)), 
                 levels = levels(cluster_id))
  stat <- factor(rep(state_names, each = length(levels(cluster_id))), 
                 levels = state_names)
  stopifnot(length(clus) == nrow(row_data), length(stat) == 
              nrow(row_data))
  row_data <- cbind(data.frame(cluster_id = clus, marker_id = stat, 
                               stringsAsFactors = FALSE), row_data)
  col_data <- colData(d_medians)
  res <- SummarizedExperiment(meds_all, rowData = row_data, 
                              colData = col_data)
  res
}
