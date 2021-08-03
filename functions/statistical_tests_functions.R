

median_test <- function(test = c("Wilcoxon", "Kruskal_Wallis", "T_Test"),
                        sce,
                        condition,
                        random_effect = NULL,
                        features = SummarizedExperiment::rowData(sce)$marker_name,
                        assay_to_use = "exprs", 
                        parallel = FALSE){
  # Controls
  test <- match.arg(test)
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  features <- match.arg(features, SummarizedExperiment::rowData(sce)$marker_name, several.ok = TRUE)
  
  #parallelization
  stopifnot(is.logical(parallel))
  
  bppar <- BiocParallel::bpparam()
  if (!parallel)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  
  # Data Preparation
  X <- t(as.data.frame(SummarizedExperiment::assay(sce, assay_to_use)))
  data <- data.table::data.table(X[, features], y = SummarizedExperiment::colData(sce)[[condition]], 
                                 sample_id = SummarizedExperiment::colData(sce)[["sample_id"]])
  if (!is.null(random_effect)) {
    match.arg(random_effect, names(SummarizedExperiment::colData(sce)))
    data$re <- SummarizedExperiment::colData(sce)[[random_effect]]
  }
  
  p.values <- data.table::rbindlist(BiocParallel::bplapply(features, function(f){
    # medians_control <- data[y == levels(data$y)[1], median(get(f)), by = sample_id]
    # medians_case <- data[y == levels(data$y)[2], median(get(f)), by = sample_id]
    # median_data <- data.table::data.table(
    #   sample_id = medians_control$sample_id,
    #   median_control = medians_control$V1,
    #   median_case = medians_case$V1
    # )
    # median_data <- data.table::melt(median_data, id.vars = "sample_id", 
    #                                 variable.name = "condition", 
    #                                 value.name = "median")
    if (!is.null(random_effect)) {
      median_data <- data[, .(sample_median = median(get(f))), by = .(sample_id, y, re)]
      paired <- TRUE
    } else {
      median_data <- data[, .(sample_median = median(get(f))), by = .(sample_id, y)]
      paired <- FALSE
    }
    if(test == "Wilcoxon"){
      test_result <- wilcox.test(sample_median ~ y, data=median_data, paired = paired)
    }else if(test == "Kruskal_Wallis"){
      test_result <- kruskal.test(sample_median ~ y, data=median_data)
    } else if (test == "T_Test") {
      test_result <- t.test(sample_median ~ y, data=median_data, paired = paired)
    } else stop(sprintf('unknown test: found "%s"', test))
    return(
      data.table::data.table(
        marker_id = f,
        p_val = ifelse(is.nan(test_result$p.value), 1.0, test_result$p.value)
      )
    )
  }, BPPARAM = bppar))
  
  return(p.values)
}

median_wilcoxon_test <- function(sce,
                                 condition,
                                 features = SummarizedExperiment::rowData(sce)$marker_name,
                                 assay_to_use = "exprs", 
                                 parallel = FALSE){
  
  return(median_test(test = "Wilcoxon", 
                     sce,
                     condition,
                     features = SummarizedExperiment::rowData(sce)$marker_name,
                     assay_to_use = "exprs", 
                     parallel = parallel))
}

median_kruskal_test <- function(sce,
                                condition,
                                features = SummarizedExperiment::rowData(sce)$marker_name,
                                assay_to_use = "exprs", 
                                parallel = FALSE){
  
  return(median_test(test = "Kruskal_Wallis", 
                     sce,
                     condition,
                     features = SummarizedExperiment::rowData(sce)$marker_name,
                     assay_to_use = "exprs", 
                     parallel = parallel))

}


exprs_test <- function(test = c("Wilcoxon, Kruskal_Wallis", "t-test"),
                        sce,
                        condition,
                        features = SummarizedExperiment::rowData(sce)$marker_name,
                        assay_to_use = "exprs", 
                        parallel = FALSE){
  
  # Controls
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  features <- match.arg(features, SummarizedExperiment:: rowData(sce)$marker_name, several.ok = TRUE)
  
  #parallelization
  stopifnot(is.logical(parallel))
  
  bppar <- BiocParallel::bpparam()
  if (!parallel)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  
  # Data Preparation
  X <- t(SummarizedExperiment::assay(sce, assay_to_use))
  data <- data.table::data.table(X[, features], y = SummarizedExperiment::colData(sce)[[condition]])
  
  p.values <- data.table::rbindlist(BiocParallel::bplapply(features, function(f){
    
    if(test == "Wilcoxon"){
      test_result <- wilcox.test(x = data[y == levels(y)[1], get(f)], y = data[y == levels(y)[2], get(f)])
    }else if(test == "t-test"){
      test_result <- t.test(x = data[y == levels(y)[1], get(f)], y = data[y == levels(y)[2], get(f)])
    }else if(test == "Kruskal_Wallis"){
      test_result <- kruskal.test(x = list(control = data[y == levels(y)[1], get(f)], case = data[y == levels(y)[2], get(f)]))
    }
    return(
      data.table::data.table(
        marker_id = f,
        p_val = ifelse(is.nan(test_result$p.value), 1.0, test_result$p.value)
      )
    )
  }, BPPARAM = bppar))
  
  return(p.values)
  
}

exprs_wilcoxon_test <- function(sce,
                                 condition,
                                 features = SummarizedExperiment::rowData(sce)$marker_name,
                                 assay_to_use = "exprs", 
                                 parallel = FALSE){
  
  return(exprs_test(test = "Wilcoxon", 
                     sce,
                     condition,
                     features = SummarizedExperiment::rowData(sce)$marker_name,
                     assay_to_use = "exprs", 
                     parallel = parallel))
}

exprs_kruskal_test <- function(sce,
                                condition,
                                features = SummarizedExperiment::rowData(sce)$marker_name,
                                assay_to_use = "exprs", 
                                parallel = FALSE){
  
  return(exprs_test(test = "Kruskal_Wallis", 
                     sce,
                     condition,
                     features = SummarizedExperiment::rowData(sce)$marker_name,
                     assay_to_use = "exprs", 
                     parallel = parallel))
  
}

exprs_t_test <- function(sce,
                               condition,
                               features = SummarizedExperiment::rowData(sce)$marker_name,
                               assay_to_use = "exprs", 
                               parallel = FALSE){
  
  return(exprs_test(test = "t-test", 
                    sce,
                    condition,
                    features = SummarizedExperiment::rowData(sce)$marker_name,
                    assay_to_use = "exprs", 
                    parallel = parallel))
  
}


calculate_effect_size <- function(sce,
                                  condition,
                                  features = SummarizedExperiment::rowData(sce)$marker_name,
                                  assay_to_use = "exprs", 
                                  parallel = FALSE){
  # Controls
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  features <- match.arg(features, SummarizedExperiment::rowData(sce)$marker_name, several.ok = TRUE)
  
  #parallelization
  stopifnot(is.logical(parallel))
  
  bppar <- BiocParallel::bpparam()
  if (!parallel)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  
  # Data Preparation
  X <- t(as.data.frame(SummarizedExperiment::assay(sce, assay_to_use)))
  data <- data.table::data.table(X[, features], y = SummarizedExperiment::colData(sce)[[condition]])
  #effect sizes
  effect_sizes <- data.table::rbindlist(BiocParallel::bplapply(features, function(f){
    meanA <- data[y == levels(data$y)[1], mean(get(f))]
    meanB <- data[y == levels(data$y)[2], mean(get(f))]
    std <- data[, sd(get(f))]
    return(data.table::data.table(
      marker_id = f,
      effect_size = abs(meanA - meanB) / std
    ))
  }, BPPARAM = bppar))
  return(effect_sizes)
}

