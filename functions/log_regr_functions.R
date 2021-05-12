

logistic_regression <- function(sce,
                                condition,
                                random_effect = NULL,
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
  X <- t(as.data.frame(SummarizedExperiment::assay(sce, assay_to_use)))
  data <- data.frame(X[, features], y = SummarizedExperiment::colData(sce)[[condition]])
  
  if(is.null(random_effect)){
    p.values <- data.table::rbindlist(BiocParallel::bplapply(features, function(f){
      fit <- glm(as.formula(paste("y ~ ", f)), family=binomial, data)
      return(data.table::data.table(
        marker_id = f,
        p_val = coef(summary(fit))[2,4]
      ))
    }, BPPARAM = bppar))
  }else{
    for (re in random_effect){
      match.arg(re, names(SummarizedExperiment::colData(sce)))
    }
    R <- SummarizedExperiment::colData(sce)[random_effect]
    message(paste("...including random effects:", toString(random_effect)))
    data <- cbind(data, R)
    p.values <- data.table::rbindlist(BiocParallel::bplapply(features, function(f){
      rand_formula <- paste(paste0("( 1 | ", random_effect, ")"), 
                                          collapse = " + ")
      formula <- as.formula(paste("y ~ ", f, "+", rand_formula))
      fit <- lme4::glmer(formula = formula, family=binomial, data = data)
      return(data.table::data.table(
        marker_id = f,
        p_val = coef(summary(fit))[2,4]
      ))
    }, BPPARAM = bppar))
  }
  return(p.values)
}


