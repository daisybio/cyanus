zibSeq <- function (sce, 
                    condition, 
                    random_effect = NULL, 
                    weighted = FALSE, 
                    assay_to_use = "exprs",
                    features = SummarizedExperiment::rowData(sce)$marker_name){
  
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  
  features <- match.arg(features, SummarizedExperiment:: rowData(sce)$marker_name, several.ok = TRUE)
  
  if (!is.null(random_effect)){
    match.arg(random_effect, names(SummarizedExperiment::colData(sce)))
    R = SummarizedExperiment::colData(sce)[[random_effect]]
    message("Fitting a zero-inflated beta mixed model with random effects...")
  } else {
    message("Fitting a zero-inflated beta model...")
  }
  
  X <- t(as.data.frame(SummarizedExperiment::assay(sce, assay_to_use)))
  X <- X[, features]
  Y <- colData(sce)[[condition]]
  
  # inverse weights
  if (weighted){
    my_weights <- NULL
    message("Fitting model with weights")
  } else{
    my_weights <- rep(1/ei(sce)$n_cells, ei(sce)$n_cells)
    message("Fitting model without weights")
  }

  beta = matrix(data = NA, length(features), 2)
  for (i in seq(length(features))) {
    message(paste("Fitting marker", colnames(X)[i]))

    x.prop <- X[, i]
    x.prop[x.prop < 0] <- 0
    x.prop <- (x.prop - min(x.prop))/(max(x.prop)-min(x.prop))
    x.prop[which(x.prop==1)] <- x.prop[which(x.prop==1)] - 2.225074e-10
    data <- data.table(
      exprs = x.prop,
      condition = Y
    )
    
    if (!is.null(random_effect)){ 
      # with random effect
      data$random = R
      bereg = gamlss::gamlss(exprs ~ condition + re(random=(~1|random)), family = BEZI(), 
                             trace = FALSE, control = gamlss.control(n.cyc = 100), data = data, weights = my_weights)
    } else {
      # without random effect
      bereg = gamlss::gamlss(exprs ~ condition, family = BEZI(), 
                             trace = FALSE, control = gamlss.control(n.cyc = 100), data = data, weights = my_weights)
    }
    out = summary(bereg)
    beta[i, ] = out[2, c(1, 4)]
  }
  
  pvalues = beta[, 2]
  padj <- p.adjust(pvalues, method = "BH")
  res <- data.table(
    marker_id = colnames(X), 
    p_val = pvalues, 
    p_adj = padj
  )
  return (res)
}
