sceGAMLSS <- function (sce, 
                    method = c("BEZI", "ZAGA", "ZAIG"),
                    condition, 
                    random_effect = NULL, 
                    assay_to_use = "exprs",
                    features = SummarizedExperiment::rowData(sce)$marker_name){
  
  library(gamlss)
  library(gamlss.dist)
  
  # Controls
  method <- match.arg(method, several.ok = FALSE)
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  features <- match.arg(features, SummarizedExperiment:: rowData(sce)$marker_name, several.ok = TRUE)
  
  # Data Preparation
  X <- t(as.data.frame(SummarizedExperiment::assay(sce, assay_to_use)))
  X <- X[, features]
  Y <- colData(sce)[[condition]]
  
  
  # Weights Incorporation
  # if (weighted){
  #   my_weights <- NULL
  #   message("Fitting model with weights...")
  # } else{
  #   my_weights <- rep(1/ei(sce)$n_cells, ei(sce)$n_cells)
  #   message("Fitting model without weights...")
  # }
  
  
  # Method Choice
  if (method == "BEZI"){
    family <- BEZI()
    message("Fitting a zero-inflated beta distribution...")
  } else if (method == "ZAGA"){
    family <- ZAGA()
    message("Fitting a zero adjusted Gamma distribution...")
  } else if (method == "ZAIG"){
    family <- ZAIG()
    message("Fitting a zero adjusted Inverse Gaussian distribution...")
  }
  
  # Random Effects
  if (!is.null(random_effect)){
    for (re in random_effect){
      match.arg(re, names(SummarizedExperiment::colData(sce)))
    }
    R <- SummarizedExperiment::colData(sce)[random_effect]
    message(paste("...including random effects:", toString(random_effect)))
  }
  
  # Iteration over Markers
  beta = matrix(data = NA, length(features), 2)
  for (i in seq(length(features))) {
    
    message(paste("...fitting marker", colnames(X)[i]))
    x.prop <- X[, i]
    x.prop[x.prop < 0] <- 0
    x.prop <- (x.prop - min(x.prop))/(max(x.prop)-min(x.prop))
    x.prop[which(x.prop==1)] <- x.prop[which(x.prop==1)] - 2.225074e-10
    data <- data.frame(
      exprs = x.prop,
      condition = Y
    )
    tryCatch({
    if (!is.null(random_effect)){ 
      # with random effect
      data <- cbind(data, R)
      random_terms <- list()
      for (re in random_effect){
        random_terms[[re]] <- ~ 1
      }
      bereg = gamlss::gamlss(exprs ~ condition + re(random=random_terms), family = family, 
                             trace = FALSE, control = gamlss.control(n.cyc = 200), data = data, weights = NULL)
    } else {
      # without random effect
      bereg = gamlss::gamlss(exprs ~ condition, family = family, 
                             trace = FALSE, control = gamlss.control(n.cyc = 200), data = data, weights = NULL)
    }
    out = summary(bereg)
    beta[i, ] = out[2, c(1, 4)]
    }, error=function(e){
      message(e)
      beta[i, ] = c(NA, NA)
    }
    )
  }
  
  pvalues = beta[, 2]
  res <- data.frame(
    marker_id = colnames(X), 
    p_val = pvalues 
  )
  return (res)
}
