sceGAMLSS <- function (sce, 
                    method = c("BEZI", "ZAGA", "ZAIG"),
                    condition, 
                    random_effect = NULL, 
                    assay_to_use = "exprs",
                    parallel = FALSE,
                    features = SummarizedExperiment::rowData(sce)$marker_name){
  
  library(gamlss)
  library(gamlss.dist)
  
  bppar <- BiocParallel::bpparam()
  if (!parallel)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  else 
    message('running in parallel')
  # Controls
  method <- match.arg(method, several.ok = FALSE)
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  features <- match.arg(features, SummarizedExperiment:: rowData(sce)$marker_name, several.ok = TRUE)
  
  # Method Choice
  if (method == "BEZI"){
    family <- gamlss.dist::BEZI()
    message("Fitting a zero-inflated beta distribution...")
  } else if (method == "ZAGA"){
    family <- gamlss.dist::ZAGA()
    message("Fitting a zero adjusted Gamma distribution...")
  } else if (method == "ZAIG"){
    family <- gamlss.dist::ZAIG()
    message("Fitting a zero adjusted Inverse Gaussian distribution...")
  } else stop("invalid family")
  
  # Data Preparation
  X <- t(SummarizedExperiment::assay(sce, assay_to_use))
  X <- X[, features, drop = FALSE]
  X <- (abs(X) + X) / 2
  X <-
    apply(
      X,
      MARGIN = 2,
      FUN = function(x)
        (x - min(x)) / diff(range(x))
    )
  X[which(X == 1)] <- X[which(X == 1)] - 2.225074e-10
  DF <- as.data.frame(X)
  DF$Y <- sce[[condition]]
  
  # Weights Incorporation
  # if (weighted){
  #   my_weights <- NULL
  #   message("Fitting model with weights...")
  # } else{
  #   my_weights <- rep(1/ei(sce)$n_cells, ei(sce)$n_cells)
  #   message("Fitting model without weights...")
  # }
  
  
  # Random Effects
  if (!is.null(random_effect)){
    for (re in random_effect){
      match.arg(re, names(SummarizedExperiment::colData(sce)))
    }
    R <- SummarizedExperiment::colData(sce)[random_effect]
    message(paste("...including random effects:", toString(random_effect)))
    DF <- cbind(DF, R)
  }
  
  # # Iteration over Markers
  data.table::rbindlist(BiocParallel::bplapply(features, function(marker, DF) {
    message(paste('fitting marker', marker))
    pval <- NA
    tryCatch({
      if (!is.null(random_effect)){
        random_terms <- list()
        for (re in random_effect){
          random_terms[[re]] <- ~ 1
        }
        bereg <- gamlss::gamlss(get(marker) ~ Y + re(random=random_terms), family = family,
                                trace = FALSE, control = gamlss.control(n.cyc = 200), data = DF, weights = NULL)
      }
      else
        bereg <- gamlss::gamlss(get(marker) ~ Y, family = family,
                                trace = FALSE, control = gamlss.control(n.cyc = 200), data = DF, weights = NULL)
      out <- summary(bereg)
      pval <-  out[2, 4]
    }, error=function(e) {
      message(e$message)
    })
    data.table::data.table(marker_id = marker, 
                           p_val = pval)
  }, DF, BPPARAM = bppar))
}
