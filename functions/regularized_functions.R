library(glmnet)

sce <- readRDS("data/downsampled_files/sce_spiked_clustered_full_ds_15000.rds")

regularizedSCE <- function(sce, 
                           condition, 
                           features = SummarizedExperiment::rowData(sce)$marker_name,
                           assay_to_use = "exprs",
                           weighted = TRUE,
                           alpha=0){
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  features <-
    match.arg(features,
              SummarizedExperiment::rowData(sce)$marker_name,
              several.ok = TRUE)
  X <- t(SummarizedExperiment::assay(sce, assay_to_use))
  X <- X[, features]
  X <- (abs(X) + X) / 2
  y <- sce[["base_spike"]]
  fit <- glmnet::cv.glmnet(X, y, family = "binomial", alpha=0)  
  
}
