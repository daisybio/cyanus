hurdleBeta <- function(sce,
                       condition,
                       random_effect = NULL,
                       features = SummarizedExperiment::rowData(sce)$marker_name,
                       assay_to_use = "exprs",
                       weighted = TRUE,
                       parallel = FALSE) {
  
  bppar <- BiocParallel::bpparam()
  if (!parallel)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
  match.arg(condition, names(SummarizedExperiment::colData(sce)))
  if (!is.null(random_effect)) {
    match.arg(random_effect, names(SummarizedExperiment::colData(sce)))
  }
  features <-
    match.arg(features,
              SummarizedExperiment::rowData(sce)$marker_name,
              several.ok = TRUE)
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
  if (!is.null(random_effect)) {
    DF$G <- sce[[random_effect]]
  }
  
  # if (weighted)
  #   my_weights <- NULL
  # else
  #   my_weights <- rep(1 / CATALYST::ei(sce)$n_cells, CATALYST::ei(sce)$n_cells)
  
  data.table::rbindlist(BiocParallel::bplapply(features, function(marker, DF) {
    message(paste('fitting marker', marker))
    if (is.null(random_effect))
      my_formula <- as.formula(paste0(marker, " ~ Y"))
    else
      my_formula <- as.formula(paste0(marker, " ~ Y + (1|G)"))
    bereg <-
      glmmTMB::glmmTMB(
        formula = my_formula,
        ziformula = ~ .,
        data = DF,
        family = glmmTMB::beta_family()
        # weights = my_weights
      )
    out = summary(bereg)
    # res[i, 2:3] <- c(out$coefficients$cond[2, 4],
    #                  out$coefficients$zi[2, 4])#c(out$coef_table["Ytreatment", "p-value"], out$coef_table_zi["Ytreatment", "p-value"])
    data.table::data.table(marker_id = marker, 
                            p_val = out$coefficients$cond[2, 4],
                           p_val_zi = out$coefficients$zi[2, 4])
  }, DF, BPPARAM = bppar))
}

# library(CATALYST)
# library(glmmTMB)
# 
# sapply(list.files("functions", full.names = TRUE),
#        source,
#        environment())()
# 
# sce <- simulateSCE()
# res <- hurdleBeta(sce, "condition", "patient_id") # , "patient_id")
