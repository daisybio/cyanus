hurdleBeta <- function(sce,
                       condition,
                       random_effect = NULL,
                       features = SummarizedExperiment::rowData(sce)$marker_name,
                       assay_to_use = "exprs",
                       weighted = TRUE) {
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
  X <- X[, features]
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
  res <- matrix(data = NA, length(features), 3)
  
  if (weighted)
    my_weights <- NULL
  else
    my_weights <- rep(1 / ei(sce)$n_cells, ei(sce)$n_cells)
  
  for (i in seq(length(features))) {
    if (is.null(random_effect))
      my_formula <- as.formula(paste0(names(DF)[i], " ~ Y"))
    else
      my_formula <- as.formula(paste0(names(DF)[i], " ~ Y + (1|G)"))
    bereg <-
      glmmTMB::glmmTMB(
        formula = my_formula,
        ziformula = ~ .,
        data = DF,
        family = beta_family(),
        weights = my_weights
      )
    out = summary(bereg)
    res[i, 2:3] <- c(out$coefficients$cond[2, 4],
                     out$coefficients$zi[2, 4])#c(out$coef_table["Ytreatment", "p-value"], out$coef_table_zi["Ytreatment", "p-value"])
  }
  res[, 1] <- features
  colnames(res) <- c("marker", "p-value", "p-value-zi")
  res
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
