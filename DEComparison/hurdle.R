library(CATALYST)
library(glmmTMB)

sapply(list.files("functions", full.names = TRUE),
       source,
       environment())()

hurdleBeta <- function(sce,
                       condition,
                       group = NULL,
                       features = SummarizedExperiment::rowData(sce)$marker_name,
                       assay_to_use = "exprs") {
  match.arg(assay_to_use, SummarizedExperiment::assayNames(sce))
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
  DF$G <- sce[[group]]
  res <- matrix(data = NA, length(features), 3)
  for (i in seq(length(features))) {
    # bereg <-
    #   GLMMadaptive::mixed_model(
    #     fixed = as.formula(paste0(names(DF)[i], " ~ Y")),
    #     random = ~ 1 | G,
    #     data = DF,
    #     family = hurdle.beta.fam(),
    #     zi_fixed = ~ Y
    #   )
    bereg <-
      glmmTMB::glmmTMB(
        formula = as.formula(paste0(names(DF)[i], " ~ Y + (1|G)")),
        ziformula = ~ .,
        data = DF,
        family = beta_family()
      )
    #browser()
    #bereg2 <- update(bereg, zi_random = ~ 1 | G)
    #anova(bereg, bereg2)
    #bereg_nohurdle <- GLMMadaptive::mixed_model(fixed = as.formula(paste0(names(DF)[i]," ~ Y")), random = ~ 1 | G, data = DF,
    #                                            family = beta.fam(), n_phis = 1)
    #browser()
    out = summary(bereg)
    res[i, 2:3] <- c(out$coefficients$cond[2, 4],
                     out$coefficients$zi[2, 4])#c(out$coef_table["Ytreatment", "p-value"], out$coef_table_zi["Ytreatment", "p-value"])
  }
  res[, 1] <- features
  colnames(res) <- c("marker", "p-value", "p-value-zi")
  res
}

sce <- simulateSCE()
#res <- hurdleBeta(sce, "condition", "patient_id")

X <- t(SummarizedExperiment::assay(sce, "exprs"))
#X <- X[, features]
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
DF$Y <- sce$condition
DF$G <- sce$patient_id
glmmTMB_res <-
  glmmTMB::glmmTMB(
    formula = m01 ~ Y + (1 |
                           G),
    ziformula = ~ .,
    data = DF,
    family = beta_family()
  )
bereg <- gamlss::gamlss(x.prop ~ Y, family = BEZI(), 
                       trace = FALSE, control = gamlss.control(n.cyc = 100))
