zibSeq <- function (sce, condition, random_effect=NULL, alpha = 0.05)
{
  X = t(as.data.frame(assays(sce)$exprs))
  Y = colData(sce)[[condition]]
  
  if (!is.null(random_effect)){
    R = colData(sce)[[random_effect]]
    message("Fitting a zero-inflated beta mixed model with random effects...")
  } else {
    message("Fitting a zero-inflated beta model...")
  }
  P = dim(X)[2]
  beta = matrix(data = NA, P, 2)
  for (i in 1:P) {
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
                             trace = FALSE, control = gamlss.control(n.cyc = 100), data = data)
    } else {
      # without random effect
      bereg = gamlss::gamlss(exprs ~ condition, family = BEZI(), 
                             trace = FALSE, control = gamlss.control(n.cyc = 100), data = data)
    }
  
    out = summary(bereg)
    beta[i, ] = out[2, c(1, 4)]
  }
  pvalues = beta[, 2]
  qvalues = ZIBseq:::calc_qvalues(pvalues)
  sig = which(qvalues < alpha)
  sigFeature = colnames(X)[sig]
  list(sigFeature = sigFeature, useFeature = P, qvalues = qvalues, 
       pvalues = pvalues)
}