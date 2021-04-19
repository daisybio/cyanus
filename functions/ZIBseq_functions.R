zibSeq <- function (sce, condition, alpha = 0.05)
{
  X = t(as.data.frame(assays(sce)$exprs))
  Y = colData(sce)[[condition]]
  P = dim(X)[2]
  beta = matrix(data = NA, P, 2)
  for (i in 1:P) {
    x.prop <- X[, i]
    x.prop[x.prop < 0] <- 0
    x.prop <- (x.prop - min(x.prop))/(max(x.prop)-min(x.prop))
    x.prop[which(x.prop==1)] <- x.prop[which(x.prop==1)] - 2.225074e-10
    bereg = gamlss(x.prop ~ Y, family = BEZI(), 
                   trace = FALSE, control = gamlss.control(n.cyc = 100))
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
