library(CATALYST)
library(gamlss.dist)
library(data.table)

data(PBMC_panel, PBMC_md, PBMC_fs)
sce_pbmc<- prepData(PBMC_fs, PBMC_panel, PBMC_md, transform=T)

# zero inflated distribution
nn <- ncol(sce_pbmc) #- 1000
set.seed(1)
bdata <- rbeta(nn, shape1 = 1, shape2 = 5)
#bdata <- c(bdata, as.vector(rep(0,1000)))



# distribution of marker CD3
CD3 <- assays(sce_pbmc)$exprs["CD3",]
CD3 <- (CD3 - min(CD3))/(max(CD3)-min(CD3))

dist <- data.table(beta = bdata, CD3 = CD3)
dist <- melt(dist, measure.vars = c("beta", "CD3"), variable.name = "distribution", value.name = "value")

ggplot(dist, aes(x=value, fill=distribution)) + geom_density(alpha=0.4)

# ZIBSeq
exprs <- as.data.frame(assays(sce_pbmc)$exprs)
group <-colData(sce_pbmc)$condition
exprs <- t(exprs)

result <- zibSeq(data = exprs, outcome = group, transform=F)
padj <- p.adjust(result$pvalues, method="BH")
names(padj) <- colnames(exprs)

zibSeq <- function (data, outcome, transform = F, alpha = 0.05) 
{
  X = data
  Y = outcome
  P = dim(X)[2]
  beta = matrix(data = NA, P, 2)
  for (i in 1:P) {
    x.prop = X[, i]
    x.prop <- (x.prop - min(x.prop))/(max(x.prop)-min(x.prop))
    x.prop[which(x.prop==1)] <- x.prop[which(x.prop==1)] - 2.225074e-10
    if (transform == T) {
      x.prop = sqrt(x.prop)
    }
    bereg = gamlss(x.prop ~ Y, family = BEZI(sigma.link = "identity"), 
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





# gamlss
exprs <- as.data.frame(assays(sce_pbmc)$exprs)
group <-colData(sce_pbmc)$condition

num_markers <- nrow(exprs)

for (i in 1:num_markers){
  exprs_i <- as.vector(exprs[i,])
  zib <- gamlss(exprs_i ~ group, family = "BEZI")

  zinb_try <- try(gamlssML(c(exprs_1, exprs_2), 
                           family = "BEZI"), silent = TRUE)
}

?gamlssML
gamlssML(, family = "BEZI")
