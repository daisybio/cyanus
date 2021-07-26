Rcpp::cppFunction('double emdC(NumericVector a, NumericVector b) {
  int n = a.size();
  NumericVector dist = NumericVector(n);
  double emd = 0;
  for(int i = 0; i < (n - 1); ++i) {
    dist[i + 1] = a[i] - b[i] + dist[i];
  }
  dist = abs(dist);
  for (auto& d : dist)
    emd += d;
  return emd;
}')

myEMD <-  function(A, B, binSize = NULL) {
  stopifnot(is.numeric(A) & is.numeric(B))#  & (length(A) == length(B))) they do not have to be the same length, but will it have a great effect for large numbers?
  if (is.null(binSize)) binSize <- 2 * IQR(c(A[A!=0], B[B!=0])) / length(c(A[A!=0], B[B!=0]))^(1/3)
  
  bins <- seq(floor(min(c(A, B))),
              ceiling(max(c(A, B))),
              by=binSize )
  if (max(bins) < max(A,B)) bins <- c(bins, bins[length(bins)] + binSize)
  
  histA <- hist(A, breaks=bins, plot=FALSE)
  histB <- hist(B, breaks=bins, plot=FALSE)
  
  densA <- histA$density
  densA <- densA/sum(densA)
  densB <- histB$density
  densB <- densB/sum(densB)
  
  emdC(densA, densB)
}

rowwiseEMD <- function(mat, condition, binSize = NULL) {
  stopifnot(is.matrix(mat), is.numeric(mat), nlevels(as.factor(condition)) == 2, ncol(mat)==length(condition))
  
  condition <- as.factor(condition)
  
  result <- apply(mat, 1, function(marker) {
    grouped <- split(marker, condition)
    myEMD(grouped[[1]], grouped[[2]])
  })
  out_dt <- data.table::as.data.table(result, keep.rownames="marker_id")
  out_dt[, marker_id := as.factor(marker_id)]
  out_dt
}

cyEMD <- function(sce, condition, binSize=NULL, nperm=100, assay="exprs", seed=1, parallel=FALSE) {
  # suppressPackageStartupMessages(library(data.table))
  bppar <- BiocParallel::bpparam()
  
  if (!parallel)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  
  set.seed(1)
  assay <- match.arg(assay, names(SummarizedExperiment::assays(sce)))
  
  # get data matrix
  data <- SummarizedExperiment::assay(sce, assay)
  
  # compute actual emd
  emd_real <- rowwiseEMD(mat = data, condition = sce[[condition]], binSize = binSize)
  data.table::setnames(emd_real, "result", "real_emd")
  data.table::setkey(emd_real, marker_id)
  
  # compute permutations of sample conditions
  sceEI <- CATALYST::ei(sce)
  perms <- RcppAlgos::permuteSample(sceEI[[condition]], n = nperm, seed = seed)
  perm_res <- BiocParallel::bplapply(as.data.frame(t(unclass(perms))), function(perm, sceEI, data, binSize) {
    condition_permutation_cells <- rep(perm, times=sceEI$n_cells)
    rowwiseEMD(mat = data, condition = condition_permutation_cells, binSize = binSize)
  }, sceEI, data, binSize, BPPARAM = bppar)
  
  # gather results and compute empirical p-value
  all_perms <- data.table::rbindlist(perm_res, idcol = "permutation")
  data.table::setkey(all_perms, marker_id)
  res_agg <- all_perms[emd_real][, .(p_val = (sum(result >= real_emd) + 1)/(nperm + 1)), by = c("marker_id", "real_emd")]
  data.table::setnames(res_agg, "real_emd", "emd")
  # res_agg[, p_adj := p.adjust(p_val, "BH")]
  
  return(res_agg)
}