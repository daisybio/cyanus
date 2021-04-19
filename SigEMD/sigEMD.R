
SigEMD <- function(sce, k, condition, Hur_gene=NULL, binSize=NULL, nperm=100, assay="exprs", seed=1, parallel=FALSE, permute_samples=FALSE) {
  library(aod)
  library(arm)
  library(fdrtool)
  library(lars)
  library(emdist)
  library(data.table)
  source("SigEMD/FunImpute.R")
  source("SigEMD/SigEMDHur.R")
  source("SigEMD/SigEMDnonHur.R")
  source("SigEMD/plot_sig.R")
  
  set.seed(1)
  
  CATALYST:::.check_sce(sce, TRUE)
  k <- CATALYST:::.check_k(sce, k)
  CATALYST:::.check_cd_factor(sce, condition)
  assay <- match.arg(assay, names(SummarizedExperiment::assays(sce)))
  
  
  cluster_ids <- cluster_ids(sce, k)
  res <- lapply(levels(cluster_ids), function(curr_cluster_id) {
    print(sprintf("calculating SigEMD for cluster %s", curr_cluster_id))
    
    sce_cluster <- filterSCE(sce, cluster_id == curr_cluster_id, k = k)
    data <- assay(sce_cluster, assay)
    
    
    data <- dataclean((abs(data)+data)/2)
    colnames(data) <- as.character(seq.int(to = ncol(data)))
    
    condition_cluster <- colData(sce_cluster)[[condition]]
    names(condition_cluster) <- colnames(data)
    results <- calculate_single(data =  data,condition =  condition_cluster,Hur_gene = Hur_gene, binSize, nperm=ifelse(permute_samples, 1, nperm), parallel = parallel)
    
    results$emdall <- as.data.frame(results$emdall)
    data.table::setnames(results$emdall, old = c("pvalue", "padjust"), new = c("p_val", "p_adj"))
    results$emdall$cluster_id <- curr_cluster_id
    results$emdall$marker_id <- rownames(results$emdall)
    
    if (permute_samples) {
      # gather real result
      res_real <- as.data.table(results$emdall)[, .(real_emd = emd, marker_id, cluster_id)]
      setkey(res_real, marker_id)
      
      # permute samplewise
      used_permutations <- list()
      used_permutations[[nperm + 1]] <- metadata(sce_cluster)$experiment_info[[condition]]
      all_results <- list()
      ei <- metadata(sce_cluster)$experiment_info
      
      for (i in 1:nperm){
        message(sprintf("Permutation number: %d", i))
        
        # sample condition
        set.seed(i)
        repeat {
          condition_permutation <- sample(ei[[condition]])
          if (Position(function(x) identical(x, condition_permutation), used_permutations, nomatch = 0) == 0) {
            used_permutations[[i]] <- condition_permutation
            break
          }
        }
        
        condition_permutation_cells <- rep(condition_permutation, times=ei$n_cells)
        names(condition_permutation_cells) <- colnames(data)
        
        all_results[[i]] <- as.data.frame(calculate_single(data =  data,condition =  condition_permutation_cells,Hur_gene = Hur_gene, binSize, nperm=1, parallel = parallel)$emdall)
        
      }
      all_perms <- rbindlist(lapply(all_results, as.data.table, keep.rownames = "marker_id"), idcol = "permutation")
      all_perms[, pvalue := NULL]
      all_perms[, padjust := NULL]
      setkey(all_perms, marker_id)
      
      res_agg <- all_perms[res_real][, .(p_val = (sum(emd >= real_emd) + 1)/(nperm + 1)), by = c("marker_id", "real_emd", "cluster_id")]
      setnames(res_agg, "real_emd", "emd")
      res_agg[, p_adj := p.adjust(p_val, "BH")]
      results$emdall <- res_agg
      
    }
    
    
    
    results
  })
  
  names(res) <- levels(cluster_ids)
  
  return(res)
}

