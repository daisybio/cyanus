library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)
library(ggvenn)

############### Preprocessing ###############

makeSampleSelection <- function(sce=sce, deselected_samples){
  # All Samples
  all_samples <- as.character(unique(colData(sce)$sample_id))
  
  # Samples to keep
  samples <- all_samples[!all_samples %in% deselected_samples]
  
  # Filter Single Cell Experiment
  sce <- filterSCE(sce, sample_id %in% samples)
  
  return (sce)
}

makePatientSelection <- function(sce, deselected_patients){
  # All patients
  all_patients <- as.character(unique(colData(sce)$patient_id))
  
  # Patients to keep
  patients <- all_patients[!all_patients %in% deselected_patients]
  
  # Filter Single Cell Experiment
  sce <- filterSCE(sce, patient_id %in% patients)

  return (sce)
}

############### Visualization ###############

#run the dimensionality reduction; function either calls CATALYST::runDR or runIsomap
runDimRed <- function(sce, dr_chosen = c("UMAP", "TSNE", "PCA", "MDS", "DiffusionMap", "Isomap"), 
                      cells_chosen = NULL, feature_chosen = "type", assay_chosen = "exprs", scale = T, k = 3){
  match.arg(dr_chosen)
  if (dr_chosen == "Isomap") {
    sce <- runIsomap(
      sce,
      cells = cells_chosen,
      features = feature_chosen,
      assay = assay_chosen,
      scale = scale,
      k = k
    )
    
  } else{
    sce <- CATALYST::runDR(
      sce,
      dr = dr_chosen,
      cells = cells_chosen,
      features = feature_chosen,
      assay = assay_chosen,
      scale = scale
    )
  }
}

#plot the dimensionality reduction by calling:
#CATALYST::plotDR(sce, dr = dr_chosen, color_by = color_chosen, facet_by = facet_chosen, assay = assay_chosen, scale = scale)

#adapted from CATALYST::runDR
runIsomap <- function (x, cells = NULL, features = "type", assay = "exprs", scale = TRUE, k = 5) 
{
  stopifnot(is(x, "SingleCellExperiment"))
  #sample cells from each sample if cells is specified, if not take the whole SCE
  if (is.null(cells)) {
    cs <- TRUE
  }
  else {
    if (is.null(x$sample_id)) 
      stop("colData column sample_id not found,\n ", " but is required to downsample cells.")
    
    stopifnot(is.numeric(cells), length(cells) == 1, as.integer(cells) == 
                cells, cells > 0)
    
    cs <- split(seq_len(ncol(x)), x$sample_id)
    cs <- unlist(lapply(cs, function(u) sample(u, min(cells, length(u)))))
  }
  #subset SCE according to selected markers
  if(is.null(features)){
    selectedMarkers <- TRUE
  }else if(features %in% c("type", "state")){
    selectedMarkers <- rowData(x)[rowData(x)$marker_class == features, "marker_name"]
  }else if(features %in% rownames(x)){
    selectedMarkers <- rowData(x)[rowData(x)$marker_name %in% features, "marker_name"]
  }
  
  selectedCounts <- assays(x[selectedMarkers, ])[[assay]]
  selectedCounts <- selectedCounts[,  cs]
  if(scale){
    selectedCounts <- scale(selectedCounts, center = TRUE, scale = TRUE)
  }
  res_isomap <- embed(t(selectedCounts), "Isomap", knn = k)
  
  if (is.null(cells)){
    reducedDim(x, "Isomap") <- res_isomapb@data@data
    return(x)
  }
  
  m <- matrix(NA, nrow = ncol(x), ncol = ncol(res_isomap@data@data))
  m[cs, ] <- res_isomap@data@data
  reducedDim(x, "Isomap") <- m
  return(x)
}

############### Clustering ###############

source("../server/clusterFun.R")

############### Differential Expression ###############

prepDiffExp <- function(sce, contrastVars, colsDesign, colsFixed, colsRandom,
                        method = c( "diffcyt-DA-edgeR", "diffcyt-DA-voom","diffcyt-DA-GLMM", "diffcyt-DS-limma", "diffcyt-DS-LMM") ){
  match.arg(method)
  parameters <- list()
  parameters[["ei"]] <- metadata(sce)$experiment_info
  
  if(method %in% c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DS-limma")){
    
    parameters[["design"]] <- diffcyt::createDesignMatrix(parameters[["ei"]], cols_design = colsDesign)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, parameters[["design"]], designMatrix = T)
    
  }else if(method == "diffcyt-DS-LMM"){
    parameters[["formula"]] <- diffcyt::createFormula(parameters[["ei"]], cols_fixed = colsFixed, cols_random = colsRandom)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, diffcyt::createDesignMatrix(parameters[["ei"]], cols_design = colsFixed), designMatrix = T)
  }else{
    
    parameters[["formula"]] <- diffcyt::createFormula(parameters[["ei"]], cols_fixed = colsFixed, cols_random = colsRandom)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, colsFixed, designMatrix = F)
    
  }
  return(parameters) 
}

#Then run: 
#diffcyt::diffcyt(d_input = sce,
#                 design = parameters[["design"]],
#                 formula = parameters[["formula"]],
#                 contrast = parameters[["contrast"]],
#                 analysis_type = c("DA", "DS"),
#                 method_DS = c("diffcyt-DS-limma", "diffcyt-DS-LMM"),
#                 method_DA = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"),
#                 clustering_to_use = "all",
#                 markers_to_test = markersToTest, 
#                 trend = trend, 
#                 normalize = T, ...)

#Then run: CATALYST::plotDiffHeatmap(...) and CATALYST::plotPbExprs(...)

createCustomContrastMatrix <- function(sce, contrastVars, matrix, designMatrix = T){
  if(designMatrix){
    #the entries have to correspond to the columns of the design matrix
    cnames <- colnames(matrix)
    bool <- getBools(cnames, contrastVars)
    contrast <- createContrast(bool)
    print(contrast)
    return(contrast)
  }else{
    #the entries have to correspond to the levels of the fixed effect terms in the model formula
    lvlList <- lapply(matrix, function(x){levels(colData(sce)[[x]])})
    names(lvlList) <- matrix
    bool <- getBools(matrix, contrastVars)
    names(bool) <- matrix
    contrast <- unlist(lapply(names(lvlList), function(x){
      return( c( 0, rep(bool[x], length(lvlList[[x]]) -1 )) ) 
    }))
    print(contrast)
    return(createContrast(unname(contrast)))
  }
}

getBools <- function(names, contrastVars){
  bool <- unlist(lapply(names, function(x){
    any(lapply(contrastVars, function(y){
      grepl(y,x, fixed = T )
    }))
  }))
  bool <- as.numeric(bool)
  return(bool)
}

SigEMD <- function(sce, k, condition, Hur_gene=NULL, binSize=NULL, nperm=100, assay="exprs", seed=1, parallel=FALSE, permute_samples=FALSE) {
  library(aod)
  library(arm)
  library(fdrtool)
  library(lars)
  library(emdist)
  library(data.table)
  source("../SigEMD/FunImpute.R")
  source("../SigEMD/SigEMDHur.R")
  source("../SigEMD/SigEMDnonHur.R")
  source("../SigEMD/plot_sig.R")
  
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


# SigEMD_perm_per_sample <- function(sce) {
#   
#   res_real <- SigEMD()
#   res_real <- res_real$all$emdall[, .(real_emd = emd, marker_id)]
#   setkey(res_real, marker_id)
#   res_real
#   
#   sce_dual_sampled <- sce_dual
#   
#   nperm <- 100
#   used_permutations <- list()
#   used_permutations[[nperm + 1]] <- metadata(sce_dual_sampled)$experiment_info$activated_baseline
#   all_results <- list()
#   
#   for (i in 1:nperm){
#     ei <- metadata(sce_dual_sampled)$experiment_info
#     
#     # sample A and B
#     set.seed(i)
#     repeat {
#       condition <- sample(ei$activated_baseline)
#       if (Position(function(x) identical(x, condition), used_permutations, nomatch = 0) == 0) {
#         used_permutations[[i]] <- condition
#         break
#       }
#     }
#     
#     # edit activated_baseline condition in experiment_info
#     ei$activated_baseline <- condition
#     metadata(sce_dual_sampled)$experiment_info <- ei
#     
#     # edit activated_baseline in colData (not corrected in sample name)
#     coldata_condition <- rep(condition, times=ei$n_cells)
#     colData(sce_dual_sampled)$activated_baseline <- as.factor(coldata_condition)
#     
#     all_results[[i]] <- SigEMD(sce = sce_dual_sampled,
#                                k = "all",
#                                condition = "activated_baseline",
#                                Hur_gene=rownames(sce_dual_sampled),
#                                nperm = 1,
#                                parallel = TRUE)
#     
#   }
#   all_perms <- rbindlist(lapply(all_results, function(x) x$all$emd), idcol = "permutation")
#   all_perms[, p_val := NULL]
#   all_perms[, p_adj := NULL]
#   setkey(all_perms, marker_id)
#   
#   res_agg <- all_perms[res_dual_real][, .(p_val = (sum(emd >= real_emd) + 1)/(nperm + 1)), by = marker_id]
#   res_agg[, p_adj := p.adjust(p_val, "BH")]
#   res_agg[order(p_adj)]
# }


DEsingleSCE <- function(sce, condition, k, assay="exprs", parallel=FALSE){
  library(DEsingle)
  
  CATALYST:::.check_sce(sce, TRUE)
  k <- CATALYST:::.check_k(sce, k)
  CATALYST:::.check_cd_factor(sce, condition)
  assay <- match.arg(assay, names(SummarizedExperiment::assays(sce)))
  
  sce_desingle <- sce
  new_counts <- assay(sce_desingle, "exprs")
  new_counts <- (abs(new_counts) + new_counts)/2
  assay(sce_desingle, "counts") <- new_counts
  set.seed(1)
  
  cluster_ids <- cluster_ids(sce_desingle, k)
  res <- lapply(levels(cluster_ids), function(cluster_id) {
    print(sprintf("calculating DEsingle for cluster %s", cluster_id))
    
    group <- colData(sce_desingle[, cluster_ids == cluster_id])[[condition]]

    
    results <- DEsingle::DEsingle(sce_desingle[, cluster_ids == cluster_id], group, parallel = parallel)
    results.classified <- DEsingle::DEtype(results = results, threshold = 0.05)
    
    data.table::setnames(results.classified, old = c("pvalue", "pvalue.adj.FDR"), new = c("p_val", "p_adj"))
    results.classified$cluster_id <- cluster_id
    results.classified$marker_id <- rownames(results.classified)
    
    results.classified
  })
  
  return(data.table::rbindlist(res))
}

runDS <- function(sce, condition, de_methods = c("limma", "LMM", "SigEMD", "DEsingle"), k = "all", parallel = TRUE, features = c("all", "type","state"), ...) {
  de_methods <- match.arg(de_methods, several.ok = TRUE)
  features <- match.arg(features, several.ok = FALSE)
  
  # if (is.null(marker)) marker <- rownames(sce)
  # else marker <- match.arg(marker, rownames(sce), several.ok = TRUE)
  extra_args <- list(...)
  
  result <- list()
  if ("limma" %in% de_methods) {
    message("Using limma")
    parameters <- prepDiffExp(sce, contrastVars = c(condition), colsDesign = c(condition), method = "diffcyt-DS-limma")
    
    #blockID <- metadata(sce_dual_ab)$experiment_info[["patient_id"]]
    
    markers_to_test <- getMarkersToTest(sce,"limma",features)
    limma_res <- diffcyt::diffcyt(d_input = sce,
                                      design = parameters[["design"]],
                                      contrast = parameters[["contrast"]],
                                      analysis_type = "DS",
                                      method_DS = "diffcyt-DS-limma",
                                      clustering_to_use = k,
                                      markers_to_test = markers_to_test)
  
    result$limma <- limma_res

  }
  if ("LMM" %in% de_methods) {
    message("Using LMM")
    parameters <-  prepDiffExp(sce, contrastVars = c(condition), colsFixed = c(condition), colsRandom=c("patient_id"), method = "diffcyt-DS-LMM")
    
    markers_to_test <- getMarkersToTest(sce,"LMM",features)
    LMM_res <- diffcyt::diffcyt(d_input = sce,
                                    formula = parameters[["formula"]],
                                    contrast = parameters[["contrast"]],
                                    analysis_type = "DS",
                                    method_DS = "diffcyt-DS-LMM",
                                    clustering_to_use = k,
                                    markers_to_test = markers_to_test)
    
    result$LMM <- LMM_res
  }
  if ("SigEMD" %in% de_methods) {
    message("Using SigEMD")
    nperm <- ifelse(is.null(extra_args$nperm), 100, extra_args$nperm)
    permute_samples <- ifelse(is.null(extra_args$permute_samples), FALSE, extra_args$permute_samples)
    
    markers_to_test <- getMarkersToTest(sce,"SigEMD",features)
    
    # make subselection
    sce_SigEMD <- filterSCE(sce, marker_name %in% markers_to_test)
    
    SigEMD_res <- SigEMD(sce_SigEMD, k, condition, Hur_gene=rownames(sce_SigEMD), nperm=nperm, parallel=parallel, permute_samples=permute_samples)
    SigEMD_res$overall <- data.table::rbindlist(lapply(SigEMD_res, function(x) x$emdall))
    
    result$SigEMD <- SigEMD_res
  }
  if ("DEsingle" %in% de_methods) {
    message("Using DEsingle")
    
    markers_to_test <- getMarkersToTest(sce,"DEsingle",features)
    
    # make subselection
    sce_DEsingle <- filterSCE(sce, marker_name %in% markers_to_test)
    
    DEsingle_res <- DEsingleSCE(sce_DEsingle, condition, k, parallel=parallel)
    
    
    result$DEsingle <- DEsingle_res
  }
  result
}

# create venn diagram of all methods that were performed in runDS
createVennDiagram <- function(res) {
  input_venn <- list()
  
  # take of each method the data table containing the pvalues
  for (ds_method in names(res)) {
    if (ds_method == "DEsingle") {
      output <- res[[ds_method]]
      result <- data.frame(output)[c("marker_id", "p_val", "p_adj")]
    }
    if (ds_method == "SigEMD") {
      output <- res[[ds_method]]$overall
      result <- data.frame(output)[c("marker_id", "p_val", "p_adj")]
    }
    if (ds_method %in% c("LMM", "limma")) {
      result <-
        data.frame(rowData(res[[ds_method]]$res))[c("marker_id", "p_val", "p_adj")]
    }
    result$significant <- result$p_adj < 0.05
    significants <-
      unlist(subset(
        result,
        significant == TRUE,
        select = c(marker_id),
        use.ames = FALSE
      ))
    input_venn[[ds_method]] <- significants
  }
  
  venn <- ggvenn::ggvenn(input_venn, show_elements = TRUE, label_sep ="\n", fill_alpha = 0.3, set_name_size = 8, text_size = 4)
  return(venn)
}

# get appropriate vector for each method containing the markers that want to be tested
getMarkersToTest <- function(sce, ds_method, features){
  if (ds_method %in% c("limma", "LMM")){
    # vector for diffcyt methods to test on specific markers
    is_marker <- rowData(sce)$marker_class %in% c("type","state")
    if (features == "all"){
      markers_to_test <- (rowData(sce)$marker_class %in% c("type","state"))[is_marker]
    } else if (features == "state"){
      markers_to_test <- (rowData(sce)$marker_class == "state")[is_marker]
    } else if (features == "type") {
      markers_to_test <- (rowData(sce)$marker_class == "type")[is_marker]
    }
  } else if (ds_method %in% c("SigEMD", "DEsingle")){
    if (features == "all"){
      markers_to_test <- rowData(sce)[rowData(sce)$marker_class %in% c("type","state"),]$marker_name
    } else if (features == "state"){
      markers_to_test <- rowData(sce)[rowData(sce)$marker_class=="state",]$marker_name
    } else if (features == "type"){
      markers_to_test <- rowData(sce)[rowData(sce)$marker_class=="type",]$marker_name
    }
  }
  return (markers_to_test)
}

