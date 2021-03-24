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
runIsomap <- function (x, cells = NULL, features = "type", assay = "exprs", scale = TRUE, k = 5, dimensions = 2) 
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
  res_isomap <- embed(t(selectedCounts), "Isomap", knn = k, ndim = dimensions)
  
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

#source("../server/clusterFun.R")


############### Differential Expression ###############

prepDiffExp <- function(sce, contrastVars, colsDesign, colsFixed, colsRandom=NULL,
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

library(Rcpp)
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
  stopifnot(is.matrix(mat) & is.numeric(mat) & (nlevels(as.factor(condition)) == 2) & (ncol(mat)==length(condition)))
  
  condition <- as.factor(condition)
  
  result <- apply(mat, 1, function(marker) {
    grouped <- split(marker, condition)
    myEMD(grouped[[1]], grouped[[2]])
  })
  data.table::as.data.table(result, keep.rownames="marker")
}

sceEMD <- function(sce, k, condition, binSize=NULL, nperm=100, assay="exprs", seed=1, parallel=FALSE) {
  library(data.table)
  
  bppar <- BiocParallel::bpparam()
  
  if (parallel == FALSE)
    bppar <- BiocParallel::SerialParam(progressbar = TRUE)
  
  set.seed(1)
  
  CATALYST:::.check_sce(sce, TRUE)
  k <- CATALYST:::.check_k(sce, k)
  CATALYST:::.check_cd_factor(sce, condition)
  assay <- match.arg(assay, names(SummarizedExperiment::assays(sce)))
  
  
  cluster_ids <- cluster_ids(sce, k)
  res <- lapply(levels(cluster_ids), function(curr_cluster_id) {
    print(sprintf("calculating sceEMD for cluster %s", curr_cluster_id))
    
    sce_cluster <- filterSCE(sce, cluster_id == curr_cluster_id, k = k)
    data <- assay(sce_cluster, assay)
    
    condition_cluster <- colData(sce_cluster)[[condition]]
    emd_real <- rowwiseEMD(mat = data, condition = condition_cluster, binSize = binSize)
    emd_real$cluster_id <- curr_cluster_id
    setnames(emd_real, "result", "real_emd")
    setkey(emd_real, marker)
    
    ei <- metadata(sce_cluster)$experiment_info
    perms <- RcppAlgos::permuteSample(ei[[condition]], n = nperm, seed = seed)
    
    perm_res <- BiocParallel::bplapply(as.data.frame(t(unclass(perms))), function(perm, ei, data, binSize) {
      condition_permutation_cells <- rep(perm, times=ei$n_cells)
      rowwiseEMD(mat = data, condition = condition_permutation_cells, binSize = binSize)
    }, ei, data, binSize, BPPARAM = bppar)
    
    all_perms <- rbindlist(perm_res, idcol = "permutation")
    setkey(all_perms, marker)
    
    res_agg <- all_perms[emd_real][, .(p_val = (sum(result >= real_emd) + 1)/(nperm + 1)), by = c("marker", "real_emd", "cluster_id")]
    setnames(res_agg, "real_emd", "emd")
    res_agg[, p_adj := p.adjust(p_val, "BH")]
    res_agg
  })
  
  return(data.table::rbindlist(res))
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

runDS <- function(sce, condition, random_effect = NULL, de_methods = c("limma", "LMM", "SigEMD", "DEsingle"), k = "all", parallel = TRUE, features = c("all", "type","state"), ...) {
  de_methods <- match.arg(de_methods, several.ok = TRUE)
  features <- match.arg(features, several.ok = FALSE)
  
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
    parameters <-  prepDiffExp(sce, contrastVars = c(condition), colsFixed = c(condition), colsRandom=c(random_effect), method = "diffcyt-DS-LMM")
    
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


############# CATALYST reprogrammed functions #############


# function anno_features from CATALYST with different color palette
.anno_factors <- function (x, ids, which, type = c("row", "column")) 
{
  type <- match.arg(type)
  cd <- colData(x)
  df <- data.frame(cd, check.names = FALSE)
  df <- select_if(df, ~!is.numeric(.))
  df <- mutate_all(df, ~droplevels(factor(.x)))
  m <- match(ids, df$sample_id)
  ns <- split(df, df$sample_id) %>% lapply(mutate_all, droplevels) %>% 
    lapply(summarize_all, nlevels) %>% do.call(what = "rbind")
  keep <- names(which(colMeans(ns) == 1))
  keep <- setdiff(keep, c("sample_id", "cluster_id"))
  if (is.character(which)) 
    keep <- intersect(keep, which)
  if (length(keep) == 0) 
    return(NULL)
  df <- df[m, keep, drop = FALSE]
  lvls <- lapply(as.list(df), levels)
  nlvls <- vapply(lvls, length, numeric(1))
  #pal <- brewer.pal(8, "Set2")
  pal <-  c("#999999","#009E73","#E69F00", "#56B4E9", "#F0E442","#D55E00","#0072B2", "#CC79A7")
  names(is) <- is <- colnames(df)
  cols <- lapply(is, function(i) {
    if (nlvls[i] > length(pal)) 
      pal_i <- colorRampPalette(pal)(max(nlvls))
    else pal_i <- pal
    u <- pal_i[seq_len(nlvls[i])]
    names(u) <- lvls[[i]]
    u
  })
  ComplexHeatmap::HeatmapAnnotation(which = type, df = df, col = cols, gp = gpar(col = "white"))
}

plotExprHeatmap <- function (x, features = NULL, by = c("sample_id", "cluster_id", 
                                                        "both"), k = "meta20", m = NULL, assay = "exprs", fun = c("median", 
                                                                                                                  "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, 
                             row_anno = TRUE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
                             row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, 
                             bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
                             k_pal = CATALYST:::.cluster_cols, m_pal = k_pal, distance = c("euclidean", 
                                                                                           "maximum", "manhattan", "canberra", "binary", "minkowski"), 
                             linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
                                         "median", "centroid", "ward.D2")) 
{
  args <- as.list(environment())
  CATALYST:::.check_args_plotExprHeatmap(args)
  distance <- match.arg(distance)
  linkage <- match.arg(linkage)
  scale <- match.arg(scale)
  fun <- match.arg(fun)
  by <- match.arg(by)
  x <- x[unique(CATALYST:::.get_features(x, features)), ]
  if (by != "sample_id") {
    CATALYST:::.check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
  }
  if (by == "both") 
    by <- c("cluster_id", "sample_id")
  .do_agg <- function() {
    z <- CATALYST:::.agg(x, by, fun, assay)
    if (length(by) == 1) 
      return(z)
    set_rownames(do.call("rbind", z), levels(x$cluster_id))
  }
  .do_scale <- function() {
    if (scale == "first") {
      z <- assay(x, assay)
      z <- CATALYST:::.scale_exprs(z, 1, q)
      assay(x, assay, FALSE) <- z
      return(x)
    }
    else CATALYST:::.scale_exprs(z, 1, q)
  }
  z <- switch(scale, first = {
    x <- .do_scale()
    .do_agg()
  }, last = {
    z <- .do_agg()
    .do_scale()
  }, never = {
    .do_agg()
  })
  if (length(by) == 1) 
    z <- t(z)
  if (scale != "never" && !(assay == "counts" && fun == "sum")) {
    qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  }
  else lgd_aes <- list()
  lgd_aes$title_gp <- grid::gpar(fontsize = 10, fontface = "bold", 
                                 lineheight = 0.8)
  if (!isFALSE(row_anno)) {
    left_anno <- switch(by[1], sample_id = .anno_factors(x, 
                                                         levels(x$sample_id), row_anno, "row"), .anno_clusters(x, 
                                                                                                               k, m, k_pal, m_pal))
  }
  else left_anno <- NULL
  if (!isFALSE(col_anno) && length(by) == 2) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, 
                              "colum")
  }
  else top_anno <- NULL
  if (bars) {
    right_anno <- CATALYST:::.anno_counts(x[[by[1]]], perc)
  }
  else right_anno <- NULL
  if (bin_anno) {
    cell_fun <- function(j, i, x, y, ...) grid.text(gp = grid::gpar(fontsize = 8), 
                                                    sprintf("%.2f", z[i, j]), x, y)
  }
  else cell_fun <- NULL
  a <- ifelse(assay == "exprs", "expression", assay)
  f <- switch(fun, median = "med", fun)
  hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
                                            "scaled", a), last = sprintf("%s %s\n%s", "scaled", 
                                                                         fun, a), never = paste(fun, a, sep = "\n"))
  if (length(by) == 2) {
    col_title <- features
  }
  else if (length(features) == 1 && features %in% c("type", 
                                                    "state")) {
    col_title <- paste0(features, "_markers")
  }
  else col_title <- ""
  ComplexHeatmap::Heatmap(matrix = z, name = hm_title, col = circlize::colorRamp2(seq(min(z), 
                                                                                      max(z), l = n <- 100), grDevices::colorRampPalette(hm_pal)(n)), 
                          column_title = col_title, column_title_side = ifelse(length(by) == 
                                                                                 2, "top", "bottom"), cell_fun = cell_fun, cluster_rows = row_clust, 
                          cluster_columns = col_clust, show_row_dend = row_dend, 
                          show_column_dend = col_dend, clustering_distance_rows = distance, 
                          clustering_method_rows = linkage, clustering_distance_columns = distance, 
                          clustering_method_columns = linkage, show_row_names = (is.null(left_anno) || 
                                                                                   isTRUE(by == "sample_id")) && !perc, row_names_side = ifelse(by[1] == 
                                                                                                                                                  "cluster_id" || isFALSE(row_anno) && !row_dend || 
                                                                                                                                                  isFALSE(row_clust), "left", "right"), top_annotation = top_anno, 
                          left_annotation = left_anno, right_annotation = right_anno, 
                          rect_gp = gpar(col = "white"), heatmap_legend_param = lgd_aes)
}

plotDiffHeatmap <- function (x, y, k = NULL, top_n = 20, fdr = 0.05, lfc = 1, all = FALSE, 
                             sort_by = c("padj", "lfc", "none"), y_cols = list(padj = "p_adj", 
                                                                               lfc = "logFC", target = "marker_id"), assay = "exprs", 
                             fun = c("median", "mean", "sum"), normalize = TRUE, col_anno = TRUE, 
                             row_anno = TRUE, hm_pal = NULL, fdr_pal = c("lightgrey", 
                                                                         "lightgreen"), lfc_pal = c("blue3", "white", "red3")) 
{
  fun <- match.arg(fun)
  sort_by <- match.arg(sort_by)
  args <- as.list(environment())
  defs <- as.list(formals("plotDiffHeatmap")$y_cols[-1])
  miss <- !names(defs) %in% names(args$y_cols)
  if (any(miss)) 
    y_cols <- args$y_cols <- c(args$y_cols, defs[miss])[names(defs)]
  .check_args_plotDiffHeatmap(args)
  stopifnot(y_cols[[sort_by]] %in% names(y))
  y_cols <- y_cols[y_cols %in% names(y)]
  if (is.null(k)) {
    kids <- levels(y$cluster_id)
    same <- vapply(cluster_codes(x), function(u) identical(levels(u), 
                                                           kids), logical(1))
    if (!any(same)) 
      stop("Couldn't match any clustering", " in input data 'x' with results in 'y'.")
    k <- names(cluster_codes(x))[same][1]
  }
  else {
    k <- .check_k(x, k)
  }
  x$cluster_id <- cluster_ids(x, k)
  y <- data.frame(y, check.names = FALSE)
  y <- mutate_if(y, is.factor, as.character)
  if (any(rownames(x) %in% unlist(y))) {
    features <- intersect(rownames(x), y[[y_cols$target]])
    if (length(features) == 0) 
      stop("Couldn't match features between", " results 'y' and input data 'x'.")
    i <- y[[y_cols$target]] %in% rownames(x)
    type <- "ds"
  }
  else {
    i <- TRUE
    type <- "da"
  }
  y <- rename(y, target = y_cols$target, padj = y_cols$padj, 
              lfc = y_cols$lfc)
  i <- i & !is.na(y$padj) & y$cluster_id %in% levels(x$cluster_id)
  if (!all) {
    i <- i & y$padj < fdr
    if (!is.null(y$lfc)) 
      i <- i & abs(y$lfc) > lfc
  }
  y <- y[i, , drop = FALSE]
  if (nrow(y) == 0) 
    stop("No results remaining;", " perhaps 'x' or 'y' has been filtered,", 
         " or features couldn't be matched.")
  if (sort_by != "none") {
    o <- order(abs(y[[sort_by]]), decreasing = (sort_by == 
                                                  "lfc"))
    y <- y[o, , drop = FALSE]
  }
  if (top_n > nrow(y)) 
    top_n <- nrow(y)
  top <- y[seq_len(top_n), ]
  if (!isFALSE(col_anno)) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, 
                              "column")
  }
  else top_anno <- NULL
  if (is.null(hm_pal)) 
    hm_pal <- rev(brewer.pal(11, ifelse(type == "ds", "RdYlBu", 
                                        "RdBu")))
  if (row_anno) {
    s <- factor(ifelse(top$padj < fdr, "yes", "no"), levels = c("no", 
                                                                "yes"))
    if (!is.null(top$lfc)) {
      lfc_lims <- range(top$lfc, na.rm = TRUE)
      if (all(lfc_lims > 0)) {
        lfc_brks <- c(0, lfc_lims[2])
        lfc_pal <- lfc_pal[-1]
      }
      else if (all(lfc_lims < 0)) {
        lfc_brks <- c(lfc_lims[1], 0)
        lfc_pal <- lfc_pal[-3]
      }
      else lfc_brks <- c(lfc_lims[1], 0, lfc_lims[2])
      lfc_anno <- top$lfc
      anno_cols <- list(logFC = colorRamp2(lfc_brks, lfc_pal))
    }
    else {
      lfc_anno <- NULL
      anno_cols <- list()
    }
    names(fdr_pal) <- levels(s)
    anno_cols$significant <- fdr_pal
    right_anno <- rowAnnotation(logFC = lfc_anno, significant = s, 
                                foo = row_anno_text(scientific(top$padj, 2), gp = gpar(fontsize = 8)), 
                                col = anno_cols, gp = gpar(col = "white"), show_annotation_name = FALSE, 
                                simple_anno_size = unit(4, "mm"))
  }
  else right_anno <- NULL
  switch(type, da = {
    ns <- table(x$cluster_id, x$sample_id)
    fq <- prop.table(ns, 2)
    fq <- fq[top$cluster_id, ]
    y <- as.matrix(unclass(fq))
    if (normalize) y <- .z_normalize(asin(sqrt(y)))
    Heatmap(matrix = y, name = paste0("normalized\n"[normalize], 
                                      "frequency"), col = hm_pal, na_col = "lightgrey", 
            rect_gp = gpar(col = "white"), cluster_rows = FALSE, 
            cluster_columns = FALSE, row_names_side = "left", 
            top_annotation = top_anno, right_annotation = right_anno)
  }, ds = {
    y <- assay(x, assay)
    cs <- .split_cells(x, c("cluster_id", "sample_id"))
    z <- t(mapply(function(k, g) vapply(cs[[k]], function(cs) {
      if (length(cs) == 0) return(NA)
      get(fun)(y[g, cs, drop = FALSE])
    }, numeric(1)), k = top$cluster_id, g = top$target))
    rownames(z) <- sprintf("%s(%s)", top$target, top$cluster_id)
    if (normalize) z <- .z_normalize(z)
    Heatmap(matrix = z, name = paste0("z-normalized\n"[normalize], 
                                      "expression"), col = hm_pal, cluster_rows = FALSE, 
            cluster_columns = FALSE, top_annotation = top_anno, 
            row_names_side = "left", rect_gp = gpar(col = "white"), 
            right_annotation = right_anno, heatmap_legend_param = list(title_gp = gpar(fontsize = 10, 
                                                                                       fontface = "bold", lineheight = 0.8)))
  })
}



