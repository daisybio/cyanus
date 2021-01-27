library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)

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
makePatientSelection <- function(sce,deselected_patients){
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

SigEMD <- function(sce, k, condition, Hur_gene=NULL, binSize=NULL, nperm=100, assay="exprs", seed=1, parallel=FALSE) {
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
  res <- lapply(levels(cluster_ids), function(cluster_id) {
    print(sprintf("calculating SigEMD for cluster %s", cluster_id))
    data <- assay(sce[, cluster_ids == cluster_id], assay)
    
    
    data <- dataclean((abs(data)+data)/2)
    colnames(data) <- as.character(seq.int(to = ncol(data)))
    
    condition_cluster <- colData(sce[, cluster_ids == cluster_id])[[condition]]
    names(condition_cluster) <- colnames(data)
    
    results <- calculate_single(data =  data,condition =  condition_cluster,Hur_gene = Hur_gene, binSize, nperm=nperm, parallel = parallel)
    
    results$cluster_id <- cluster_id
    
    results
  })
  
  names(res) <- levels(cluster_ids)
  
  return(res)
}


############### Differential Expression ###############

prepDiffExp <- function(sce, contrastVars, colsDesign, colsFixed, colsRandom,
                        method = c( "diffcyt-DA-edgeR", "diffcyt-DA-voom","diffcyt-DA-GLMM", "diffcyt-DS-limma", "diffcyt-DS-LMM") ){
  match.arg(method)
  parameters <- list()
  parameters[["ei"]] <- metadata(sce)$experiment_info
  
  if(method %in% c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DS-limma")){
    
    parameters[["design"]] <- diffcyt::createDesignMatrix(parameters[["ei"]], cols_design = colsDesign)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, parameters[["design"]], designMatrix = T)
    
  }else{
    parameters[["formula"]] <- diffcyt::createFormula(parameters[["ei"]], cols_fixed = colsFixed, cols_random = colsRandom)
    parameters[["contrast"]] <- createCustomContrastMatrix(sce, contrastVars, diffcyt::createDesignMatrix(parameters[["ei"]], cols_design = colsFixed), designMatrix = T)
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


