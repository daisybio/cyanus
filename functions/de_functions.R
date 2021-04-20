#the EMD functions are in sceEMD.R
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
  bool <- sapply(names, function(x){
    any(sapply(contrastVars, function(y){
      grepl(y,x, fixed = T )
    }))
  })
  bool <- as.numeric(bool)
  return(bool)
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
    source("SigEMD/sigEMD.R")
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