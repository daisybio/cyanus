library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)

###Visualization

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
