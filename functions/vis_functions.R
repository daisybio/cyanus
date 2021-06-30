runCatalystDR <-
  function(sce, 
           dr_chosen = c("UMAP", "TSNE", "PCA", "MDS", "DiffusionMap", "Isomap"),
           cells_chosen = NULL,
           seed = 42,
           feature_chosen = "type",
           assay_chosen = "exprs",
           scale = T,
           k = 5,
           dimensions = 2) {
    match.arg(dr_chosen)
    if (scale == "yes") {
      scale <- T
    } else{
      scale <- F
    }
    if (dr_chosen == "Isomap") {
      set.seed(seed)
      sce <- runIsomap(
        sce,
        cells = cells_chosen,
        features = feature_chosen,
        assay = assay_chosen,
        scale = scale,
        k = k,
        dimensions = dimensions
      )
      
    } else{
      set.seed(seed)
      sce <- CATALYST::runDR(
        x = sce,
        dr = dr_chosen,
        cells = cells_chosen,
        features = feature_chosen,
        assay = assay_chosen,
        scale = scale,
        ncomponents = dimensions
      )
    }
    return(sce)
  }

makeDR <-
  function(sce,
           dr_chosen,
           color_chosen,
           facet_chosen,
           assay_chosen,
           scale,
           dims = c(1,2)) {
    if (color_chosen == "") {
      color_chosen <- NULL
    }
    if (facet_chosen == "") {
      facet_chosen <- NULL
    }
    if (scale == "yes") {
      scale <- T
    } else{
      scale <- F
    }
    g <-
      CATALYST::plotDR(
        sce,
        dr = dr_chosen,
        color_by = color_chosen,
        facet_by = facet_chosen,
        assay = assay_chosen,
        scale = scale,
        dims = dims
      )
    return(g)
  }

#adapted from CATALYST::runDR
runIsomap <- function (x, cells = NULL, features = "type", assay = "exprs", scale = TRUE, k = 5, dimensions = 2) 
{
  library(dimRed)
  library(RANN)
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
  res_isomap <- dimRed::embed(t(selectedCounts), "Isomap", knn = k, ndim = dimensions)
  
  if (is.null(cells)){
    reducedDim(x, "Isomap") <- res_isomapb@data@data
    return(x)
  }
  
  m <- matrix(NA, nrow = ncol(x), ncol = ncol(res_isomap@data@data))
  m[cs, ] <- res_isomap@data@data
  reducedDim(x, "Isomap") <- m
  return(x)
}


renameColorColumn <- function(columnNames, color_by = T) {
  returnVector <- c()
  if ("sample_id" %in% columnNames) {
    returnVector <- c(returnVector, "Sample ID" = "sample_id")
    columnNames <- columnNames[columnNames != "sample_id"]
  }
  if ("patient_id" %in% columnNames) {
    returnVector <- c(returnVector, "Patient ID" = "patient_id")
    columnNames <- columnNames[columnNames != "patient_id"]
  }
  if ("cluster_id" %in% columnNames) {
    returnVector <- c(returnVector, "flowSOM ID" = "cluster_id")
    columnNames <- columnNames[columnNames != "cluster_id"]
    if (color_by) {
      returnVector <-
        c(returnVector, names(cluster_codes(reactiveVals$sce))[-1])
    }
  }
  returnVector <- c(returnVector, columnNames)
  return(returnVector)
}