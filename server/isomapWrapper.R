library(dimRed)
library(RANN)

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
  #old package
  #res_isomap <- do.isomap(t(selectedCounts), preprocess = preprocess, type = type)
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
