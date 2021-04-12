
runCytoGLMM <- function(sce, condition, group){
  data <- as.data.frame(t(assays(sce)$exprs))
  marker_names <- rowData(sce)$marker_name
  marker_names <- sapply(marker_names, function(marker) {
    gsub("[^[:alnum:]]", "", marker)
  })
  colnames(data) <- marker_names
  data$donor <- colData(sce)[[group]]
  data$condition <- colData(sce)[[condition]]
  glmm_fit <- CytoGLMM:: cytoglmm(data,
                                  protein_names = marker_names,
                                  condition = "condition",
                                  group = "donor")
  return(glmm_fit)
}


runCytoGLM <- function(sce, condition, group){
  data <- as.data.frame(t(assays(sce)$exprs))
  marker_names <- rowData(sce)$marker_name
  marker_names <- sapply(marker_names, function(marker) {
    gsub("[^[:alnum:]]", "", marker)
  })
  colnames(data) <- marker_names
  data$donor <- colData(sce)[[group]]
  data$condition <- colData(sce)[[condition]]
  glm_fit <- CytoGLMM:: cytoglm(data,
                                protein_names = marker_names,
                                condition = "condition", 
                                group = "donor",
                                num_boot = 1000,
                                num_cores = 22)
  return(glm_fit)
}
