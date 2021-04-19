library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)
library(ggvenn)

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