library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)
library(ggvenn)

#create tile heatmap of all methods that were performed in runDA / runDS

createVennHeatmap <- function(res, DS = T, fdr_threshold = 0.05, columns = NULL){
  colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  results <-
    data.table::rbindlist(sapply(res[names(res) != "effect_size"], as.data.table),
      idcol = "method", fill = TRUE)
  if(DS){
    eff_r <- res[["effect_size"]]
    results[, feature := marker_id]#paste0(marker_id, "(", cluster_id, ")")]
    eff_r[, marker_id := sapply(strsplit(eff_r$group2,'::'), "[", 1)]
    eff_r[, feature := marker_id]#paste0(marker_id, "(", cluster_id, ")")]
    eff_r <- eff_r[feature %in% results$feature]
    label <- "Marker"# (Cluster)"
  }else{
    results[, feature := cluster_id]
    label <- "Cluster ID"
  }
  results[, significant := p_adj < fdr_threshold]
  g <- ggplot(results, aes(feature, method))+ 
    geom_tile(aes(fill = significant), color = "white", size = 1)+
    xlab(label = label) +
    ylab("Method") +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    scale_fill_manual(values = colorBlindBlack8[c(7, 3, 1)],
                      name = "Significant",
                      na.value = "transparent")
  if(DS){
    g <- g+ facet_wrap(~ cluster_id) +
      ggside::geom_xsidetile(data=eff_r, aes(y=overall_group, xfill=magnitude), color="white", size=0.2) + 
      ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,2,6)], name='effect size\nmagnitude', na.value="transparent")
  }
  return(g)
}

# create venn diagram of all methods that were performed in runDA / runDS
createVennDiagram <- function(res, DS = T, fdr_threshold = 0.05, columns = NULL) {
  input_venn <- list()
  if(DS){
    feature <- "marker_id"
  }else{
    feature <- "cluster_id"
  }
  # take of each method the data table containing the pvalues
  for (ds_method in names(res)) {
    
    if(feature == "cluster_id"){
      result <-
        data.frame(res[[ds_method]])[c(feature, "p_val", "p_adj")]
      featureNew <- "cluster_id"
    }else{
      result <-
        data.frame(res[[ds_method]])[c("cluster_id", feature, "p_val", "p_adj")]
      if(result$cluster_id[1] != "all"){
        result$marker_id_joined <- paste0(result$marker_id, "(", result$cluster_id, ")")
        featureNew <- "marker_id_joined"
      }else{
        featureNew <- "marker_id"
      }
      
    }
    result$significant <- result$p_adj < fdr_threshold
    significants <-
      unlist(subset(
        result,
        significant == TRUE,
        select = c(get(featureNew)),
        use.names = FALSE
      ))
    input_venn[[ds_method]] <- significants
  }
  venn <- ggvenn::ggvenn(input_venn, columns = columns, show_elements = TRUE, label_sep ="\n", fill_alpha = 0.3, set_name_size = 6, text_size = 4)
  return(venn)
}

createVennDiagramOld <- function(res) {
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
