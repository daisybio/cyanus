#PBMC

sapply(list.files("functions", full.names = TRUE), source)

library(BiocParallel)
param <- MulticoreParam(workers = 20, progressbar = T)
register(param)

sce <- readRDS("/nfs/home/students/l.arend/data/pbmc_all/bigPBMC_SCE.rds")
sce <- clusterSCE(sce)

CATALYST::plotPbExprs(sce, color_by = "condition")

conditions <- unique(ei(sce)$condition)
combinations <- as.data.table(combn(conditions,2))

for (i in 1:ncol(combinations)){
  c1 <- combinations[1,..i]
  c2 <- combinations[2,..i]
  sce_tmp <- filterSCE(sce, condition %in% c(c1,c2))
  bpPlot <- plotPbExprs(sce_tmp, color_by = "condition", features = c("pp38", "pStat5", "pSHP2", "pZap70", "pSlp76", "pBtk", "pS6"))
  exprsPlot <- plotExprs(sce_tmp, color_by = "condition", features = c("pp38", "pStat5", "pSHP2", "pZap70", "pSlp76", "pBtk", "pS6"))
  plt <- grid.arrange(bpPlot, exprsPlot, ncol=2)
  ggplot2::ggsave(paste0("/nfs/home/students/l.arend/data/pbmc_all/Plots/",c1,"_", c2,".png"),plt, width=12, height=8)
}
 



comparisons <- data.table(
  condition1 = c("BCR-XL", "BCR-XL", "IFN-a", "IFN-a", "IFN-a", "IFN-a", "IFN-a", "IFN-a", "IFN-g", "IFN-g", "IL-2", "IL-12"),
  condition2 = c("Il-3", "LPS", "IFN-g", "Il-3", "IL-12", "LPS", "PMA","Reference", "LPS", "PMA", "PMA", "LPS"),
  marker = c("pBtk", "pBtk", "pSHP2", "pSHP2" ,"pSHP2", "pSHP2", "pSHP2", "pSHP2", "pSTAT5", "pSTAT5", "pSHP2", "pBtk")
)

for (i in 1:nrow(comparisons)){
  c1 <- as.character(comparisons[i,1])
  c2 <- as.character(comparisons[i,2])
  marker <- comparisons[i,3]
  sce_tmp <- filterSCE(sce, condition %in% c(c1,c2))
  
  results <- runDS(sce_tmp, 
                   clustering_to_use = "all", 
                   contrast_vars = "condition", 
                   markers_to_test = c("type", "state"), 
                   ds_methods = c("diffcyt-DS-limma", "diffcyt-DS-LMM","sceEMD"), 
                   fixed_effects = "condition",
                   random_effects = "patient_id",
                   parallel = TRUE,
                   design_matrix_vars = c("patient_id", "condition"))
  
  #bpPlot <- plotPbExprs(sce_tmp, color_by = "condition", features = c("pp38", "pStat5", "pSHP2", "pZap70", "pSlp76", "pBtk", "pS6"))
  #exprsPlot <- plotExprs(sce_tmp, color_by = "condition", features = c("pp38", "pStat5", "pSHP2", "pZap70", "pSlp76", "pBtk", "pS6"))
  venn <- createVennDiagram(results) + ggtitle(paste(c1, "vs", c2, "checking for "), marker)
  #plt <- grid.arrange(bpPlot, exprsPlot,venn, ncol=3)
  ggplot2::ggsave(paste0("/nfs/home/students/l.arend/data/pbmc_all/Plots/",c1,"_", c2,"_venn.png"),venn, width=12, height=8)
}

