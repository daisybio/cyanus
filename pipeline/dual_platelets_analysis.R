# DUAL PLATELETS

sapply(list.files("functions", full.names = TRUE), source)



sce <- readRDS("/nfs/home/students/l.arend/data/platelets_dual/sce.rds")

library(BiocParallel)
param <- MulticoreParam(workers = 22, progressbar = T)
register(param)
results <- runDS(sce, clustering_to_use = "all", contrast_vars = "platelets", markers_to_test = c("type", "state"), 
                 ds_methods = c("diffcyt-DS-limma",
                                "diffcyt-DS-LMM",
                                # "ZAIG",
                                # "hurdleBeta",
                                "sceEMD",
                                "CytoGLMM",
                                "logRegression",
                                "wilcoxon_median",
                                "kruskal_median"
                 ),
                 design_matrix_vars = c("patient_id", "platelets"), 
                 fixed_effects = "platelets", random_effects = "patient_id", parallel = T)

results <- sapply(results, as.data.table)
results_df <- rbindlist(results, fill=TRUE)
