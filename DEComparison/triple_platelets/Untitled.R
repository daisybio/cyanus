# library(CATALYST)
# sce <- readRDS("/nfs/home/students/jbernett/cytof/PBMC/bigPBMC_SCE.rds")
# 
# sceIFNg <- filterSCE(sce, condition %in% c('Reference', 'IFN-g'))
# 
# dt <- as.data.table(cbind(t(assay(sceIFNg, 'exprs')), as.data.frame(colData(sceIFNg))))
# 
# library(ggplot2)
# ggplot(dt, aes(x = pSHP2, fill = condition, color = condition)) + 
#   geom_density(alpha = .2)
# 
# fisher.test(dt[, table(.(condition, pSHP2 == 0))])
setwd("/nfs/home/students/ga89koc/hiwi/cytof/")
library(CATALYST)

sapply(list.files("functions", full.names = TRUE), source)

x <- readRDS("/nfs/home/students/l.arend/data/platelets/sce_transformed.rds")
triple <- filterSCE(x, therapy == "triple")
triple <- clusterSCE(triple)


old_classes <- CATALYST::marker_classes(triple)
SummarizedExperiment::rowData(triple)$marker_class[old_classes == "none"] <- "state"

condition <- 'platelets'
random_effect <- 'patient_id'
CATALYST::plotPbExprs(triple, color_by = condition, features = c('CD141', 'CD154', 'CD3')) + theme(text = element_text(size = 18))
CATALYST::plotExprs(triple, color_by = condition, features = c('CD141', 'CD154', 'CD3')) + theme(text = element_text(size = 18))


# library(BiocParallel)
# param <- MulticoreParam(workers = 22, progressbar = T)
# register(param)

results <- runDS(triple_downsample,
                 clustering_to_use = "all",
                 contrast_vars = condition,
                 markers_to_test = c("state", "type"),
                 ds_methods = c("diffcyt-DS-limma",
                                "diffcyt-DS-LMM",
                                # "BEZI",
                                "ZAGA",
                                # "ZAIG",
                                # "hurdleBeta",
                                "sceEMD",
                                "CytoGLMM",
                                # "CytoGLM",
                                "logRegression",
                                "wilcoxon_median",
                                "kruskal_median",
                                "t_test"
                 ),
                 design_matrix_vars = c(random_effect, condition),
                 fixed_effects = condition,
                 random_effects = random_effect,
                 parallel = FALSE)
# saveRDS(results, '/nfs/home/students/ga89koc/hiwi/cytof/DEComparison/triple_platelets/triple_results_nonparallel.rds')

results <- data.table::rbindlist(sapply(results, function(x) data.table::as.data.table(x)[,method := NULL], simplify = FALSE), fill = TRUE, idcol = 'method')
results[marker_id == "CD154", .(method, p_adj)]

