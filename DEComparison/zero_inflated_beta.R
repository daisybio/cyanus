library(CATALYST)
library(gamlss)
library(data.table)

source("functions/diffcyt_functions.R")
source("functions/cluster_functions.R")
source("functions/de_functions.R")
source("functions/cytoGLMM_functions.R")
source("functions/ZIBseq_functions.R")

# cytoGLMM simulator (5 of 20 markers significant)
sce_cytoGLMM <- simulateSCE(n_true = 5, n_markers = 10)

CATALYST::plotExprs(sce_cytoGLMM, color_by = "condition")

CATALYST::plotExprHeatmap(sce_cytoGLMM, features = NULL, scale="never")

result <- zibSeq(sce = sce_cytoGLMM, condition = "condition", random_effect = "patient_id")


result_no_random <- zibSeq(sce = sce_cytoGLMM, condition = "condition")
