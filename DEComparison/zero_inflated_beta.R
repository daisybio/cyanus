library(CATALYST)
library(gamlss)
library(data.table)

source("functions/diffcyt_functions.R")
source("functions/cluster_functions.R")
source("functions/de_functions.R")
source("functions/cytoGLMM_functions.R")
source("functions/ZIBseq_functions.R")
source("functions/prep_functions.R")

# cytoGLMM simulator (5 of 20 markers significant)
sce_cytoGLMM <- simulateSCE(n_true = 5, n_markers = 10)


################## PBMC DATA ####################
data(PBMC_panel, PBMC_md, PBMC_fs)
sce_pbmc<- prepData(PBMC_fs, PBMC_panel, PBMC_md, transform=T)
sce_pbmc <- clusterSCE(sce_pbmc)
rowData(sce_pbmc)$marker_class <- rowData(sce_pbmc)$marker_class[rowData(sce_pbmc)$marker_class == "none"] <- "state"


# check diffcyt LMM no weights
markers_to_test <- getMarkersToTest(sce_pbmc,"LMM","all")
LMM_res_no_weights <- diffcyt_method(d_input = sce_pbmc,
                                     formula = createFormula(ei(sce_pbmc), cols_fixed = c("condition")),
                                     contrast = createContrast(c(0,1)),
                                     analysis_type = "DS",
                                     method_DS = "diffcyt-DS-LMM",
                                     clustering_to_use = "all",
                                     use_weights = FALSE,
                                     markers_to_test = markers_to_test)

results_LMM <- as.data.table(rowData(LMM_res_no_weights$res))$p_adj
names(results_LMM) <- as.data.table(rowData(LMM_res_no_weights$res))$marker_id

# ZIBSeq
result <- zibSeq(sce = sce_pbmc, condition = "condition")
padj <- p.adjust(result$pvalues, method="BH")
names(padj) <- colnames(exprs)



################# Simulated data set ########################
sce_cytoGLMM <- simulateSCE()

result <- zibSeq(sce = sce_cytoGLMM, condition = "condition")
padj <- p.adjust(result$pvalues, method="BH")
names(padj) <- colnames(exprs)

results_cytoGLMM <- runCytoGLMM(sce_cytoGLMM, "condition", "patient_id")


#################### Covid spiked data #########################
sce_covid_spiked <- readRDS("/nfs/home/students/l.arend/data/covid_spiked/sce_spiked_full.rds")

sce_covid_spiked <- transformData(sce_covid_spiked)
sce_covid_spiked <- clusterSCE(sce_covid_spiked)

rowData(sce_covid_spiked)$marker_class <- rowData(sce_covid_spiked)$marker_class[rowData(sce_covid_spiked)$marker_class == "none"] <- "type"

zib_results <- zibSeq(sce = sce_covid_spiked, condition = "base_spike")

markers_to_test <- getMarkersToTest(sce_covid_spiked,"LMM","all")
LMM_results <- diffcyt_method(d_input = sce_covid_spiked,
                                     formula = createFormula(ei(sce_covid_spiked), cols_fixed = c("base_spike")),
                                     contrast = createContrast(c(0,1)),
                                     analysis_type = "DS",
                                     method_DS = "diffcyt-DS-LMM",
                                     clustering_to_use = "all",
                                     use_weights = FALSE,
                                     markers_to_test = markers_to_test)
