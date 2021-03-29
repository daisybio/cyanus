devtools::install_github("ChristofSeiler/CytoGLMM")

library("CytoGLMM")
library("magrittr")
library(dplyr)

## GENERATED DATA
set.seed(23)
df <- generate_data()
df

protein_names <- names(df)[3:12]
df %<>% dplyr:: mutate_at(protein_names, function(x) asinh(x/5))
glmm_fit <- CytoGLMM::cytoglmm(df,
                               protein_names = protein_names,
                               condition = "condition",
                               group = "donor",
                               num_cores=1)
glmm_fit
plot(glmm_fit)
summary(glmm_fit)
summary(glmm_fit) %>% dplyr::filter(pvalues_adj < 0.05)


glm_fit <- CytoGLMM::cytoglm(df,
                             protein_names = protein_names,
                             condition = "condition",
                             group = "donor",
                             num_boot = 1000)
glm_fit
summary(glm_fit)



## TRY WITH PBMC DATA
data(PBMC_panel, PBMC_md, PBMC_fs)
sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md, transform=T)

pbmc_data <- assays(sce)$exprs

marker_names <- rowData(sce)$marker_name
marker_names <- sapply(marker_names, function(marker) {
  gsub("[^[:alnum:]]", "", marker)
})

pbmc_data <- as.data.frame(t(pbmc_data))
colnames(pbmc_data) <- marker_names

pbmc_data$donor <- colData(sce)$patient_id
pbmc_data$condition <- colData(sce)$condition


pbmc_data$donor <- as.character(pbmc_data$donor)

glmm_fit <- CytoGLMM:: cytoglmm(pbmc_data,
                                protein_names = marker_names,
                                condition = "condition",
                                group = "donor",
)

summary(glmm_fit)
plot(glmm_fit)

glm_fit <- CytoGLMM:: cytoglm(pbmc_data,
                              protein_names = marker_names,
                              condition = "condition",
                              group = "donor",
                              num_boot = 1000)

summary(glm_fit)
plot(glm_fit)


## TRY WITH PLATELETS DATA
sce_platelets <- readRDS("/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/cytof/data/platelets_small/sce.rds")
sce_platelets
data <- assays(sce_platelets)$exprs
marker_names <- rowData(sce_platelets)$marker_name
marker_names <- sapply(marker_names, function(marker) {
  gsub("[^[:alnum:]]", "", marker)
})
data <- as.data.frame(t(data))
colnames(data) <- marker_names
data$donor <- colData(sce_platelets)$patient_id
data$condition <- colData(sce_platelets)$platelets

data$donor <- as.character(data$donor)
glmm_fit <- CytoGLMM:: cytoglmm(data,
                                protein_names = marker_names,
                                condition = "condition",
                                group = "donor")
summary(glmm_fit)
results <- summary(glmm_fit)
p.adjust(results$pvalues_unadj, method = "BH")


glm_fit <- CytoGLMM::cytoglm(data,
                             protein_names = marker_names,
                             condition = "condition",
                             group = "donor",
                             num_boot = 1000)
summary(glm_fit)

## SAMPLE PLATELETS DATA
ei_platelets <- ei(sce_platelets)
ei_platelets$platelets <- sample(ei_platelets$platelets)
platelets <- rep(ei_platelets$platelets, times = as.numeric(ei_platelets$n_cells))
data$condition <- platelets

glmm_fit <- CytoGLMM:: cytoglmm(data,
                                protein_names = marker_names,
                                condition = "condition",
                                group = "donor")
summary(glmm_fit)
results <- summary(glmm_fit)
p.adjust(results$pvalues_unadj, method = "BH")

