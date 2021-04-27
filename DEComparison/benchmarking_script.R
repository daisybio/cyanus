sce <- readRDS("/nfs/home/students/l.arend/data/covid_spiked/sce_spiked_clustered_ds_25.rds")

source("functions/de_functions.R")
source("functions/sceEMD.R")
source("functions/cytoGLMM_functions.R")
source("functions/ZIBseq_functions.R")
source("functions/venn_functions.R")
library(BiocParallel)
param <- MulticoreParam(workers = 18, progressbar = T)
register(param)
results <- runDS(sce, 
                 ds_methods = c("diffcyt-DS-limma","diffcyt-DS-LMM", "sceEMD", "ZIBseq"), 
                 clustering_to_use = "all", contrast_vars = "base_spike", 
                 design_matrix_vars = c("patient_id", "base_spike"), fixed_effects = "base_spike", 
                 random_effects = "patient_id", markers_to_test = c("type", "state"), 
                 sceEMD_condition = "base_spike", binSize = 0, nperm = 500, parallel=T, time_methods = T)

res <- results[["results"]]
times <- results[["times"]]
createVennDiagram(res, DS=T, 0.05)
library(ggplot2)
dt <- data.table::data.table(
  method = character(),
  user = numeric(),
  system = numeric(),
  elapsed = numeric()
)
for(method in names(times)){
  dt <- rbind(dt, 
              data.table::data.table(
                method = method,
                user = as.vector(times[[method]])[1],
                system = as.vector(times[[method]])[2],
                elapsed = as.vector(times[[method]])[3]
              ))
}

ggplot(dt, aes(x = user, y = system, color = method))+
  geom_point()+
  theme_bw()
