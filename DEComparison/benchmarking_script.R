sce <- readRDS("/nfs/home/students/l.arend/data/covid_spiked/sce_spiked_clustered_ds_full.rds")
sapply(list.files("./functions", full.names = TRUE), source)
library(BiocParallel)
param <- MulticoreParam(workers = 18, progressbar = T)
register(param)
results <- runDS(sce, 
                 clustering_to_use = "all",
                 contrast_vars = "base_spike", 
                 markers_to_test = "state", 
                 ds_methods = c("diffcyt-DS-limma",
                                "diffcyt-DS-LMM", 
                                "sceGAMLSS",
                                "sceEMD",
                                "hurdleBeta",
                                "CytoGLMM"), 
                 design_matrix_vars = c("patient_id", "base_spike"), 
                 fixed_effects = "base_spike", 
                 random_effects = "patient_id", 
                 parallel = T,
                 sceEMD_nperm = 500,
                 sceEMD_binsize = 0, 
                 time_methods = T)

res <- results[["results"]]
times <- results[["times"]]
#only possible for max. 4 sets
#createVennDiagram(res, DS=T, 0.05)
dt <- data.table::rbindlist(sapply(res, data.table::as.data.table), fill = T, idcol="method")
dt[p_adj <= 0.05]
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
