library(data.table)
library(CATALYST)
library(diffcyt)
library(flowCore)

panel <- fread("~/cytof/panelNew.csv")
#fwrite(panel, file = "panelNew.csv", quote = F)

fcs_files <- list.files("~/cytof/FlowRepository_FR-FCM-Z222_files/", full.names = T)
fcs_list <- sapply(fcs_files, function(x){
  fcs_file <- read.FCS(x)
  fcs_file <- fcs_file[, !colnames(fcs_file) %in% c("Event-length", "Event_length")]
  if(length(colnames(fcs_file)) == 59){
    fcs_file <- fcs_file[, !colnames(fcs_file) %in% c("Width")]
  }
  fcs_file <- fcs_file[, order(colnames(fcs_file))]
  return(fcs_file)
})
fcs_set <- flowSet(fcs_list)
names <- gsub(" ", "_", tstrsplit(list.files("~/cytof/FlowRepository_FR-FCM-Z222_files/"), ".fcs", keep=1)[[1]], fixed = T)
write.flowSet(fcs_set, outdir = "~/cytof/fcs_files", filename = names)

metadata <- data.table(
  "file_name" = list.files("~/cytof/fcs_files", pattern = "fcs"),
  "patient_id" = c("37", "37", "10", "10", "11", "11", "3", "3", "3", "4", "4", "4", "5", "5", "5", "6", "6", "7", "8", "9", "9", "1", "2", "3", "4", "5", "6", "7", "2", "50", "50", "93", "93", "404", "404", "639", "639"),
  "sample_id" = tstrsplit(list.files("~/cytof/fcs_files", pattern = ".fcs"), ".fcs", keep=1)[[1]],
  "condition" = c("Breast", "Breast", rep("PBMC", 27), rep("Melanoma", 4), rep("Ovarian", 4)),
  "cancer" = c("Cancer", "Cancer", rep("Healthy", 27), rep("Cancer", 8)),
  "fresh" = c(rep("fresh", 23), rep("frozen", 5), rep("fresh", 9))
)

sce <- prepData(list.files("~/cytof/fcs_files", full.names = T, pattern = "fcs"), panel=panel, md = metadata, md_cols = list(file="file_name", id="sample_id", 
                                                                    factors = c("condition", "patient_id", "cancer", "fresh")))

assay(sce, "counts") <- (assay(sce, "counts") + abs(assay(sce, "counts"))) / 2
assay(sce, "exprs") <- (assay(sce, "exprs") + abs(assay(sce, "exprs"))) / 2
source("functions/cluster_functions.R")
sce <- clusterSCE(sce)
saveRDS(sce, "~/cytof/cytof_vs_flow_sce.rds")

plotCounts(sce, group_by = "patient_id", color_by = "condition")
pbMDS(sce)
delta_area(sce)
plotExprs(sce, features = c("CD3", "Perforin", "CD103", "CD39"), color_by = "cancer")

sce_dual_triple <- readRDS("/nfs/home/students/l.arend/data/platelets/sce_transformed.rds")
sce_triple <- filterSCE(sce_dual_triple, therapy == "triple")
sce_triple <- filterSCE(sce_triple, !patient_id %in% c("RPS 108", "RPS 099"))

pbMDS(sce_triple, color_by = "platelets")
source("functions/de_functions.R")
source("functions/sceEMD.R")
source("functions/cluster_functions.R")
sce_triple <- addClusterAll(sce_triple)
library(BiocParallel)
param <- MulticoreParam(workers = 22, progressbar = T)
register(param)
results <- runDS(sce_triple, clustering_to_use = "all", contrast_vars = "platelets", markers_to_test = c("type", "state"), 
      ds_methods = c("diffcyt-DS-limma","diffcyt-DS-LMM","sceEMD"), design_matrix_vars = c("patient_id", "platelets"), 
      fixed_effects = "platelets", random_effects = "patient_id", parallel = T, sceEMD_nperm = 2000)
source("functions/venn_functions.R")
createVennDiagram(results)



