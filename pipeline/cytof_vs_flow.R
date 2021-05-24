library(data.table)
library(CATALYST)
library(diffcyt)
library(flowCore)

panel <- fread("~/Downloads/panelNew.csv")
#fwrite(panel, file = "panelNew.csv", quote = F)

fcs_files <- list.files("~/Downloads/FlowRepository_FR-FCM-Z222_files/", full.names = T)
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
metadata <- data.table(
  "file_name" = fcs_files,
  "patient_id" = c("37", "37", "10", "10", "11", "11", "3", "3", "3", "4", "4", "4", "5", "5", "5", "6", "6", "7", "8", "9", "9", "1", "2", "3", "4", "5", "6", "7", "2", "50", "50", "93", "93", "404", "404", "639", "639"),
  "sample_id" = tstrsplit(list.files("~/Downloads/FlowRepository_FR-FCM-Z222_files/"), ".fcs", keep=1)[[1]],
  "condition" = c("Breast", "Breast", rep("PBMC", 27), rep("Melanoma", 4), rep("Ovarian", 4)),
  "cancer" = c("Cancer", "Cancer", rep("Healthy", 27), rep("Cancer", 8)),
  "fresh" = c(rep("fresh", 23), rep("frozen", 5), rep("fresh", 9))
)

sce <- prepData(fcs_set, panel=panel, md = metadata, md_cols = list(file="file_name", id="sample_id", 
                                                                    factors = c("condition", "patient_id", "cancer", "fresh")))
#saveRDS(sce, "data/cytof_vs_flow_sce.rds")
assay(sce, "exprs") <- (assay(sce, "exprs") + abs(assay(sce, "exprs"))) / 2
plotExprs(sce, features = c("CD3", "IgM"))



