# BENCHMARKING SCRIPT

# run the script like this for example: 
#       Rscript benchmarking_script.R /nfs/home/students/l.arend/data/covid_spiked/sce_spiked_clustered_ds_full.rds /nfs/home/students/l.arend/data/covid_spiked TRUE FALSE

# first argument -> path of the sce object or objects -> can be a file or a directory
# second argument -> output directory
# third argument -> if it should be timed
# fourth argument -> if run parallel


setwd("..")
source("renv/activate.R")

library(fs)
library(CATALYST)

# source all files
sapply(list.files("functions", full.names = TRUE), source)

# save arguments
args <- commandArgs(TRUE)
# args <- c("/nfs/home/students/l.arend/data/cytoGLMM_simulated/simulated_cytoGLMM_1000_cells.rds", "/nfs/home/students/ga89koc/hiwi/cytof/DEComparison/simulatedCytoGLMM", "condition", "patient_id", "FALSE", "TRUE")
scePath <-  args[1] #"/nfs/home/students/l.arend/data/covid_spiked/downsampled_files/"
outputPath <- args[2] # "DEComparison/"
condition <- args[3]
random_effect <- args[4]
timed <- as.logical(args[5]) # TRUE
runParallel <- as.logical(args[6]) # FALSE

# check if scePath is file or directory
if (file.exists(scePath) && !dir.exists(scePath)){
  sceFiles <- scePath
} else if (dir.exists(scePath)){
  sceFiles <- list.files(scePath, pattern = , "\\.rds$", full.names=TRUE)
} else {
  stop("There is no file or directory with the given name!")
}

# check if the outputPath is a directory
if (!dir.exists(outputPath)){
 stop("You must specify a output directory and not a file or something else!") 
}

# run parallel
if(runParallel){
  library(BiocParallel)
  param <- MulticoreParam(workers = 20, progressbar = T)
  register(param)
}

for (sceFile in sceFiles){
  
  
  # create output file name
  fileName <- strsplit(path_file(sceFile), ".rds")[[1]]
  if (timed) {
    add <- "_res_timed.rds"
  } else {
    add <- "_res.rds"
  }
  outputFile <- paste0(outputPath, "/", fileName, add)
  if (file.exists(outputFile)) next
  
  #read SCE
  sce <- readRDS(sceFile)
  
  # set all none markers to state
  old_classes <- CATALYST::marker_classes(sce)
  SummarizedExperiment::rowData(sce)$marker_class[old_classes == "none"] <- "state"
  
  # run all methods
  results <- runDS(sce,
                   clustering_to_use = "all",
                   contrast_vars = condition,
                   markers_to_test = c("state", "type"),
                   ds_methods = c("diffcyt-DS-limma",
                                  "diffcyt-DS-LMM",
                                  "BEZI",
                                  "ZAGA",
                                  # "ZAIG",
                                  # "hurdleBeta",
                                  "sceEMD",
                                  "CytoGLMM",
                                  "CytoGLM",
                                  "logRegression", 
                                  "wilcoxon_median", 
                                  "kruskal_median"
                   ),
                   design_matrix_vars = c(random_effect, condition),
                   fixed_effects = condition,
                   random_effects = random_effect,
                   parallel = runParallel,
                   sceEMD_nperm = 500,
                   sceEMD_binsize = 0,
                   time_methods = timed)

  # save the results of the methods
  res <- results[["results"]]

  #only possible for max. 4 sets
  #createVennDiagram(res, DS=T, 0.05, columns = c("diffcyt-DS-limma","diffcyt-DS-LMM","sceEMD", "hurdleBeta")
  res <- data.table::rbindlist(sapply(res, data.table::as.data.table, simplify = FALSE), fill = T, idcol="method")

  objectToSave <- list(results = res)

  #res[p_adj <= 0.05]
  #library(ggplot2)

  # if timed, also save the times of the methods
  if (timed) {
    times <- results[["times"]]
    times <- data.table::rbindlist(sapply(times, function(x) as.list(x), simplify = FALSE), idcol = "method")
    objectToSave$times <- times
  }

  # save file
  saveRDS(objectToSave, outputFile)

}






#ggplot(dt, aes(x = user, y = system, color = method))+
#  geom_point()+
#  theme_bw()
