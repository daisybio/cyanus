# BENCHMARKING SCRIPT

# run the script like this for example: 
# Rscript benchmarking_script.R DataGeneration/covid/sce_spiked_clustered_ds_full.rds DataGeneration/covid/ TRUE FALSE

# first argument -> path of the sce object or objects -> can be a file or a directory
# second argument -> output directory
# third argument -> condition
# fourth argument -> grouping variable; can be NULL
# fifth argument -> if it should be timed
# sixth argument -> if it should be run in parallel
# seventh argument -> which clustering should be used; if none is specified, cluster "all" is taken


setwd("../..")
source("renv/activate.R")

library(fs)
library(CATALYST)

# source all files
sapply(list.files("functions", full.names = TRUE), source)

# save arguments
args <- commandArgs(TRUE)
#args <- c("DataGeneration/cytoGLMM_simulated", "DEComparison/simulatedCytoGLMM", "condition", "patient_id", "TRUE", "FALSE")
scePath <-  args[1] #"DataGeneration/covid/"
outputPath <- args[2] # "DEComparison/"
condition <- args[3] # condition
random_effect <- args[4] # patient_id
if(tolower(random_effect) == "null"){
  random_effect <- NULL
}
timed <- as.logical(args[5]) # TRUE
runParallel <- as.logical(args[6]) # FALSE
if(length(args) == 7){
  clustering_to_use <- args[7]
}else{
  clustering_to_use <- "all"
}

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
  message(paste('processing sce', sceFile))
  
  # create output file name
  fileName <- strsplit(path_file(sceFile), ".rds")[[1]]
  if (timed) {
    add <- "_res_timed.rds"
  } else {
    add <- "_res.rds"
  }
  outputFile <- paste0(outputPath, "/", fileName, add)
  if (file.exists(outputFile)) {
    # try not to compute everything multiple times
    message('output file found')
    savedObj <- readRDS(outputFile)
    if (is.null(savedObj$eff)) {
      message('adding effect size to output')
      sce <- readRDS(sceFile)
      savedObj$eff <- effectSize(sce, condition, random_effect, clustering_to_use)
      saveRDS(savedObj, outputFile)
    }
    next
  }
  
  #read SCE
  sce <- readRDS(sceFile)
  
  # set all none markers to state
  old_classes <- CATALYST::marker_classes(sce)
  SummarizedExperiment::rowData(sce)$marker_class[old_classes == "none"] <- "state"
  
  # run all methods
  results <- runDS(sce,
                   clustering_to_use = clustering_to_use,
                   contrast_vars = condition,
                   markers_to_test = c("state", "type"),
                   ds_methods = c("diffcyt-DS-limma",
                                  "diffcyt-DS-LMM",
                                  "BEZI",
                                  "ZAGA",
                                  "CyEMD",
                                  "CytoGLMM",
                                  "CytoGLM",
                                  "logRegression",
                                  "wilcoxon_median",
                                  "kruskal_median",
                                  "t_test"
                   ),
                   design_matrix_vars = c(random_effect, condition),
                   fixed_effects = condition,
                   random_effects = random_effect,
                   parallel = runParallel,
                   cyEMD_nperm = 500,
                   cyEMD_binsize = 0,
                   time_methods = timed)

  objectToSave <- list()
  # lets compute and save the effect size:
  objectToSave$eff <- effectSize(sce, condition, random_effect, clustering_to_use)
  
  # if timed, also save the times of the methods
  if (timed) {
    # save the results of the methods
    res <- results[["results"]]
    times <- results[["times"]]
    times <- data.table::rbindlist(sapply(times, function(x) as.list(x), simplify = FALSE), idcol = "method")
    objectToSave$times <- times
  } else 
    res <- results


  #only possible for max. 4 sets
  res <- data.table::rbindlist(sapply(res, data.table::as.data.table, simplify = FALSE), fill = T, idcol="method")
  objectToSave$results <- res
  
  # save file
  saveRDS(objectToSave, outputFile)

}






#ggplot(dt, aes(x = user, y = system, color = method))+
#  geom_point()+
#  theme_bw()
