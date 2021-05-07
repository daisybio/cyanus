# BENCHMARKING SCRIPT

# run the script like this for example: 
#       Rscript benchmarking_script.R /nfs/home/students/l.arend/data/covid_spiked/sce_spiked_clustered_ds_full.rds /nfs/home/students/l.arend/data/covid_spiked TRUE FALSE

# first argument -> path of the sce object or objects -> can be a file or a directory
# second argument -> output directory
# third argumetn -> if it should be timed
# fourth argument -> if run parallel


library.path <- .libPaths()[1]
print(library.path)

library(fs)
library(CATALYST, lib.loc=library.path)

# source all files
sapply(list.files("../functions", full.names = TRUE), source)

# save arguments
args <- commandArgs(TRUE)
scePath <- args[1]
outputPath <- args[2]
timed <- as.logical(args[3])
runParallel <- as.logical(args[4])

# check if scePath is file or directory
if (file.exists(scePath) && !dir.exists(scePath)){
  sceFiles <- list(scePath)
} else if (dir.exists(scePath)){
  sceFiles <- list.files(scePath, full.names=TRUE)
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
  param <- MulticoreParam(workers = 18, progressbar = T)
  register(param)
}

for (sceFile in sceFiles){
  #read SCE
  sce <- readRDS(sceFile)

  # run all methods
  results <- runDS(sce,
                   clustering_to_use = "all",
                   contrast_vars = "base_spike",
                   markers_to_test = "state",
                   ds_methods = c("diffcyt-DS-limma"),
                                  #"diffcyt-DS-LMM"),
                                  #"BEZI",
                                  #"ZAGA",
                                  #"ZAIG",
                                  #"sceEMD",
                                  #"hurdleBeta",
                                  #"CytoGLMM"),
                   design_matrix_vars = c("patient_id", "base_spike"),
                   fixed_effects = "base_spike",
                   random_effects = "patient_id",
                   parallel = runParallel,
                   sceEMD_nperm = 500,
                   sceEMD_binsize = 0,
                   time_methods = timed)

  # save the results of the methods
  res <- results[["results"]]

  #only possible for max. 4 sets
  #createVennDiagram(res, DS=T, 0.05, columns = c("diffcyt-DS-limma","diffcyt-DS-LMM","sceEMD", "hurdleBeta")
  res <- data.table::rbindlist(sapply(res, data.table::as.data.table), fill = T, idcol="method")

  objectToSave <- list(results = res)

  #res[p_adj <= 0.05]
  #library(ggplot2)

  # if timed, also save the times of the methods
  if (timed){
    times <- data.table::data.table(
      method = character(),
      user = numeric(),
      system = numeric(),
      elapsed = numeric()
    )

    times <- results[["times"]]
    for(method in names(times)){
      times <- rbind(times,
                  data.table::data.table(
                    method = method,
                    user = as.vector(times[[method]])[1],
                    system = as.vector(times[[method]])[2],
                    elapsed = as.vector(times[[method]])[3]
                  ))
    }

    objectsToSave <- append(objectToSave, times = times)
  }


  # create output file name
  fileName <- strsplit(path_file(sceFile), ".rds")[[1]]
  if (timed){
    add <- "_res_timed.rds"
  } else {
    add <- "_res.rds"
  }
  outputFile <- paste0(outputFile, "/", fileName, add)

  # save file
  saveRDS(objectToSave, outputFile)

}






#ggplot(dt, aes(x = user, y = system, color = method))+
#  geom_point()+
#  theme_bw()
