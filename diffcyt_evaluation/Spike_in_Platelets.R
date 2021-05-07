spikeInMarkers <- function(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                           subsetCondition = "covid",
                           subsetConditionValue = "healthy",
                           baselineVar = "platelets",
                           baselineValue = "B",
                           activatedValue = "A",
                           markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"),
                           alpha = 0) {
  library(CATALYST)
  library(data.table)
  sce <- readRDS(pathToSCE)
  sce <- filterSCE(sce, get(subsetCondition) == subsetConditionValue)
  sceBaseline <- filterSCE(sce, get(baselineVar) == baselineValue)
  sceActivated <- filterSCE(sce, get(baselineVar) == activatedValue)
  remove(sce)
  
  # ---------------------------------------------------------------
  # Split each baseline sample into two halves: 'base' and 'spike'
  # ---------------------------------------------------------------
  
  n_cells_baseline <- ei(sceBaseline)$n_cells
  start_next_sample <- c(0, n_cells_baseline)
  start_next_sample <-
    cumsum(start_next_sample[-length(start_next_sample)])
  
  #how many do I have in activated?
  n_cells_activated <- ei(sceActivated)$n_cells
  n_cells_activated - floor(n_cells_baseline / 2)
  # --> enough ? no -> take all there is from the activated condition for patient CVD003A and CVD013A
  
  set.seed(100)
  # generate random indices
  inds <- lapply(c(1:length(n_cells_baseline)), function(i) {
    if (n_cells_activated[i] - floor(n_cells_baseline[i] / 2) < 0) {
      n <- n_cells_activated[i]
      i_spike <- sort(sample(seq_len(n_cells_baseline[i]), n))
      i_base <- setdiff(seq_len(n_cells_baseline[i]), i_spike)
    } else{
      n <- n_cells_baseline[i]
      i_base <- sort(sample(seq_len(n), floor(n / 2)))
      i_spike <- setdiff(seq_len(n), i_base)
    }
    list(base = i_base, spike = i_spike)
  })
  
  inds2 <- lapply(seq(1:length(inds)), function(x) {
    base_inds <- inds[[x]][["base"]] + start_next_sample[x]
    spike_inds <- inds[[x]][["spike"]] + start_next_sample[x]
    list(base = base_inds, spike = spike_inds)
  })
  
  inds_base <- unlist(lapply(inds2, function(l)
    l[[1]]))
  inds_spike <- unlist(lapply(inds2, function(l)
    l[[2]]))
  
  conditionCol <-
    data.table::data.table(index = c(1:ncol(sceBaseline)), sample_id = sceBaseline$sample_id)
  conditionCol <-
    conditionCol[, base_spike := ifelse(index %in% inds_base, "base", "spike")]
  conditionCol <-
    conditionCol[, sample_id := ifelse(index %in% inds_base,
                                       paste0(sample_id, "_base"),
                                       paste0(sample_id, "_spike"))]
  colData(sceBaseline)$base_spike <-
    as.factor(conditionCol$base_spike)
  colData(sceBaseline)$sample_id <- as.factor(conditionCol$sample_id)
  # -------------------------------------------------------------------------------------
  # Replace cells in 'spike' samples with cells from stimulated condition but
  # only CD63 and CD62P measurements
  # -------------------------------------------------------------------------------------
  
  # generate random indices from activated
  
  #how many are needed per sample?
  needed_per_sample <- sapply(inds2, function(x) {
    return(length(x[["spike"]]))
  })
  
  start_of_next_sample_A <- c(0, n_cells_activated)
  start_of_next_sample_A <-
    cumsum(start_of_next_sample_A[-length(start_of_next_sample_A)])
  
  ## select correct number of activated cells for each sample and generate random indices
  set.seed(100)
  indsActivated <-
    unlist(sapply(c(1:length(
      needed_per_sample
    )), function(i) {
      needed <- needed_per_sample[i]
      got <- n_cells_activated[i]
      if ((got - needed) < 0) {
        needed <- got
      }
      sort(sample(seq_len(got), needed)) + start_of_next_sample_A[i]
    }))
  for (marker in markersToSpike) {
    message(paste("Replacing marker", marker, "..."))
    if(alpha == 0){
      assays(sceBaseline)$counts[marker, inds_spike] <-
        assays(sceActivated)$counts[marker, indsActivated]
    }else{
      assays(sceBaseline)$counts[marker, inds_spike] <-
        5 * sinh(
          asinh(assays(sceActivated)$counts[marker, indsActivated]/5) - 
          alpha * (
            asinh(assays(sceActivated)$counts[marker, indsActivated] / 5) - 
              asinh(assays(sceBaseline)$counts[marker, inds_spike] / 5)
            )
          )
    }
  }
  
  n_cells <- as.vector(n_cells(sceBaseline))
  
  exp_info <- data.table(
    sample_id = levels(sceBaseline$sample_id),
    subsetCondition = subsetConditionValue,
    patient_id = rep(ei(sceBaseline)$patient_id, each = 2),
    baselineVar = baselineValue,
    base_spike = tstrsplit(
      levels(sceBaseline$sample_id),
      "_",
      fixed = T,
      keep = 2
    )[[1]],
    batch = rep(ei(sceBaseline)$batch, each = 2),
    acquisition_date = rep(ei(sceBaseline)$acquisition_date, each = 2),
    n_cells = n_cells
  )
  setnames(exp_info,
           c("subsetCondition", "baselineVar"),
           c(subsetCondition, baselineVar))
  
  metadata(sceBaseline)$experiment_info <- exp_info
  return(sceBaseline)
}

sce <- spikeInMarkers(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                      subsetCondition = "covid",
                      subsetConditionValue = "healthy",
                      baselineVar = "platelets",
                      baselineValue = "B",
                      activatedValue = "A",
                      markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"), 
                      alpha = 0)
sce25 <- spikeInMarkers(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                      subsetCondition = "covid",
                      subsetConditionValue = "healthy",
                      baselineVar = "platelets",
                      baselineValue = "B",
                      activatedValue = "A",
                      markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"), 
                      alpha = 0.25)
sce50 <- spikeInMarkers(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                      subsetCondition = "covid",
                      subsetConditionValue = "healthy",
                      baselineVar = "platelets",
                      baselineValue = "B",
                      activatedValue = "A",
                      markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"), 
                      alpha = 0.5)
sce75 <- spikeInMarkers(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                      subsetCondition = "covid",
                      subsetConditionValue = "healthy",
                      baselineVar = "platelets",
                      baselineValue = "B",
                      activatedValue = "A",
                      markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"), 
                      alpha = 0.75)
sce100 <- spikeInMarkers(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                      subsetCondition = "covid",
                      subsetConditionValue = "healthy",
                      baselineVar = "platelets",
                      baselineValue = "B",
                      activatedValue = "A",
                      markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"), 
                      alpha = 1.0)
#outputFile <- "~/cytof/covid/sce_spiked_full.rds"
#saveRDS(sce, outputFile)


# ---------------------------------------------------------------
# See if it worked
# ---------------------------------------------------------------
#sce <- readRDS("~/cytof/covid/sce_spiked_full.rds")

source("functions/prep_functions.R")
sce <- transformData(sce = sce)
sce25 <- transformData(sce = sce25)
sce50 <- transformData(sce = sce50)
sce75 <- transformData(sce = sce75)
sce100 <- transformData(sce = sce100)

source("functions/cluster_functions.R")
sce <- clusterSCE(sce)
sce25 <- clusterSCE(sce25)
sce50 <- clusterSCE(sce50)
sce75 <- clusterSCE(sce75)
sce100 <- clusterSCE(sce100)

saveRDS(sce, "~/cytof/covid/sce_spiked_clustered_full.rds")
saveRDS(sce25, "~/cytof/covid/sce_spiked_clustered_25.rds")
saveRDS(sce50, "~/cytof/covid/sce_spiked_clustered_50.rds")
saveRDS(sce75, "~/cytof/covid/sce_spiked_clustered_75.rds")
saveRDS(sce100, "~/cytof/covid/sce_spiked_clustered_100.rds")


source("functions/de_functions.R")
source("functions/sceEMD.R")
allResults <- list()
i <- 1
names <- c("sce", "sce25", "sce50", "sce75", "sce100")
library(BiocParallel)
param <- MulticoreParam(workers = 18, progressbar = T)
register(param)
for(sceObj in list(sce, sce25, sce50, sce75, sce100)){
  print(i)
  results <- runDS(sceObj, 
                   ds_methods = c("diffcyt-DS-limma","diffcyt-DS-LMM","sceEMD"), 
                   clustering_to_use = "all", contrast_vars = "base_spike", 
                   design_matrix_vars = c("patient_id", "base_spike"), fixed_effects = "base_spike", 
                   random_effects = "patient_id", markers_to_test = c("type", "state"), 
                   sceEMD_condition = "base_spike", binSize = 0, nperm = 500, parallel=T)
  allResults[[ names[i] ]] <- results
  i <- i + 1
}

source("functions/venn_functions.R")
createVennDiagram(allResults[["sce25"]], DS=T, 0.05)
createVennDiagram(allResults[["sce50"]], DS=T, 0.05)
createVennDiagram(allResults[["sce75"]], DS=T, 0.05)
createVennDiagram(allResults[["sce100"]], DS=T, 0.05)
source("functions/cytoGLMM_functions.R")
makeDF <- function(model){
  df <- data.frame(
    cluster_id = factor("all"),
    marker_id = factor(summary(model)$protein_name),
    p_val <- summary(model)$pvalues_unadj,
    p_adj <- summary(model)$pvalues_adj
  )
  colnames(df) <- c("cluster_id", "marker_id", "p_val", "p_adj")
  df
}
cytoModelSCE <- runCytoGLMM(sce=sce,
            condition="base_spike",
            group="patient_id")
plot(cytoModelSCE)
allResults[["sce"]]$cytoGLMM <- makeDF(cytoModelSCE)
createVennDiagram(allResults[["sce"]], DS=T, 0.05)


cytoModelSCE25 <- runCytoGLMM(sce=sce25,
                              condition="base_spike",
                              group="patient_id")
plot(cytoModelSCE25)
allResults[["sce25"]]$cytoGLMM <- makeDF(cytoModelSCE25)
createVennDiagram(allResults[["sce25"]], DS=T, 0.05)


cytoModelSCE50 <- runCytoGLMM(sce=sce50,
                              condition="base_spike",
                              group="patient_id")
plot(cytoModelSCE50)
allResults[["sce50"]]$cytoGLMM <- makeDF(cytoModelSCE50)
createVennDiagram(allResults[["sce50"]], DS=T, 0.05)


cytoModelSCE75 <- runCytoGLMM(sce=sce75,
                              condition="base_spike",
                              group="patient_id")
plot(cytoModelSCE75)
allResults[["sce75"]]$cytoGLMM <- makeDF(cytoModelSCE75)
createVennDiagram(allResults[["sce75"]], DS=T, 0.05)


cytoModelSCE100 <- runCytoGLMM(sce=sce100,
                               condition="base_spike",
                               group="patient_id")
plot(cytoModelSCE100)
allResults[["sce100"]]$cytoGLMM <- makeDF(cytoModelSCE100)
createVennDiagram(allResults[["sce100"]], DS=T, 0.05)



