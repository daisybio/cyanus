spikeInMarkers <- function(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                           subsetCondition = "covid",
                           subsetConditionValue = "healthy",
                           baselineVar = "platelets",
                           baselineValue = "B",
                           activatedValue = "A",
                           markersToSpike = c("CD63", "CD62P", "CD107a", "CD154")) {
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
    assays(sceBaseline)$counts[marker, inds_spike] <-
      assays(sceActivated)$counts[marker, indsActivated]
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
                      markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"))
#outputFile <- "~/cytof/covid/sce_spiked_full.rds"
#saveRDS(sce, outputFile)


# ---------------------------------------------------------------
# See if it worked
# ---------------------------------------------------------------
#sce <- readRDS("~/cytof/covid/sce_spiked_full.rds")

source("functions/prep_functions.R")
sce <- transformData(sce = sce)

source("functions/cluster_functions.R")
sce <- clusterSCE(sce)

source("functions/de_functions.R")

parameters <- prepDiffExp(sce = sce, 
                          contrastVars = c("base_spike"), 
                          colsFixed= c("base_spike"),
                          colsRandom = c("patient_id"),
                          method = "diffcyt-DS-LMM")
markersToTest <- c("type", "state")
is_marker <- rowData(sce)$marker_class %in% c("type", "state")
markersToTest <- (rowData(sce)$marker_class %in% markersToTest)[is_marker]

library(diffcyt)
out <- diffcyt::diffcyt(
  d_input = sce,
  formula = parameters[["formula"]],
  contrast = parameters[["contrast"]],
  analysis_type = "DS",
  method_DS = "diffcyt-DS-LMM",
  clustering_to_use = "all",
  markers_to_test = markersToTest
)
topTable <- as.data.table(diffcyt::topTable(out$res,all=TRUE,format_vals=TRUE))
CATALYST::plotDiffHeatmap(
  x = sce, 
  y = rowData(out$res),
  all = T
)

plotDiffHeatmapCustom(
  x = sce, 
  y = rowData(out$res),
  all = T,
  normalize = F
)

