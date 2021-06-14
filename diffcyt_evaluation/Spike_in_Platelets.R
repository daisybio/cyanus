spikeInMarkers <- function(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                           subsetCondition = "covid",
                           subsetConditionValue = "healthy",
                           baselineVar = "platelets",
                           baselineValue = "B",
                           activatedValue = "A",
                           markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"),
                           alpha = NULL) {
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
  
  if(is.null(alpha)){
    sce0 <- sceBaseline
    sce25 <- sceBaseline
    sce50 <- sceBaseline
    sce75 <- sceBaseline
    sce100 <- sceBaseline
    sceList <- list(sce0, sce25, sce50, sce75, sce100)
    names(sceList) <- c("0", "0.25", "0.5", "0.75", "1.0")
  }
  
  
  for (marker in markersToSpike) {
    message(paste("Replacing marker", marker, "..."))
    if(is.null(alpha)){
      for (alpha2 in names(sceList)){
        assays(sceList[[alpha2]])$counts[marker, inds_spike] <-
          5 * sinh(
            asinh(assays(sceActivated)$counts[marker, indsActivated]/5) - 
              as.numeric(alpha2) * (
                asinh(assays(sceActivated)$counts[marker, indsActivated] / 5) - 
                  asinh(assays(sceList[[alpha2]])$counts[marker, inds_spike] / 5)
              )
          )
      }
    }else if(alpha == 0){
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
  if(is.null(alpha)){
    sceList <- lapply(sceList, function(x){
      metadata(x)$experiment_info <- exp_info
      return(x)
    })
    return(sceList)
  }else{
    metadata(sceBaseline)$experiment_info <- exp_info
    return(sceBaseline)
  }
}

sceList <- spikeInMarkers(pathToSCE = "/nfs/home/students/l.arend/data/covid/sce_untransformed.rds",
                          subsetCondition = "covid",
                          subsetConditionValue = "healthy",
                          baselineVar = "platelets",
                          baselineValue = "B",
                          activatedValue = "A",
                          markersToSpike = c("CD63", "CD62P", "CD107a", "CD154"), 
                          alpha = NULL)

sce <- sceList[["0"]]
sce25 <- sceList[["0.25"]]
sce50 <- sceList[["0.5"]]
sce75 <- sceList[["0.75"]]
sce100 <- sceList[["1.0"]]
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

######tests
sce <- readRDS("~/cytof/covid/sce_spiked_clustered_full.rds")
sce25 <- readRDS("~/cytof/covid/sce_spiked_clustered_25.rds")
sce50 <- readRDS("~/cytof/covid/sce_spiked_clustered_50.rds")
sce75 <- readRDS("~/cytof/covid/sce_spiked_clustered_75.rds")
sce100 <- readRDS("~/cytof/covid/sce_spiked_clustered_100.rds")

library(data.table)
exprsDT <- as.data.table(t(assays(sce)$counts))
coldataDT <- as.data.table(colData(sce))[, c("patient_id", "base_spike")]
exprsDT <- cbind(exprsDT, coldataDT)
CD63_means <- exprsDT[, mean(CD63), by = .(patient_id, base_spike)]
colnames(CD63_means) <- c("patient_id", "base_spike", "mean_full")
CD62P_means <- exprsDT[, mean(CD62P), by = .(patient_id, base_spike)]
colnames(CD62P_means) <- c("patient_id", "base_spike", "mean_full")
CD107a_means <- exprsDT[, mean(CD107a), by = .(patient_id, base_spike)]
colnames(CD107a_means) <- c("patient_id", "base_spike", "mean_full")
CD154_means <- exprsDT[, mean(CD154), by = .(patient_id, base_spike)] 
colnames(CD154_means) <- c("patient_id", "base_spike", "mean_full")
CD63_means[, sd_full := exprsDT[, sd(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, sd_full := exprsDT[, sd(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, sd_full := exprsDT[, sd(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, sd_full := exprsDT[, sd(CD154), by = .(patient_id, base_spike)]$V1]


exprsDT <- as.data.table(t(assays(sce25)$counts))
coldataDT <- as.data.table(colData(sce25))[, c("patient_id", "base_spike")]
exprsDT <- cbind(exprsDT, coldataDT)
CD63_means[, mean_25 := exprsDT[, mean(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, mean_25 := exprsDT[, mean(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, mean_25 := exprsDT[, mean(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, mean_25 := exprsDT[, mean(CD154), by = .(patient_id, base_spike)]$V1]
CD63_means[, sd_25 := exprsDT[, sd(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, sd_25 := exprsDT[, sd(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, sd_25 := exprsDT[, sd(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, sd_25 := exprsDT[, sd(CD154), by = .(patient_id, base_spike)]$V1]

exprsDT <- as.data.table(t(assays(sce50)$counts))
coldataDT <- as.data.table(colData(sce50))[, c("patient_id", "base_spike")]
exprsDT <- cbind(exprsDT, coldataDT)
CD63_means[, mean_50 := exprsDT[, mean(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, mean_50 := exprsDT[, mean(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, mean_50 := exprsDT[, mean(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, mean_50 := exprsDT[, mean(CD154), by = .(patient_id, base_spike)]$V1]
CD63_means[, sd_50 := exprsDT[, sd(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, sd_50 := exprsDT[, sd(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, sd_50 := exprsDT[, sd(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, sd_50 := exprsDT[, sd(CD154), by = .(patient_id, base_spike)]$V1]

exprsDT <- as.data.table(t(assays(sce75)$counts))
coldataDT <- as.data.table(colData(sce75))[, c("patient_id", "base_spike")]
exprsDT <- cbind(exprsDT, coldataDT)
CD63_means[, mean_75 := exprsDT[, mean(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, mean_75 := exprsDT[, mean(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, mean_75 := exprsDT[, mean(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, mean_75 := exprsDT[, mean(CD154), by = .(patient_id, base_spike)]$V1]
CD63_means[, sd_75 := exprsDT[, sd(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, sd_75 := exprsDT[, sd(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, sd_75 := exprsDT[, sd(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, sd_75 := exprsDT[, sd(CD154), by = .(patient_id, base_spike)]$V1]

exprsDT <- as.data.table(t(assays(sce100)$counts))
coldataDT <- as.data.table(colData(sce100))[, c("patient_id", "base_spike")]
exprsDT <- cbind(exprsDT, coldataDT)
CD63_means[, mean_100 := exprsDT[, mean(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, mean_100 := exprsDT[, mean(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, mean_100 := exprsDT[, mean(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, mean_100 := exprsDT[, mean(CD154), by = .(patient_id, base_spike)]$V1]
CD63_means[, sd_100 := exprsDT[, sd(CD63), by = .(patient_id, base_spike)]$V1]
CD62P_means[, sd_100 := exprsDT[, sd(CD62P), by = .(patient_id, base_spike)]$V1]
CD107a_means[, sd_100 := exprsDT[, sd(CD107a), by = .(patient_id, base_spike)]$V1]
CD154_means[, sd_100 := exprsDT[, sd(CD154), by = .(patient_id, base_spike)]$V1]

all_stats <- rbindlist(list("CD63" = CD63_means, 
                            "CD62P" = CD62P_means, 
                            "CD107a" = CD107a_means, 
                            "CD154" = CD154_means), idcol = "marker")

all_stats <- melt(all_stats, id.vars = c("marker", "patient_id", "base_spike"), 
                  value.name = "value", variable.name="variable")
all_stats[, model := tstrsplit(variable, "_", keep=2)]
all_stats[, model := factor(model, levels = c("full", "25", "50", "75", "100"))]
all_stats[, mean_sd := tstrsplit(variable, "_", keep=1)]
all_stats[, mean_sd := as.factor(mean_sd)]
all_stats[, variable := NULL]
library(ggplot2)
ggplot(all_stats, aes(x = model, y = value, fill = base_spike))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(mean_sd ~ marker, scales = "free")+
  scale_fill_manual(values = c("indianred1", "lightblue1"))

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



