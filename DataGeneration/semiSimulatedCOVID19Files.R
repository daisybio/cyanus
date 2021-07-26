source("functions/prep_functions.R")
source("functions/cluster_functions.R")


spikeInMarkers <- function(pathToSCE = "DataGeneration/covid/sce_untransformed.rds",
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

#### Read in original data

panel <-
  read_excel(
    "DataGeneration/covid/panel_umap.xlsx"
  )
md <- 
  read_excel(
    "DataGeneration/covid/meta_11vs8_batch.xlsx"
  )
#download the data from FlowRepository
exp <-
  list.files(
    "DataGeneration/covid/covid_platelets",
    pattern = "\\.fcs$",
    full.names = T
  )
sce_covid <-
  prepData(
    exp,
    panel,
    md,
    transform = FALSE,
    md_cols = list(
      file = "file_name",
      id = "sample_id",
      factors = c("covid","patient_id", "platelets", "batch", "acquisition_date")
    )
  )

sce_covid <- clusterSCE(sce_covid)
saveRDS(sce_covid, "DataGeneration/covid/sce_untransformed.rds")

#### Spike in CD63, CD62P, CD107a and CD154 expressions

sceList <- spikeInMarkers(pathToSCE = "DataGeneration/covid/sce_untransformed.rds",
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

sce <- transformData(sce = sce)
sce25 <- transformData(sce = sce25)
sce50 <- transformData(sce = sce50)
sce75 <- transformData(sce = sce75)
sce100 <- transformData(sce = sce100)

sce <- clusterSCE(sce)
sce25 <- clusterSCE(sce25)
sce50 <- clusterSCE(sce50)
sce75 <- clusterSCE(sce75)
sce100 <- clusterSCE(sce100)

saveRDS(sce, "DataGeneration/covid/sce_spiked_clustered_full_ds_full.rds")
saveRDS(sce25, "DataGeneration/covid/sce_spiked_clustered_25_ds_full.rds")
saveRDS(sce50, "DataGeneration/covid/sce_spiked_clustered_50_ds_full.rds")
saveRDS(sce75, "DataGeneration/covid/sce_spiked_clustered_75_ds_full.rds")
saveRDS(sce100, "DataGeneration/covid/sce_spiked_clustered_100_ds_full.rds")

# ---------------------------------------------------------------
# Perform downsampling
# ---------------------------------------------------------------

scefull <- readRDS("DataGeneration/covid/sce_spiked_clustered_full_ds_full.rds")
sce25 <- readRDS("DataGeneration/covid/sce_spiked_clustered_25_ds_full.rds")
sce50 <- readRDS("DataGeneration/covid/sce_spiked_clustered_50_ds_full.rds")
sce75 <- readRDS("DataGeneration/covid/sce_spiked_clustered_75_ds_full.rds")
sce100 <- readRDS("DataGeneration/covid/sce_spiked_clustered_100_ds_full.rds")

source("functions/prep_functions.R")
for(scetmp in c("scefull", "sce25", "sce50", "sce75", "sce100")){
  last_sce <- get(scetmp)
  for (n in rev(c(1000, 2000, 5000, 10000, 15000, 20000))) {
    sampling <- tstrsplit(scetmp, "sce", keep=2)[[1]]
    print(CATALYST::ei(last_sce))
    downsampled_sce <- downSampleSCE(sce=last_sce, cells = n, per_sample = T, seed = 1234)
    saveRDS(downsampled_sce, paste0("DataGeneration/covid/sce_spiked_clustered_", sampling, "_ds_", n, ".rds"))
    last_sce <- downsampled_sce
  }
}

### plot spiked-in data 


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
spike_in_sces <- lapply(list.files("DataGeneration/covid/downsampled_files", full.names = T), readRDS)
names(spike_in_sces) <- tstrsplit(list.files("DataGeneration/covid/downsampled_files"), ".rds", keep=1)[[1]]
countsDT <- lapply(spike_in_sces, function(x){as.data.table(t(assays(x)$exprs))})
countsDT <- rbindlist(countsDT, idcol = "filename")
coldataDT <- lapply(spike_in_sces, function(x){as.data.table(colData(x))[, c("patient_id", "base_spike")]})
coldataDT <- rbindlist(coldataDT, idcol = "filename")
countsDT <- cbind(coldataDT[, c("patient_id", "base_spike")], countsDT)
countsDT[, c("alpha", "n_cells") := tstrsplit(filename, "_", keep = c(4,6))]

statemarkerDT <- countsDT[, c("CD62P", "CD63", "CD154", "CD107a", "patient_id", "base_spike", "alpha", "n_cells")]
statemarkerDT <- data.table::melt(statemarkerDT, id.vars = c("patient_id","base_spike", "alpha", "n_cells"), variable.name = "marker", value.name = "expression")
meanDT <- statemarkerDT[, mean(expression), by = c("marker", "patient_id", "base_spike", "alpha", "n_cells")]
colnames(meanDT) <- c("marker", "patient", "base_spike", "alpha", "n_cells", "mean")
meanDT[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
meanDT[, base_spike := factor(base_spike, levels = c("spike", "base"))]
meanDT[, n_cells := factor(n_cells, levels = c(1000, 2000, 5000, 10000, 15000, 20000))]

spike <- ggplot(meanDT, aes(x=alpha, y = mean, fill = base_spike))+
  geom_boxplot()+
  facet_wrap( ~ marker, scales = "free")+
  scale_fill_manual(values = colorBlindBlack8[c(3,8)], name = "Condition", labels = c("Spike", "Base"))+
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=16),
    axis.text = element_text(color = "black", size=14), 
    axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16)) +
  labs(x = expression(alpha*" x 100"),y = "Mean Expression")

