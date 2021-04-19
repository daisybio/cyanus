library(CATALYST)
library(data.table)
sce <- readRDS("/nfs/home/students/l.arend/data/covid/sce_untransformed.rds")
sce <- filterSCE(sce, covid == "healthy")
sceBaseline <- filterSCE(sce, platelets == "B")
sceActivated <- filterSCE(sce, platelets == "A")
remove(sce)

# ---------------------------------------------------------------
# Split each baseline sample into two halves: 'base' and 'spike'
# ---------------------------------------------------------------

n_cells_baseline <- ei(sceBaseline)$n_cells
start_next_sample <- c(0, n_cells_baseline)
start_next_sample <- cumsum(start_next_sample[-length(start_next_sample)])

#how many do I have in activated?
n_cells_activated <- ei(sceActivated)$n_cells
n_cells_activated - floor(n_cells_baseline / 2 )
# --> enough ? no -> take all there is from the activated condition for patient CVD003A and CVD013A

set.seed(100)
# generate random indices
inds <- lapply(c(1:length(n_cells_baseline)), function(i) {
  if(n_cells_activated[i] - floor(n_cells_baseline[i] / 2 ) < 0){
    n <- n_cells_activated[i]
    i_spike <- sort(sample( seq_len(n_cells_baseline[i]), n) )
    i_base <- setdiff(seq_len(n_cells_baseline[i]), i_spike)
  }else{
    n <- n_cells_baseline[i]
    i_base <- sort(sample(seq_len(n), floor(n / 2)))
    i_spike <- setdiff(seq_len(n), i_base)
  }
  list(base = i_base, spike = i_spike)
})

inds2 <- lapply(seq(1:length(inds)), function(x){
  base_inds <- inds[[x]][["base"]] + start_next_sample[x]
  spike_inds <- inds[[x]][["spike"]] + start_next_sample[x]
  list(base = base_inds, spike = spike_inds)
})

inds_base <- unlist(lapply(inds2, function(l) l[[1]]))
inds_spike <- unlist(lapply(inds2, function(l) l[[2]]))

conditionCol <- data.table(index = c(1:sum(ei(sceBaseline)$n_cells)), sample_id = colData(sceBaseline)$sample_id)
conditionCol <- conditionCol[ , base_spike := ifelse(index %in% inds_base, "base", "spike")]
conditionCol <- conditionCol[ , sample_id := ifelse(index %in% inds_base, paste0(sample_id, "_base"), paste0(sample_id, "_spike"))]
colData(sceBaseline)$base_spike <- as.factor(conditionCol$base_spike)
colData(sceBaseline)$sample_id <- as.factor(conditionCol$sample_id)
# -------------------------------------------------------------------------------------
# Replace cells in 'spike' samples with cells from stimulated condition but 
# only CD63 and CD62P measurements
# -------------------------------------------------------------------------------------

# generate random indices from activated

#how many are needed per sample?
needed_per_sample <- sapply(inds2, function(x){
  return(length(x[["spike"]]))
})

start_of_next_sample_A <- c(0, n_cells_activated)
start_of_next_sample_A <- cumsum(start_of_next_sample_A[-length(start_of_next_sample_A)])

## select correct number of activated cells for each sample and generate random indices
set.seed(100)
indsActivated <- unlist(sapply(c(1:length(needed_per_sample)), function(i) {
  needed <- needed_per_sample[i]
  got <- n_cells_activated[i]
  if((got - needed) < 0){
    needed <- got
  }
  sort(sample(seq_len(got), needed)) + start_of_next_sample_A[i]
}))

#CD63
assays(sceBaseline)$counts["CD63", inds_spike] <- assays(sceActivated)$counts["CD63", indsActivated]
#CD62P
assays(sceBaseline)$counts["CD62P", inds_spike] <- assays(sceActivated)$counts["CD62P", indsActivated]
#CD107a
assays(sceBaseline)$counts["CD107a", inds_spike] <- assays(sceActivated)$counts["CD107a", indsActivated]
#CD154
assays(sceBaseline)$counts["CD154", inds_spike] <- assays(sceActivated)$counts["CD154", indsActivated]

n_cells <- as.vector(n_cells(sceBaseline))

exp_info <- data.table(
  sample_id = levels(colData(sceBaseline)$sample_id), 
  covid = "healthy",
  patient_id = tstrsplit(levels(colData(sceBaseline)$sample_id), "B_", fixed = T, keep = 1)[[1]], 
  platelets = "B", 
  base_spike = tstrsplit(levels(colData(sceBaseline)$sample_id), "B_", fixed = T, keep = 2)[[1]],
  n_cells = n_cells
)

metadata(sceBaseline)$experiment_info <- exp_info
saveRDS(sceBaseline, "~/cytof/covid/sce_spiked_full.rds")

# ---------------------------------------------------------------
# See if it worked
# ---------------------------------------------------------------
sce <- readRDS("~/cytof/covid/sce_spiked_full.rds")

source("functions/prep_functions.R")
sce <- transformData(sce = sce)

source("functions/cluster_functions.R")
sce <- clusterSCE(sce)

source("functions/de_functions.R")

parameters <- prepDiffExp(sce = sce, 
                          contrastVars = c("base_spike"), 
                          colsDesign = c("base_spike", "patient_id"), 
                          method = "diffcyt-DS-limma")
markersToTest <- c("type", "state")
is_marker <- rowData(sce)$marker_class %in% c("type", "state")
markersToTest <- (rowData(sce)$marker_class %in% markersToTest)[is_marker]

library(diffcyt)
out <- diffcyt::diffcyt(
  d_input = sce,
  design = parameters[["design"]],
  contrast = parameters[["contrast"]],
  analysis_type = "DS",
  method_DS = "diffcyt-DS-limma",
  clustering_to_use = "all",
  markers_to_test = markersToTest
)
topTable <- as.data.table(diffcyt::topTable(out$res,all=TRUE,format_vals=TRUE))
CATALYST::plotDiffHeatmap(
  x = sce, 
  y = rowData(out$res),
  all = T
)
