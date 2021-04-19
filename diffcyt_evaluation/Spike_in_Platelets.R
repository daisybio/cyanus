library(CATALYST)
library(data.table)
sce <- readRDS("data/platelets/sce.rds")

sceBaseline <- filterSCE(sce, activated_baseline == "B")
sceActivated <- filterSCE(sce, activated_baseline == "A")

# ---------------------------------------------------------------
# Split each baseline sample into two halves: 'base' and 'spike'
# ---------------------------------------------------------------

n_cells_baseline <- ei(sceBaseline)$n_cells
start_next_sample <- c(0, n_cells_baseline)
start_next_sample <- cumsum(start_next_sample[-length(start_next_sample)])
set.seed(100)

# generate random indices
inds <- lapply(n_cells_baseline, function(n) {
  i_base <- sort(sample(seq_len(n), floor(n / 2)))
  i_spike <- setdiff(seq_len(n), i_base)
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

#how many do I have in activated?
n_cells_activated <- ei(sceActivated)$n_cells
n_cells_activated - needed_per_sample
# --> enough
start_of_next_sample_A <- c(0, n_cells_activated)
start_of_next_sample_A <- cumsum(start_of_next_sample_A[-length(start_of_next_sample_A)])

## select correct number of activated cells for each sample and generate random indices
set.seed(100)
indsActivated <- unlist(sapply(c(1:length(needed_per_sample)), function(i) {
  needed <- needed_per_sample[i]
  got <- n_cells_activated[i]
  sort(sample(seq_len(got), needed)) + start_of_next_sample_A[i]
}))

#CD63
assays(sceBaseline)$counts["CD63", inds_spike] <- assays(sceActivated)$counts["CD63", indsActivated]
#CD62P
assays(sceBaseline)$counts["CD62P", inds_spike] <- assays(sceActivated)$counts["CD62P", indsActivated]

n_cells <- as.vector(n_cells(sceBaseline))

exp_info <- data.table(
  sample_id = levels(colData(sceBaseline)$sample_id), 
  activated_baseline = "B", 
  dual_triple = rep(ei(sceBaseline)$dual_triple, each = 2),
  patient_id = tstrsplit(levels(colData(sceBaseline)$sample_id), "_", fixed = T, keep = 1)[[1]], 
  n_cells = n_cells
)

metadata(sceBaseline)$experiment_info <- exp_info




