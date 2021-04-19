#produce the SCE objects for the diffcyt dataset
#The .fcs files have already been generated, now we have to read them in as a SCE object

#function for the metadata
makeMetadata <- function(filenames){
  patientIDs <- unlist(sapply(filenames, tstrsplit, "\\_", keep=4))
  condition <- unlist(sapply(filenames, tstrsplit, "\\_", keep=5))
  
  metadata <- data.table(
    file_name = filenames,
    condition = condition,
    patient_id = patientIDs
  )
  metadata <- metadata %>% tidyr::unite("sample_id", condition:patient_id, sep = "_", remove = F)
}

produceSCE <- function(flowSet, metadata, panel, pop_names, pop_assignment){
  cn <- colnames(flowSet@frames[[metadata$file_name[1]]])
  #pretend that population is a marker
  cn[cn == "population"] <- pop_assignment
  for(elem in metadata$file_name){
    colnames(flowSet@frames[[elem]]@exprs) <- cn
  }
  colnames(flowSet) <- cn
  
  sce <- prepData(flowSet, panel = panel, md = metadata)
  mapping <- setNames(pop_names$name, pop_names$label)
  colData(sce)$cluster_id <- as.factor(mapping[assays(sce)$counts["population", ]])
  metadata(sce)$cluster_codes <- data.frame(population = as.factor(levels(colData(sce)$cluster_id)))
  
  #delete population from assays
  sce <- filterSCE(sce, !(marker_name %in% "population"))
  return(sce)
}

#Main dataset: 
files <- list.files("./benchmark_data/BCR_XL_sim/data/main/", pattern = "\\.fcs", full.names = T)
flowSetMain <- read.flowSet(files, transformation = F, truncate_max_range = F)
metadata <- fread("./benchmark_data/BCR_XL_sim/data/main/metadata.csv")
panel <- fread("./benchmark_data/BCR_XL_sim/data/main/panel.csv")
panel <- rbind(panel, data.table(fcs_colname = "population(Population185)Dd", antigen = "population", marker_class = "type"))
pop_names <-  fread("benchmark_data/BCR_XL_sim/data/population_names/population_names.csv")

sceDiffcytMain <- produceSCE(
  flowSet = flowSetMain, metadata = metadata, panel = panel, pop_names = pop_names, pop_assignment = "population(Population185)Dd"
)

saveRDS(sceDiffcytMain, "./data/diffcyt/main/sce.rds")
saveRDS(panel, "./data/diffcyt/main/panel.rds")
saveRDS(metadata, "./data/diffcyt/main/md.rds")

#Now do the same thing for the less distinct data set: 
#1. less_50pc
files <- list.files("./benchmark_data/BCR_XL_sim/data/less_distinct/less_50pc/", pattern = "\\.fcs", full.names = T)
flowSet50pc <- read.flowSet(files, transformation = F, truncate_max_range = F)
metadata50pc <- makeMetadata(list.files("./benchmark_data/BCR_XL_sim/data/less_distinct/less_50pc/", pattern = "\\.fcs"))
#use same panel
sceDiffcyt50pc <- produceSCE(
  flowSet = flowSet50pc, metadata = metadata50pc, panel = panel, pop_names = pop_names, pop_assignment = "population(Population185)Dd"
)
saveRDS(sceDiffcyt50pc, "./data/diffcyt/less_50pc/sce.rds")
saveRDS(panel, "./data/diffcyt/less_50pc/panel.rds")
saveRDS(metadata50pc, "./data/diffcyt/less_50pc/md.rds")

#2. less_75pc
files <- list.files("./benchmark_data/BCR_XL_sim/data/less_distinct/less_75pc/", pattern = "\\.fcs", full.names = T)
flowSet75pc <- read.flowSet(files, transformation = F, truncate_max_range = F)
metadata75pc <- makeMetadata(list.files("./benchmark_data/BCR_XL_sim/data/less_distinct/less_75pc/", pattern = "\\.fcs"))
#use same panel
sceDiffcyt75pc <- produceSCE(
  flowSet = flowSet75pc, metadata = metadata75pc, panel = panel, pop_names = pop_names, pop_assignment = "population(Population185)Dd"
)
saveRDS(sceDiffcyt75pc, "./data/diffcyt/less_75pc/sce.rds")
saveRDS(panel, "./data/diffcyt/less_75pc/panel.rds")
saveRDS(metadata75pc, "./data/diffcyt/less_75pc/md.rds")










  

