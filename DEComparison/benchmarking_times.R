# BECHNMARKING TIMES

getTimes <- function(path){
  result_rds <-
    list.files(path = path,
               pattern = "\\.rds$",
               full.names = T)
  times <-
    data.table::rbindlist(sapply(result_rds, function(rds)
      readRDS(rds)[["times"]], simplify = FALSE),
      idcol = "dataset")
  times[, dataset := basename(dataset)]
  return(times[,c("dataset", "method", "elapsed")])
}

# Simulated CytoGLMM
nr_of_samples<- 22
cytoglmm_times <- getTimes("DEComparison/simulatedCytoGLMM")
cytoglmm_times$nr_of_cells <- sapply(strsplit(cytoglmm_times$dataset,'_'), "[", 3)
cytoglmm_times$nr_of_cells <- as.numeric(cytoglmm_times$nr_of_cells)
cytoglmm_times$nr_of_cells <- cytoglmm_times$nr_of_cells * nr_of_samples
cytoglmm_times$dataset <- "simulated_CytoGLMM"
cytoglmm_times$random_effect <- "yes"

# Dual Platelets
sce_dual <- readRDS("~/data/platelets_dual/sce_dual.rds")
nr_of_cells <- sum(ei(sce_dual)$n_cells)

path_no_r <- "/nfs/home/students/l.arend/cytof/DEComparison/dual_platelets_no_random/sce_dual_res_timed.rds"
path_with_r <- "/nfs/home/students/l.arend/cytof/DEComparison/dual_platelets/sce_dual_res_timed.rds"   #with random effects

dual_random_times <- readRDS(path_with_r)[["times"]]
dual_random_times <- dual_random_times[,c("method", "elapsed")]
dual_random_times$dataset <- "dual_platelets_random"
dual_random_times$random_effect <- "yes"
dual_random_times$nr_of_cells <- nr_of_cells


dual_no_times <- readRDS(path_no_r)[["times"]]
dual_no_times <- dual_no_times[,c("method", "elapsed")]
dual_no_times$dataset <- "dual_platelets"
dual_no_times$random_effect <- "no"
dual_no_times$nr_of_cells <- nr_of_cells


# PBMC Data
sce_pbmc <- readRDS("~/data/pbmc_all/bigPBMC_SCE.rds")
nr_of_cells <- sum(ei(sce_pbmc)$n_cells)

path_pbmc <- "~/cytof/DEComparison/pbmc_benchmarking/cytof_workflow_SCE_res_timed.rds"
pbmc_times <- readRDS(path_pbmc)[["times"]]
pbmc_times <- pbmc_times[,c("method", "elapsed")]
pbmc_times$dataset <- "pbmc"
pbmc_times$random_effect <- "yes"
pbmc_times$nr_of_cells <- nr_of_cells


# Simulated Covid-Spike
sce_covid_spike <- readRDS("/localscratch/quirinmanz/cytof_data/covid_spiked/sce_spiked_clustered_full_ds_full.rds")
nr_cells_full <- sum(ei(sce_covid_spike)$n_cells)

nr_of_samples<- 22
covid_times <- getTimes("DEComparison/downsampled_covid_spike")
covid_times$nr_of_cells <- sapply(strsplit(covid_times$dataset,'_'), "[", 6)
covid_times$nr_of_cells[covid_times$nr_of_cells == "full" ] <- nr_cells_full
covid_times$nr_of_cells <- as.numeric(covid_times$nr_of_cells)
covid_times$nr_of_cells[covid_times$nr_of_cells != nr_cells_full] <- covid_times$nr_of_cells[covid_times$nr_of_cells != nr_cells_full] * nr_of_samples
covid_times$dataset <- "covid_spiked"
covid_times$random_effect <- "yes"
covid_times$nr_of_cells <- as.factor(covid_times$nr_of_cells)

tmp <- covid_times[nr_of_cells == nr_cells_full]
tmp$elapsed <- tmp$elapsed / 60
ggplot(tmp, aes(x=reorder(method, elapsed, FUN = median), y=elapsed)) + geom_boxplot() + xlab("Method") + ylab("Elapsed (in min.)") + theme_bw() + theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))


# Plot
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

times <- rbind(dual_random_times, dual_no_times, cytoglmm_times, pbmc_times, covid_times)
times$size <- times$nr_of_cells<1000000
times$size[times$size==TRUE] <- "Small Dataset"
times$size[times$size==FALSE] <- "Big Dataset"
times$size <- factor(times$size, levels=c("Small Dataset", "Big Dataset"))
times$elapsed <- times$elapsed / 60

ggplot(times,
               aes(
                 x = nr_of_cells,
                 y = elapsed,
                 color = method,
                 shape = dataset,
  )) + 
  scale_shape_manual(values=1:times[, data.table::uniqueN(method)], name="Dataset") +
  facet_wrap(~size, scales="free") +
  geom_point(size = 4, stroke=1.5) + xlab("Number of Cells") + ylab("Elapsed (in min.)") +
  theme_bw() + theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = safe_colorblind_palette, name="Method")


        