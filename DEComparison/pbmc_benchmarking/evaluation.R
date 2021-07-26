# PBMC SIMULATION EVALUATION

library(CATALYST)
library(ggplot2)
source('DEComparison/benchmarking_plots.R')

# PLOT SIZE OF DATASET
sce_pbmc <- readRDS('/nfs/home/students/jbernett/cytof/cytof/data/cytof_workflow_SCE.rds')
pbmc_coldata <- as.data.frame(colData(sce_pbmc))
pbmc_coldata$merging1 <- cluster_ids(sce_pbmc, 'merging1')
pbmc_cell_counts <- ggplot(pbmc_coldata, aes(x=sample_id)) +
  geom_bar() + 
  facet_wrap(~merging1) + 
  labs(y = 'Cell Count', x = 'Sample ID') +
  theme_bw() + 
  theme(text = element_text(size = 20),  axis.text.x = element_text(angle = 90, vjust =.5)) 
ggsave(pbmc_cell_counts, filename = 'DEComparison/pbmc_benchmarking/pbmc_cell_counts.png', width = 14, height = 7)


# HEATMAP PAPER
path <- "DEComparison/pbmc_benchmarking"
result_rds <-
  list.files(path = path,
             pattern = "\\.rds$",
             full.names = T)
results <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["results"]], simplify = FALSE),
    idcol = "dataset")
results[, dataset := basename(dataset)]
tmp <- data.frame(method = results$method, marker_id = results$marker_id, p_adj = results$p_adj)
tmp$cluster_id <- results$cluster_id
times <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["times"]], simplify = FALSE),
    idcol = "dataset")
times[, dataset := basename(dataset)]
eff <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["eff"]], simplify = FALSE),
    idcol = "dataset")
eff[, dataset := basename(dataset)]

# plot results
tmp$significant <- results$p_adj < 0.05
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

eff$marker_id <- sapply(strsplit(eff$group2,'::'), "[", 1)
tmp$marker_id[tmp$marker_id == "HLADR"] <- "HLA_DR"

# rename and order methods
tmp[tmp == "t_test"] <- "t-test"
tmp[tmp == "wilcoxon_median"] <- "Wilcoxon test"
tmp[tmp == "kruskal_median"] <- "Kruskal-Wallis test"
tmp$method[tmp$method == "sceEMD"] <- "CyEMD"
tmp$method <- factor(tmp$method, levels=rev(c("diffcyt-DS-limma", "diffcyt-DS-LMM", "t-test", "Wilcoxon test","Kruskal-Wallis test", "CytoGLM","CytoGLMM", "logRegression", "ZAGA", "BEZI", "CyEMD")))

tmp$significant <- as.character(tmp$significant)
tmp$significant[tmp$significant == TRUE] <- "Yes"
tmp$significant[tmp$significant == FALSE] <- "No"

eff$magnitude <- as.character(eff$magnitude)
eff$magnitude[eff$magnitude == "small"] <- "Small"
eff$magnitude[eff$magnitude == "negligible"] <- "Negligible"
eff$magnitude[eff$magnitude == "large"] <- "Large"
eff$magnitude[eff$magnitude == "moderate"] <- "Moderate"
eff$magnitude <- factor(eff$magnitude, levels=c("Negligible", "Small", "Moderate", "Large"))

ggplot(tmp, aes(marker_id, method)) + 
  geom_tile(aes(fill=significant), color="white", size=1) + 
  ggtitle("") + xlab(label="Marker") + ylab('Method') + 
  facet_wrap(~cluster_id, scales = "free_x") + 
  theme(text = element_text(size = 18),  axis.text.x = element_text(angle = 90, vjust=.5))+
  scale_fill_manual(values = colorBlindBlack8[c(7,3,1)], name="Significant") +
  ggside::geom_xsidetile(data=eff, aes(y=overall_group, xfill=magnitude),  color="white", size=0.2) +
  ggside::scale_xfill_manual(values=c(colorBlindBlack8[c(8,5,4,6)]), name='Cohen’s d\nEffect size\nMagnitude')
