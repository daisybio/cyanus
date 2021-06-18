# CYTOGLMM SIMULATION EVALUATION

source('DEComparison/benchmarking_plots.R')

stats_table <- preparePlotData('simulatedCytoGLMM')
plot_cells_vs_elapsed(stats_table)
plot_sens_vs_pre(stats_table)
plot_sens_vs_spec(stats_table)
plot_f1_vs_elapsed(stats_table)
plotHeatmaps('simulatedCytoGLMM')


# Plot Heatmap for Paper
path <-  "DEComparison/simulatedCytoGLMM"
trues <- c("m01", "m02", "m03", "m04", "m05")
keep <- 3

result_rds <-
  list.files(path = path,
             pattern = "\\.rds$",
             full.names = T)
results <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["results"]], simplify = FALSE),
    idcol = "dataset")
results[, dataset := basename(dataset)]
results[, nr_of_cells := tstrsplit(dataset, "_", keep=keep)]
tmp <- data.frame(method = results$method, marker_id = results$marker_id, p_adj = results$p_adj, nr_of_cells = results$nr_of_cells)
tmp$significant <- results$p_adj < 0.05
tmp$nr_of_cells <- results$nr_of_cells
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
tmp$class <- tmp$marker_id %in% trues
tmp$class[tmp$class == TRUE] <- "Differentially Expressed"
tmp$class[tmp$class == FALSE] <- "Not Differentially Expressed"


eff <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["eff"]], simplify = FALSE),
    idcol = "dataset")
eff[, dataset := basename(dataset)]
eff$nr_of_cells <- sapply(strsplit(eff$dataset,'_'), "[", 3)
eff$marker_id <- sapply(strsplit(eff$group1,'::'), "[", 1)



colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# for paper only 200000 cells
#tmp <- tmp[tmp$nr_of_cells=="200000",]
marker_classes <- tmp[,c("marker_id", "class")]
marker_classes <- unique(marker_classes)
#eff <- eff[eff$nr_of_cells=="200000",]
eff <- eff[,c("overall_group", "magnitude", "marker_id", "nr_of_cells")]
eff <- merge(eff,marker_classes, by="marker_id")

ggplot(tmp, aes(marker_id, method)) + 
  geom_tile(aes(fill=significant), color="white", size=1) + 
  ggtitle("") + 
  xlab(label="Marker") + 
  ylab("Method") +
  #facet_wrap(~class, scales = "free_x") + 
  facet_grid(nr_of_cells~class, scales = "free_x") +
  theme(text = element_text(size = 14),  axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values = colorBlindBlack8[c(7,3)], na.value="transparent", name="Significant") + 
  ggside::geom_xsidetile(data=eff, aes(y=overall_group, xfill=magnitude), color="white", size=0.2) + 
  ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,2,6)], name='Effect Size\nMagnitude', na.value="transparent", drop=FALSE)



