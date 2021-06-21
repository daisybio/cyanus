# SPIKE COVID PLATELETS EVALUATION

source('DEComparison/benchmarking_plots.R')

stats_table <- preparePlotData('downsampled_covid_spike')
plot_cells_vs_elapsed(stats_table)
plot_sens_vs_pre(stats_table)
plot_sens_vs_spec(stats_table)
plot_f1_vs_elapsed(stats_table)

# inspection of downsampling for CD154

library(CATALYST)
library(data.table)
library(ggplot2)
full_sce_paths <- list.files('/localscratch/quirinmanz/cytof_data/covid_spiked', pattern = '\\_full.rds$', full.names = TRUE)
names(full_sce_paths) <- sapply(strsplit(full_sce_paths, split = '_', fixed = T), `[`, 6)
sces_full <- rbindlist(lapply(full_sce_paths, function(x) {
  x <- readRDS(x)
  dt <- cbind(as.data.table(t(assay(x, 'exprs'))), as.data.frame(colData(x)))
  dt_melt <- melt(dt, measure.vars = rownames(x))
}), idcol = 'dataset')

# 'CD63', 'CD62P', 'CD107a', 
CD154_medians <- sces_full[variable %in% c('CD154')][, .(`CD154 Median Marker Expression` = median(value)), by=.(dataset, patient_id, base_spike)]
CD154_medians[, Alpha:=factor(dataset, levels=c('full', '25', '50', '75', '100'))]
cd154_median_plot <- ggplot(CD154_medians, aes(x=Alpha, y=`CD154 Median Marker Expression`, color=base_spike)) +
  geom_point(size=4) +
  facet_wrap(~patient_id, scales='free') +
  scale_fill_discrete(name = "Condition") +
  theme_bw() +
  theme(text= element_text(size=20))
ggsave('DEComparison/downsampled_covid_spike/cd154_median_plot.pdf', cd154_median_plot)


# HEATMAP

path <-  "DEComparison/downsampled_covid_spike"
trues <- c("CD62P", "CD63", "CD107a", "CD154")
keep <- 6
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
results[, alpha := tstrsplit(dataset, "_", keep=4)]

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


tmp <- data.frame(method = results$method, marker_id = results$marker_id, p_adj = results$p_adj, nr_of_cells = results$nr_of_cells)
tmp$alpha <- factor(results$alpha, levels=c('full', '25', '50', '75', '100'))
tmp$significant <- as.factor(results$p_adj < 0.05)
tmp$p_adj <- NULL
tmp$class <- tmp$marker_id %in% trues
tmp$class[tmp$class == TRUE] <- "State"
tmp$class[tmp$class == FALSE] <- "Type"
tmp <- as.data.table(tmp)

tmp$nr_of_cells[tmp$nr_of_cells == "full"] <- paste("full (4052622)")
tmp$nr_of_cells <- factor(tmp$nr_of_cells, levels=c("1000", "2000", "5000", "10000", "15000", "20000", "full (4052622)"))

eff[, nr_of_cells := tstrsplit(dataset, "_", keep=keep)]
eff$nr_of_cells[eff$nr_of_cells == "full"] <- paste("full (4052622)")
eff$nr_of_cells <- factor(eff$nr_of_cells, levels=c("1000", "2000", "5000", "10000", "15000", "20000", "full (4052622)"))
eff[, marker_id:= tstrsplit(group1, "::", keep=1)]
eff[, alpha:=tstrsplit(dataset, "_", keep=4)]
eff$alpha <- factor(eff$alpha, levels=c('full', '25', '50', '75', '100'))

marker_classes <- tmp[,c("marker_id", "class")]
marker_classes <- unique(marker_classes)

eff <- merge(eff,marker_classes, by="marker_id")

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# only one marker class
ggplot(as.data.table(tmp[class == "Type"]), aes(marker_id, method)) +
  geom_tile(aes(fill = significant), color = "white", size = 1) +
  ggtitle("") +
  xlab(label = " Type Markers") +
  ylab("Method") +
  facet_grid(alpha ~ nr_of_cells, scales = "free_x") +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = colorBlindBlack8[c(7, 3, 1)],
                    name = "Significant",
                    na.value = "transparent") +
  ggside::geom_xsidetile(
    data = eff[class == "Type"],
    aes(y = overall_group, xfill = magnitude),
    color = "white",
    size = 0.2
  ) +
  ggside::scale_xfill_manual(values = colorBlindBlack8[c(8, 5, 2, 6)],
                             name = 'Effect Size\nMagnitude',
                             na.value = "transparent")

# For one specific nr of cells
ggplot(tmp[nr_of_cells == "full (4052622)"], aes(marker_id, method)) +
  geom_tile(aes(fill = significant), color = "white", size = 1) +
  ggtitle("") +
  xlab(label = "Markers") +
  ylab("Method") +
  facet_grid(alpha ~ class, scales = "free_x") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = colorBlindBlack8[c(7, 3, 1)],
                    name = "Significant",
                    na.value = "transparent") +
  ggside::geom_xsidetile(
    data = eff[nr_of_cells == "full (4052622)"],
    aes(y = overall_group, xfill = magnitude),
    color = "white",
    size = 0.2
  ) +
  ggside::scale_xfill_manual(values = colorBlindBlack8[c(8, 5, 2, 6)],
                             name = 'Effect Size\nMagnitude',
                             na.value = "transparent")

  