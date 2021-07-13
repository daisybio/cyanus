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
CD154_medians$dataset[CD154_medians$dataset == "full"] <- "0"
CD154_medians$dataset[CD154_medians$dataset == "full"] <- "0"
CD154_medians$dataset[CD154_medians$dataset == "25"] <- "0.25"
CD154_medians$dataset[CD154_medians$dataset == "50"] <- "0.5"
CD154_medians$dataset[CD154_medians$dataset == "75"] <- "0.75"
CD154_medians$dataset[CD154_medians$dataset == "100"] <- "1"
CD154_medians[, Alpha:=factor(dataset, levels=c('0', '0.25', '0.5', '0.75', '1'))]
CD154_medians$base_spike <- factor(CD154_medians$base_spike, levels=c("spike", "base"))
cd154_median_plot <- ggplot(CD154_medians, aes(x=Alpha, y=`CD154 Median Marker Expression`, color=base_spike)) +
  geom_point(size=4) +
  facet_wrap(~patient_id, scales='free') +
  scale_color_manual(values =colorBlindBlack8[c(3,8)], name="Condition", label=c("Spike", "Base")) +
  theme_bw() + labs(x = expression(alpha)) + 
  theme(text= element_text(size=20))
ggsave('DEComparison/downsampled_covid_spike/cd154_median_plot.pdf', cd154_median_plot)


# TABLE
stats_table <- preparePlotData("downsampled_covid_spike")
stats_table <- stats_table[,c("method", "nr_of_cells", "alpha", "sensitivity", "specificity", "precision", "F1")]

stats_table <- stats_table %>%  group_by(.dots = c("method", "nr_of_cells")) %>% summarise(sensitivity_mean = mean(sensitivity, na.rm = TRUE), 
                                                                                           sensitivity_sd = sd(sensitivity, na.rm = TRUE), 
                                                                                           specificity_mean = mean(specificity, na.rm = TRUE), 
                                                                                           specificity_sd = sd(specificity, na.rm = TRUE),
                                                                                           precision_mean = mean(precision, na.rm = TRUE),
                                                                                           precision_sd = sd(precision, na.rm = TRUE),
                                                                                           F1_mean = mean(F1, na.rm = TRUE),
                                                                                           F1_sd = sd(F1, na.rm = TRUE))


is.num <- sapply(stats_table, is.numeric)
stats_table[is.num] <- lapply(stats_table[is.num], round, 2)
stats_table$sensitivity <- paste(stats_table$sensitivity_mean, "+/-", stats_table$sensitivity_sd)
stats_table$specificity <- paste(stats_table$specificity_mean, "+/-", stats_table$specificity_sd)
stats_table$precision <- paste(stats_table$precision_mean, "+/-", stats_table$precision_sd)
stats_table$F1 <- paste(stats_table$F1_mean, "+/-", stats_table$F1_sd)

stats_table <- stats_table[,c("method","nr_of_cells","sensitivity", "specificity", "precision", "F1")]

f1 <- stats_table[,c("method", "nr_of_cells", "F1")]
f1 <- dcast(f1, method  ~ nr_of_cells, value.var="F1")
print(xtable(f1, type = "latex"))

spec <- stats_table[,c("method", "nr_of_cells", "specificity")]
spec <- dcast(spec, method  ~ nr_of_cells, value.var="specificity")
print(xtable(spec, type = "latex"))

sens <- stats_table[,c("method", "nr_of_cells", "sensitivity")]
sens <- dcast(sens, method  ~ nr_of_cells, value.var="sensitivity")
print(xtable(sens, type = "latex"))

prec <- stats_table[,c("method", "nr_of_cells", "precision")]
prec <- dcast(prec, method  ~ nr_of_cells, value.var="precision")
print(xtable(prec, type = "latex"))




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

results$marker_id <- as.character(results$marker_id)
results$marker_id[results$marker_id == "PAC1"] <- "GPIIbIIIa"
results$marker_id <- as.factor(results$marker_id)

tmp <- data.frame(method = results$method, marker_id = results$marker_id, p_adj = results$p_adj, nr_of_cells = results$nr_of_cells)
results$alpha[results$alpha == "full"] <- "0"
results$alpha[results$alpha == "25"] <- "0.25"
results$alpha[results$alpha == "50"] <- "0.5"
results$alpha[results$alpha == "75"] <- "0.75"
results$alpha[results$alpha == "100"] <- "1"
tmp$alpha <- factor(results$alpha, levels=c('0', '0.25', '0.5', '0.75', '1'))
tmp$significant <- as.factor(results$p_adj < 0.05)
tmp$p_adj <- NULL
tmp$class <- tmp$marker_id %in% trues
tmp$class[tmp$class == TRUE] <- "State"
tmp$class[tmp$class == FALSE] <- "Type"
tmp <- as.data.table(tmp)

tmp$nr_of_cells[tmp$nr_of_cells == "full"] <- paste("4052622 (full)")
tmp$nr_of_cells <- factor(tmp$nr_of_cells, levels=c("1000", "2000", "5000", "10000", "15000", "20000", "4052622 (full)"))

eff[, nr_of_cells := tstrsplit(dataset, "_", keep=keep)]
eff$nr_of_cells[eff$nr_of_cells == "full"] <- paste("4052622 (full)")
eff$nr_of_cells <- factor(eff$nr_of_cells, levels=c("1000", "2000", "5000", "10000", "15000", "20000", "4052622 (full)"))
eff[, marker_id:= tstrsplit(group1, "::", keep=1)]
eff$marker_id[eff$marker_id == "PAC1"] <- "GPIIbIIIa"
eff[, alpha:=tstrsplit(dataset, "_", keep=4)]
eff$alpha[eff$alpha == "full"] <- "0"
eff$alpha[eff$alpha == "25"] <- "0.25"
eff$alpha[eff$alpha == "50"] <- "0.5"
eff$alpha[eff$alpha == "75"] <- "0.75"
eff$alpha[eff$alpha == "100"] <- "1"
eff$alpha <- factor(eff$alpha, levels=c('0', '0.25', '0.5', '0.75', '1'))

marker_classes <- tmp[,c("marker_id", "class")]
marker_classes <- unique(marker_classes)

eff <- merge(eff,marker_classes, by="marker_id")


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


alpha.label <- paste("α = ", unique(eff$alpha))
names(alpha.label) <- unique(eff$alpha)
n_cells.label <- paste("#cells = ",unique(eff$nr_of_cells))
names(n_cells.label) <- unique(eff$nr_of_cells)

# rename and order methods
tmp[tmp == "t_test"] <- "t-test"
tmp[tmp == "wilcoxon_median"] <- "Wilcoxon test"
tmp[tmp == "kruskal_median"] <- "Kruskal-Wallis test"

tmp$method <- factor(tmp$method, levels=rev(c("diffcyt-DS-limma", "diffcyt-DS-LMM", "t-test", "Wilcoxon test","Kruskal-Wallis test", "CytoGLM","CytoGLMM", "logRegression", "ZAGA", "BEZI", "sceEMD")))

# only one marker class
ggplot(as.data.table(tmp[class == "State"]), aes(marker_id, method)) +
  geom_tile(aes(fill = significant), color = "white", size = 0.5) +
  ggtitle("") +
  xlab(label = " State Markers") +
  ylab("Method") +
  facet_grid(nr_of_cells ~ alpha, scales = "free", labeller = labeller(nr_of_cells = n_cells.label, alpha=alpha.label)) +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_manual(values = colorBlindBlack8[c(7, 3, 1)],
                    name = "Significant",
                    na.value = "transparent") +
  ggside::geom_xsidetile(
    data = eff[class == "State"],
    aes(y = overall_group, xfill = magnitude),
    color = "white"
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
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
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



  