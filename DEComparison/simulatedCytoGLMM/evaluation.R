# CYTOGLMM SIMULATION EVALUATION

source('DEComparison/benchmarking_plots.R')

stats_table <- preparePlotData('simulatedCytoGLMM')
plot_cells_vs_elapsed(stats_table)
plot_sens_vs_pre(stats_table)
plot_sens_vs_spec(stats_table)
plot_f1_vs_elapsed(stats_table)


# Table
stats_table <- preparePlotData("simulatedCytoGLMM")
stats_table <- stats_table[,c("method", "nr_of_cells","sensitivity", "specificity", "precision", "F1")]

library(dplyr)
stats_table <- stats_table %>%  group_by(.dots = c("method")) %>% summarise(sensitivity_mean = mean(sensitivity, na.rm = TRUE), 
                                                                                           sensitivity_sd = sd(sensitivity, na.rm = TRUE), 
                                                                                           specificity_mean = mean(specificity, na.rm = TRUE), 
                                                                                           specificity_sd = sd(specificity, na.rm = TRUE),
                                                                                           precision_mean = mean(precision, na.rm = TRUE),
                                                                                           precision_sd = sd(precision, na.rm = TRUE),
                                                                                           F1_mean = mean(F1, na.rm = TRUE),
                                                                                           F1_sd = sd(F1, na.rm = TRUE))


stats_table <- as.data.table(stats_table)
is.num <- sapply(stats_table, is.numeric)
method <- stats_table$method
stats_table <- sapply(stats_table[,..is.num], function(x){
  round(x,2)
})
stats_table <- as.data.table(stats_table)
stats_table[, method:=method]

stats_table$sensitivity <- paste(stats_table$sensitivity_mean, "+/-", stats_table$sensitivity_sd)
stats_table$specificity <- paste(stats_table$specificity_mean, "+/-", stats_table$specificity_sd)
stats_table$precision <- paste(stats_table$precision_mean, "+/-", stats_table$precision_sd)
stats_table$F1 <- paste(stats_table$F1_mean, "+/-", stats_table$F1_sd)

stats_table <- stats_table[,c("method","sensitivity", "specificity", "precision", "F1")]

library(xtable)
print(xtable(stats_table, type = "latex"))



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
tmp$nr_of_cells <- factor(tmp$nr_of_cells, levels=c("1000", "2000", "5000", "10000", "15000", "20000", "200000"))
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
eff$nr_of_cells <- factor(eff$nr_of_cells, levels=c("1000", "2000", "5000", "10000", "15000","20000", "200000"))
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


class.label <- unique(tmp$class)
names(class.label) <- unique(tmp$class)
n_cells.label <- paste("#cells = ",unique(tmp$nr_of_cells))
names(n_cells.label) <- unique(tmp$nr_of_cells)

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
  ggtitle("") + 
  xlab(label="Marker") + 
  ylab("Method") +
  #labs(tag = "Number of Cells") +
  #facet_wrap(~class, scales = "free_x") + 
  facet_grid(nr_of_cells ~ class, scales = "free", labeller = labeller(nr_of_cells = n_cells.label, class=class.label)) +
  theme(text = element_text(size = 16),  
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust=0.5, size=16),
        plot.tag=element_text(angle=-90),
        plot.tag.position=c(.85, 0.5)) +
  scale_fill_manual(values = colorBlindBlack8[c(7,3)], na.value="transparent", name="Significant") + 
  ggside::geom_xsidetile(data=eff, aes(y=overall_group, xfill=magnitude), color="white", size=0.4) + 
  ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,4,6)], name="Cohen’s d\nEffect size\nMagnitude", na.value="transparent", drop=FALSE)




# plot expressions cytoGLMM
cytoGLMMS <- lapply(list.files("/localscratch/quirinmanz/cytof_data/cytoGLMM_simulated/", full.names = T, pattern=".rds"), readRDS)
names(cytoGLMMS) <- tstrsplit(list.files("/localscratch/quirinmanz/cytof_data/cytoGLMM_simulated/", pattern=".rds"), ".rds", keep=1)[[1]]

exprsDT <- lapply(cytoGLMMS, function(x){as.data.table(t(assays(x)$exprs))})
exprsDT <- rbindlist(exprsDT, idcol = "filename")
exprsDT[, n_cells := tstrsplit(filename, "_", keep=3)]
exprsDT[, c("patient_id", "condition") := rbindlist(lapply(cytoGLMMS, function(x){as.data.table(colData(x))[, c("patient_id", "condition")]}))]
exprsDT_long <- melt(exprsDT, id.vars = c("patient_id", "condition","n_cells", "filename"), variable.name = "marker", value.name = "expression")
exprs_means <- exprsDT_long[, median(expression), by = c("marker", "condition", "n_cells", "patient_id")]
colnames(exprs_means) <- c("marker", "condition","n_cells", "patient_id", "median")
exprs_means[, n_cells := factor(n_cells, levels = c("1000", "2000", "5000", "10000", "15000", "200000"))]

exprs_medians_200000 <- exprsDT_long[n_cells == "200000", median(expression), by = c("marker", "condition")]
colnames(exprs_medians_200000) <- c("antigen", "condition", "expression")

library(CATALYST)
library(ggplot2)

g <- plotExprs(cytoGLMMS$simulated_cytoGLMM_200000_cells, features = c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09"), color_by = "condition")
g <- g +
  geom_vline(data=exprs_medians_200000[antigen %in% c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09")], aes(xintercept = expression, col = condition))+
  scale_color_manual(values = colorBlindBlack8[c(3,8)], name = "Condition:\nExpression\nand Medians", labels = c("Case", "Control"), breaks=c("case", "control")) +
  labs(x="Expression", y = "Normalized Density") +
  theme(
    panel.grid = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=16),
    axis.text = element_text(color = "black", size=14), 
    axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16))
g

