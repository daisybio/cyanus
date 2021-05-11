sceObjects <- lapply(list.files("/nfs/home/students/ga89koc/hiwi/cytof/DEComparison/first_benchmark/", full.names = T), readRDS)
library(data.table)
results <- lapply(sceObjects, function(x){return(x[["results"]])})
names(results) <- list.files("/nfs/home/students/ga89koc/hiwi/cytof/DEComparison/first_benchmark/")
results <- rbindlist(results, idcol = "file")
results[, alpha := tstrsplit(file, "_", keep=4)]
results[, nr_of_cells := tstrsplit(file, "_", keep=6)]

TPs <- results[, sum(p_adj <= 0.05 & marker_id %in% c("CD62P", "CD63", "CD107a", "CD154") & alpha != 100, na.rm = T), by=.(method, nr_of_cells, alpha)]
colnames(TPs) <- c("method", "nr_of_cells", "alpha", "TP")
FPs <- results[, sum(p_adj <= 0.05 & !(marker_id %in% c("CD62P", "CD62", "CD107a", "CD154")), na.rm = T), by=.(method, nr_of_cells, alpha)]
colnames(FPs) <- c("method", "nr_of_cells", "alpha", "FP")
TNs <- results[, sum(p_adj > 0.05 & !(marker_id %in% c("CD62P", "CD62", "CD107a", "CD154")), na.rm = T), by=.(method, nr_of_cells, alpha)]
colnames(TNs) <- c("method", "nr_of_cells", "alpha", "TN")
FNs <- results[, sum(p_adj > 0.05 & marker_id %in% c("CD62P", "CD62", "CD107a", "CD154") & alpha != 100, na.rm = T), by=.(method, nr_of_cells, alpha)]
colnames(FNs) <- c("method", "nr_of_cells", "alpha", "FN")

stats_table <- data.table(TPs, FP = FPs$FP, TN = TNs$TN, FN = FNs$FN)
stats_table[, sensitivity := TP / (TP+FN)]
stats_table[is.nan(sensitivity)]$sensitivity <- 0
stats_table[, specificity := TN / (TN+FP)]
stats_table[, nr_of_cells := factor(nr_of_cells, levels = c(15000, 10000, 5000, 2000, 1000))]
stats_table[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]

library(ggplot2)
ggplot(stats_table, aes(x = 1-specificity, y = sensitivity, color = method, shape = method))+
  geom_point(size=3)+
  facet_grid(nr_of_cells ~ alpha)+
  ylim(c(0.0, 1.0))+
  theme_bw()+
  theme(text = element_text(size=20))
