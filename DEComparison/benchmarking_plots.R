library(data.table)
library(ggplot2)

result_rds <-
  list.files(path = "DEComparison/simulatedCytoGLMM/", #downsampled_covid_spike/",
             pattern = "\\.rds$",
             full.names = T)
trues <- c("m01", "m02", "m03", "m04", "m05") #c("CD62P", "CD63", "CD107a", "CD154")


times <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["times"]], simplify = FALSE),
    idcol = "dataset")
times[, dataset := basename(dataset)]
# times[, alpha := tstrsplit(dataset, "_", keep = 4)]
# times[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
times[, nr_of_cells := tstrsplit(dataset, "_", keep = 3)] # 6)]
times[, nr_of_cells := as.numeric(nr_of_cells)]


results <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["results"]], simplify = FALSE),
    idcol = "dataset")
results[, dataset := basename(dataset)]
# results[, alpha := tstrsplit(dataset, "_", keep = 4)]
# results[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
results[, nr_of_cells := tstrsplit(dataset, "_", keep = 3)] # 6)]
results[, nr_of_cells := as.numeric(nr_of_cells)]

# NA's should be false classified

stats_table <- results[, .(
  TP = sum(
    p_adj <= 0.05 &
      marker_id %in% trues
      # & alpha != 100
    ,
    na.rm = T
  ),
  FP = sum((p_adj <= 0.05 | is.na(p_adj)) &
             !(
               marker_id %in% trues
             )),
  TN = sum(p_adj > 0.05 &
             !(
               marker_id %in% trues
             ), na.rm = T),
  FN = sum(
    (p_adj > 0.05 | is.na(p_adj)) &
      marker_id %in% trues
      # & alpha != 100
  )
), by = .(method, nr_of_cells)] # , alpha)]

stats_table[, sensitivity := TP / (TP + FN)]
# stats_table[is.nan(sensitivity)]$sensitivity <- 0
stats_table[, specificity := TN / (TN + FP)]
stats_table[, precision := TP / (TP + FP)]
stats_table[, F1 := 2 * (sensitivity * precision) / (sensitivity + precision)]
#stats_table[, nr_of_cells := factor(nr_of_cells, levels = c(15000, 10000, 5000, 2000, 1000))]

stats_table <-
  merge(stats_table, times[, .(method, nr_of_cells, elapsed)]) # , alpha)])

ggplot(stats_table,
       aes(
         x = nr_of_cells,
         y = elapsed,
         color = method
         # shape = alpha
       )) +
  geom_jitter(size = 3) +
  theme_bw() + theme(text = element_text(size = 18))

ggplot(stats_table,
       aes(
         x = 1 - specificity,
         y = sensitivity,
         color = method,
         shape = method
       )) +
  scale_shape_manual(values=1:stats_table[, uniqueN(method)]) +
  geom_point(size = 5, alpha = .6) +
  facet_wrap(~ nr_of_cells) +
  # facet_grid(nr_of_cells ~ alpha) +
  ylim(c(0.0, 1.0)) +
  xlim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggplot(stats_table,
       aes(
         x = sensitivity,
         y = precision,
         color = method,
         shape = method
       )) +
  scale_shape_manual(values=1:stats_table[, data.table::uniqueN(method)]) +
  geom_point(size = 5, alpha = .6) +
  # facet_wrap(~ nr_of_cells) +
  # facet_grid(nr_of_cells ~ alpha) +
  ylim(c(0.0, 1.0)) +
  xlim(c(0.0, 1.0)) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggplot(stats_table, aes(
  x = F1,
  y = elapsed,
  color = as.factor(nr_of_cells),# method,
  shape = as.factor(nr_of_cells)
)) + facet_wrap(~ method) +
  geom_point(size = 5, alpha = .6) +
  xlim(c(0, 1)) +
  theme_bw() + theme(text = element_text(size = 18))
