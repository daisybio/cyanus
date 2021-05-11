library(data.table)
library(ggplot2)

result_rds <- list.files(path = "DEComparison/first_benchmark/", pattern = "\\.rds$", full.names = T)
results <- data.table::rbindlist(sapply(result_rds, function(rds) readRDS(rds)[["times"]], simplify = FALSE), idcol = "dataset")
results[, dataset := basename(dataset)]
results[, alpha := tstrsplit(dataset, "_", keep=4)]
results[, nr_of_cells := tstrsplit(dataset, "_", keep=6)]
results[, nr_of_cells := as.numeric(nr_of_cells)]

ggplot(results, aes(x = nr_of_cells, y = elapsed, color = method, shape = alpha)) +
  #geom_point(size = 2) +
  geom_jitter(size = 3) +
  theme_bw() + theme(text=element_text(size=18))
