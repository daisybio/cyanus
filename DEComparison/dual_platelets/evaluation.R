# DUAL PLATELETS EVALUATION

path <- "DEComparison/dual_platelets"

result_rds <-
  list.files(path = path,
             pattern = "\\.rds$",
             full.names = T)
results <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["results"]], simplify = FALSE),
    idcol = "dataset")
results[, dataset := basename(dataset)]
times <-
  data.table::rbindlist(sapply(result_rds, function(rds)
    readRDS(rds)[["times"]], simplify = FALSE),
    idcol = "dataset")
times[, dataset := basename(dataset)]


# plot times
library(ggplot2)
times$usage <- c("median", "median", "whole expression", "whole expression", "whole expression", "whole expression", "whole expression", "whole expression","median", "median")

ggplot(times, aes(x = method, y = elapsed, fill = method)) +
  geom_bar(size = 3, stat = "identity") +
  theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 30))

ggplot(times, aes(x = method, y = elapsed, fill = usage)) +
  geom_bar(size = 3, stat = "identity") +
  theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 30))

# plot results
tmp <- data.frame(method = results$method, marker_id = results$marker_id, p_adj = results$p_adj)
tmp$significant <- results$p_adj < 0.05
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
tmp$class <- tmp$marker_id %in% c("CD62P", "CD63", "CD154", "CD107a")
tmp$class[tmp$class == TRUE] <- "state"
tmp$class[tmp$class == FALSE] <- "type"

ggplot(tmp, aes(marker_id, method, fill=significant)) + geom_tile() + ggtitle("Dual Platelets") + xlab(label="marker") + facet_wrap(~class, scales = "free_x")


