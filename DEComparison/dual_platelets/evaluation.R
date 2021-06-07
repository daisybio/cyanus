# DUAL PLATELETS EVALUATION

path <- "DEComparison/dual_platelets"
path <- "/nfs/home/students/jbernett/cytof/cytof/DEComparison/dual_platelets"

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


times_re$random_effect <- "with random_effect"
times$random_effect <- "without random_effect"

tmp <- rbind(times, times_re)

# plot times
library(ggplot2)
times$usage <- c("median", "median", "whole expression", "whole expression", "whole expression", "whole expression", "whole expression", "whole expression","median", "median")
ggplot(tmp, aes(x = reorder(method,elapsed), y = elapsed, fill = random_effect))+ #+ facet_wrap(~random_effect) + 
  geom_bar(size = 3, stat = "identity", position="dodge") +
  theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust=1)) + xlab(label="method") + ggtitle("Dual Platelets")

#ggplot(times, aes(x = method, y = elapsed, fill = usage)) +
#  geom_bar(size = 3, stat = "identity") +
#  theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 30))

# plot results
tmp <- data.frame(method = results$method, marker_id = results$marker_id, p_adj = results$p_adj)
tmp$significant <- results$p_adj < 0.05
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
tmp$class <- tmp$marker_id %in% c("CD62P", "CD63", "CD154", "CD107a")
tmp$class[tmp$class == TRUE] <- "state"
tmp$class[tmp$class == FALSE] <- "type"

ggplot(tmp, aes(marker_id, method, fill=significant)) + geom_tile(color="white", size=1) + ggtitle("Dual Platelets without random effect") + xlab(label="marker") + facet_wrap(~class, scales = "free_x") + theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))


