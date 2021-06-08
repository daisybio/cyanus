library(data.table)
library(ggplot2)
library(data.table)

preparePlotData <- function(data_type){
  if (data_type == "simulatedCytoGLMM"){
    path <-  "DEComparison/simulatedCytoGLMM"
    trues <- c("m01", "m02", "m03", "m04", "m05")
    keep <- 3
  } else if (data_type == "downsampled_covid_spike"){
    path <-  "DEComparison/downsampled_covid_spike"
    trues <- c("CD62P", "CD63", "CD107a", "CD154")
    keep <- 6
  } else if (data_type == "dual_platelets"){
    path <-  "DEComparison/dual_platelets"
    trues <- c("CD62P", "CD63", "CD107a", "CD154")
  }
  
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
  
  if (data_type=="downsampled_covid_spike"){
    times[, alpha := tstrsplit(dataset, "_", keep = 4)]
    times[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
    results[, alpha := tstrsplit(dataset, "_", keep = 4)]
    results[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
  }
  
  if (data_type!="dual_platelets"){
    # because dual has no subsampling
    times[, nr_of_cells := tstrsplit(dataset, "_", keep = keep)]
    times[, nr_of_cells := as.numeric(nr_of_cells)]
    results[, nr_of_cells := tstrsplit(dataset, "_", keep = keep)]
    results[, nr_of_cells := as.numeric(nr_of_cells)]
    
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
    
  } else {
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
    ), by = .(method)]
  }
  
  
  stats_table[, sensitivity := TP / (TP + FN)]
  # stats_table[is.nan(sensitivity)]$sensitivity <- 0
  stats_table[, specificity := TN / (TN + FP)]
  stats_table[, precision := TP / (TP + FP)]
  stats_table[, F1 := 2 * (sensitivity * precision) / (sensitivity + precision)]
  
  if (data_type=="downsampled_covid_spike"){
    stats_table[, nr_of_cells := factor(nr_of_cells, levels = c(15000, 10000, 5000, 2000, 1000))]
  }
  
  if (data_type!="dual_platelets"){
    times$nr_of_cells <- as.numeric(times$nr_of_cells)
    stats_table$nr_of_cells <- as.numeric(stats_table$nr_of_cells)
    stats_table <-
      merge(stats_table, times[, .(method, nr_of_cells, elapsed)]) # , alpha)])
  } else {
    stats_table <-
      merge(stats_table, times[, .(method,elapsed)]) # , alpha)])
  }
  
  return (stats_table)
}

plotHeatmaps <- function(data_type, nr_cells_spike=15000){
  if (data_type == "simulatedCytoGLMM"){
    path <-  "DEComparison/simulatedCytoGLMM"
    trues <- c("m01", "m02", "m03", "m04", "m05")
    keep <- 3
  } else if (data_type == "downsampled_covid_spike"){
    path <-  "DEComparison/downsampled_covid_spike"
    trues <- c("CD62P", "CD63", "CD107a", "CD154")
    keep <- 6
  } else if (data_type == "dual_platelets"){
    path <-  "DEComparison/dual_platelets"
    trues <- c("CD62P", "CD63", "CD107a", "CD154")
  } else if(data_type == "pbmc"){
    path <- "DEComparison/pbmc_benchmarking"
  }
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
  if (data_type %in% c("simulatedCytoGLMM", "downsampled_covid_spike")){
    results[, nr_of_cells := tstrsplit(dataset, "_", keep=keep)]
    tmp$nr_of_cells <- results$nr_of_cells
    if(data_type == "downsampled_covid_spike"){
      results[, alpha := tstrsplit(dataset, "_", keep=4)]
      tmp$alpha <- results$alpha
    }
  }else if(data_type == "pbmc"){
    tmp$cluster_id <- results$cluster_id
  }
  times <-
    data.table::rbindlist(sapply(result_rds, function(rds)
      readRDS(rds)[["times"]], simplify = FALSE),
      idcol = "dataset")
  times[, dataset := basename(dataset)]
  
  # plot results
  tmp$significant <- results$p_adj < 0.05
  tmp$p_adj <- NULL
  tmp <- as.data.table(tmp)
  tmp$significant <- as.factor(tmp$significant)
  tmp$class <- tmp$marker_id %in% trues
  if(data_type %in% c("downsampled_covid_spike", "dual_platelets")){
    tmp$class[tmp$class == TRUE] <- "state"
    tmp$class[tmp$class == FALSE] <- "type"
  }else if(data_type == "simulatedCytoGLMM"){
    tmp$class[tmp$class == TRUE] <- "differentially expressed"
    tmp$class[tmp$class == FALSE] <- "not differentially expressed"
  }
  
  
  if(data_type == "downsampled_covid_spike"){
    ggplot(tmp, aes(marker_id, method, fill=significant)) + 
      geom_tile(color="white", size=1) + 
      ggtitle(sprintf("Downsampled Covid Spike, %d cells per sample", nr_cells_spike)) + 
      xlab(label="marker") + 
      facet_grid(alpha~class, scales = "free_x") + 
      theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))
  } else if(data_type == "simulatedCytoGLMM"){
  ggplot(tmp, aes(marker_id, method, fill=significant)) + 
    geom_tile(color="white", size=1) + 
    ggtitle("CytoGLMM Simulation") + 
    xlab(label="marker") + 
    facet_grid(nr_of_cells~class, scales = "free_x") + 
    theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))
  } else if(data_type == "dual_platelets"){
    ggplot(tmp, aes(marker_id, method, fill=significant)) + 
      geom_tile(color="white", size=1) + 
      ggtitle("Dual Platelets") + xlab(label="marker") + 
      facet_wrap(~class, scales = "free_x") + 
      theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))
  } else if(data_type == "pbmc"){
    tmp$marker_id[tmp$marker_id == "HLADR"] <- "HLA_DR"
    ggplot(tmp, aes(marker_id, method, fill=significant)) + 
      geom_tile(color="white", size=1) + 
      ggtitle("PBMC Ref vs. BCR-XL") + xlab(label="marker") + 
      facet_wrap(~cluster_id, scales = "free_x") + 
      theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))
  }
}

plot_cells_vs_elapsed <- function(stats_table){
  if ("alpha" %in% stats_table){
    plot <- ggplot(stats_table,
                   aes(
                     x = nr_of_cells,
                     y = elapsed,
                     color = method,
                     shape = alpha
                   )) +
      geom_jitter(size = 3) +
      theme_bw() + theme(text = element_text(size = 18))
  } else {
    plot <- ggplot(stats_table,
         aes(
           x = nr_of_cells,
           y = elapsed,
           color = method
         )) +
    geom_jitter(size = 3) +
    theme_bw() + theme(text = element_text(size = 18))
  }
  return(plot)
}

plot_sens_vs_spec <- function(stats_table){
  plot <- ggplot(stats_table,
                 aes(
                   x = 1 - specificity,
                   y = sensitivity,
                   color = method,
                   shape = method
                 )) +
    scale_shape_manual(values=1:stats_table[, uniqueN(method)]) +
    geom_point(size = 5, alpha = .6) +
    ylim(c(0.0, 1.0)) +
    xlim(c(0.0, 1.0)) +
    theme_bw() +
    theme(text = element_text(size = 20))
  if ("nr_of_cells" %in% colnames(stats_table)){
    if ("alpha" %in% colnames(stats_table)){
      plot <- plot + facet_grid(nr_of_cells~alpha)
    }
    plot <- plot + facet_wrap(nr_of_cells)
  }
  return(plot)
}

plot_sens_vs_pre <- function(stats_table){
    plot <- ggplot(stats_table,
                   aes(
                     x = sensitivity,
                     y = precision,
                     color = method,
                     shape = method
                   )) +
      scale_shape_manual(values=1:stats_table[, data.table::uniqueN(method)]) +
      geom_point(size = 5, alpha = .6) +
      ylim(c(0.0, 1.0)) +
      xlim(c(0.0, 1.0)) +
      theme_bw() +
      theme(text = element_text(size = 20))
    
    if ("nr_of_cells" %in% colnames(stats_table)){
      if ("alpha" %in% colnames(stats_table)){
        plot <- plot +facet_grid(nr_of_cells ~ alpha)
      } else {
        plot <- plot +facet_wrap(nr_of_cells)
      }
    }
  return(plot)
}

plot_f1_vs_elapsed <- function(stats_table){
  if ("nr_of_cells" %in% stats_table){
    plot <- ggplot(stats_table, aes(
      x = F1,
      y = elapsed,
      color = as.factor(nr_of_cells),
      shape = as.factor(nr_of_cells)
    )) + facet_wrap(~ method) +
      geom_point(size = 5, alpha = .6) +
      xlim(c(0, 1)) +
      theme_bw() + theme(text = element_text(size = 18))
  } else {
    plot <- ggplot(stats_table, aes(
      x = F1,
      y = elapsed,
      color = as.factor(method),
    )) +
      geom_point(size = 5, alpha = .6) +
      xlim(c(0, 1)) +
      theme_bw() + theme(text = element_text(size = 18))
  }
 
 return(plot)
}

#result_rds <-
#  list.files(path = "DEComparison/simulatedCytoGLMM/", #downsampled_covid_spike/",
#             pattern = "\\.rds$",
#             full.names = T)
#trues <- c("m01", "m02", "m03", "m04", "m05") #c("CD62P", "CD63", "CD107a", "CD154")


#times <-
#  data.table::rbindlist(sapply(result_rds, function(rds)
#    readRDS(rds)[["times"]], simplify = FALSE),
#    idcol = "dataset")
#times[, dataset := basename(dataset)]
# times[, alpha := tstrsplit(dataset, "_", keep = 4)]
# times[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
#times[, nr_of_cells := tstrsplit(dataset, "_", keep = 3)] # 6)]
#times[, nr_of_cells := as.numeric(nr_of_cells)]


#results <-
#  data.table::rbindlist(sapply(result_rds, function(rds)
#    readRDS(rds)[["results"]], simplify = FALSE),
#    idcol = "dataset")
#results[, dataset := basename(dataset)]
# results[, alpha := tstrsplit(dataset, "_", keep = 4)]
# results[, alpha := factor(alpha, levels = c("full", "25", "50", "75", "100"))]
#results[, nr_of_cells := tstrsplit(dataset, "_", keep = 3)] # 6)]
#results[, nr_of_cells := as.numeric(nr_of_cells)]

# NA's should be false classified

# stats_table <- results[, .(
#   TP = sum(
#     p_adj <= 0.05 &
#       marker_id %in% trues
#       # & alpha != 100
#     ,
#     na.rm = T
#   ),
#   FP = sum((p_adj <= 0.05 | is.na(p_adj)) &
#              !(
#                marker_id %in% trues
#              )),
#   TN = sum(p_adj > 0.05 &
#              !(
#                marker_id %in% trues
#              ), na.rm = T),
#   FN = sum(
#     (p_adj > 0.05 | is.na(p_adj)) &
#       marker_id %in% trues
#       # & alpha != 100
#   )
# ), by = .(method, nr_of_cells)] # , alpha)]

# stats_table[, sensitivity := TP / (TP + FN)]
# # stats_table[is.nan(sensitivity)]$sensitivity <- 0
# stats_table[, specificity := TN / (TN + FP)]
# stats_table[, precision := TP / (TP + FP)]
# stats_table[, F1 := 2 * (sensitivity * precision) / (sensitivity + precision)]
# #stats_table[, nr_of_cells := factor(nr_of_cells, levels = c(15000, 10000, 5000, 2000, 1000))]

#stats_table <-
#  merge(stats_table, times[, .(method, nr_of_cells, elapsed)]) # , alpha)])

# ggplot(stats_table,
#        aes(
#          x = nr_of_cells,
#          y = elapsed,
#          color = method
#          # shape = alpha
#        )) +
#   geom_jitter(size = 3) +
#   theme_bw() + theme(text = element_text(size = 18))
# 
# ggplot(stats_table,
#        aes(
#          x = 1 - specificity,
#          y = sensitivity,
#          color = method,
#          shape = method
#        )) +
#   scale_shape_manual(values=1:stats_table[, uniqueN(method)]) +
#   geom_point(size = 5, alpha = .6) +
#   facet_wrap(~ nr_of_cells) +
#   # facet_grid(nr_of_cells ~ alpha) +
#   ylim(c(0.0, 1.0)) +
#   xlim(c(0.0, 1.0)) +
#   theme_bw() +
#   theme(text = element_text(size = 20))
# 
# ggplot(stats_table,
#        aes(
#          x = sensitivity,
#          y = precision,
#          color = method,
#          shape = method
#        )) +
#   scale_shape_manual(values=1:stats_table[, data.table::uniqueN(method)]) +
#   geom_point(size = 5, alpha = .6) +
#   # facet_wrap(~ nr_of_cells) +
#   # facet_grid(nr_of_cells ~ alpha) +
#   ylim(c(0.0, 1.0)) +
#   xlim(c(0.0, 1.0)) +
#   theme_bw() +
#   theme(text = element_text(size = 20))
# 
# ggplot(stats_table, aes(
#   x = F1,
#   y = elapsed,
#   color = as.factor(nr_of_cells),# method,
#   shape = as.factor(nr_of_cells)
# )) + facet_wrap(~ method) +
#   geom_point(size = 5, alpha = .6) +
#   xlim(c(0, 1)) +
#   theme_bw() + theme(text = element_text(size = 18))
