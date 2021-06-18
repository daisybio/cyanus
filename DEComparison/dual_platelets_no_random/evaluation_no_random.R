# DUAL PLATELETS EVALUATION
library(CATALYST)

sce_dual <- readRDS("/nfs/home/students/l.arend/data/platelets_dual/sce_dual.rds")

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


boxplot <- plotPbExprs(sce_dual, features = c("CD62P", "CD63", "CD154", "CD107a"), color_by = "platelets", ncol=4) + 
  xlab("") +
  ylab("Median Expression") +
  labs(color="Platelets") +
  theme(
  panel.grid = element_blank(), 
  strip.background = element_blank(),
  strip.text = element_text(face = "bold", size=16),
  axis.text = element_text(color = "black", size=14), 
  axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16)) + scale_color_manual(values =colorBlindBlack8[c(3,8)])

# exprs function from catalyst
# subset features to use
x <- sce_dual
features <- c("CD62P", "CD63", "CD154", "CD107a")
assay <- "exprs"
features <- CATALYST:::.get_features(x, features)
y <- assay(x, assay)[features, ]
color_by <- "platelets"

# construct data.frame include cell metadata
df <- data.frame(t(y), colData(x), check.names = FALSE)
value <- ifelse(assay == "exprs", "expression", assay)

gg_df <- melt(df, 
              value.name = value,
              variable.name = "antigen", 
              id.vars = names(colData(x)))

exprs <- ggplot(gg_df, fill = NULL, 
       aes_string(
         x = value, y = "..ndensity..",
         col = color_by, group = "sample_id")) + 
  facet_wrap(~ antigen, scales = "free_x", ncol=4) +
  geom_density() + 
  ylab("Normalized density") +
  xlab("Expression") +
  theme_classic() + theme(
    panel.grid = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=16),
    axis.text = element_text(color = "black", size=14), 
    axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16)) +  scale_color_manual(values =colorBlindBlack8[c(3,8)])

library(ggpubr)
ggarrange(boxplot, exprs, ncol=1,
          labels=c('A', 'B'),
          font.label=list(size=18),
          legend = "right",
          common.legend = T)


path_no_r <- "/nfs/home/students/l.arend/cytof/DEComparison/dual_platelets_no_random/sce_dual_res_timed.rds"
path_with_r <- "/nfs/home/students/l.arend/cytof/DEComparison/dual_platelets/sce_dual_res_timed.rds"   #with random effects

results_r <- readRDS(path_with_r)[["results"]]
results_no_r <- readRDS(path_no_r)[["results"]]

times_r <- readRDS(path_with_r)[["times"]]
times_no_r <- readRDS(path_no_r)[["times"]]

eff_r <- readRDS(path_with_r)[["eff"]]

# plot results
tmp <- data.frame(method = results_r$method, marker_id = results_r$marker_id, p_adj = results_r$p_adj)
tmp$significant <- results_r$p_adj < 0.05
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
tmp$class <- tmp$marker_id %in% c("CD62P", "CD63", "CD154", "CD107a")
tmp$class[tmp$class == TRUE] <- "State"
tmp$class[tmp$class == FALSE] <- "Type"
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# add marker_id column and class column
marker_classes <- tmp[,c("marker_id", "class")]
marker_classes <- unique(marker_classes)
eff_r$marker_id <- sapply(strsplit(eff_r$group2,'::'), "[", 1)
eff_r <- merge(eff_r, marker_classes, by ="marker_id")


with_random_effect <- ggplot(tmp, aes(marker_id, method)) + 
  geom_tile(aes(fill=significant),color="white", size=1) + 
  ggtitle("") + xlab(label="Marker") + ylab("Method") +
  facet_wrap(~class, scales = "free_x") + 
  theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values = colorBlindBlack8[c(7,3)], na.value="transparent", name="Significant") + 
  ggside::geom_xsidetile(data=eff_r, aes(y=overall_group, xfill=magnitude), color="white", size=0.2) + 
  ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,2,6)], name='Effect size\nMagnitude', na.value="transparent")

# plot heatmap
tmp <- data.frame(method = results_no_r$method, marker_id = results_no_r$marker_id, p_adj = results_no_r$p_adj)
tmp$significant <- results_no_r$p_adj < 0.05
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
tmp$class <- tmp$marker_id %in% c("CD62P", "CD63", "CD154", "CD107a")
tmp$class[tmp$class == TRUE] <- "State"
tmp$class[tmp$class == FALSE] <- "Type"
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


no_random_effect <- ggplot(tmp, aes(marker_id, method)) + 
  geom_tile(aes(fill=significant), color="white", size=1) + 
  ggtitle("") + xlab(label="Marker") + ylab("Method") +
  facet_wrap(~class, scales = "free_x") + 
  theme(text = element_text(size = 16),  axis.text.x = element_text(angle = 45, hjust=1))+
  scale_fill_manual(values = colorBlindBlack8[c(7,3)], na.value="transparent", name="Significant") + 
  ggside::geom_xsidetile(data=eff_r, aes(y=overall_group, xfill=magnitude), color="white", size=0.2) + 
  ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,2,6)], name='Effect size\nMagnitude', na.value="transparent")


library(ggpubr)
ggarrange(with_random_effect, no_random_effect, nrow=2,
          labels=c('A Random Effect', 'B No Random Effect'),
          font.label=list(size=16),
          legend = "right",
          common.legend = T)


# plot times
times_r$random_effect <- "Yes"
times_no_r$random_effect <- "No"

times <- rbind(times_r, times_no_r)

ggplot(times, aes(x = reorder(method,elapsed), y = elapsed, fill = random_effect))+ #+ facet_wrap(~random_effect) + 
  geom_bar(size = 3, stat = "identity", position="dodge") + labs(fill="Random Effect") + ylab("elapsed (in sec.)") + scale_fill_manual(values =colorBlindBlack8[c(2,4)]) +
  theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust=1)) + xlab(label="method")





