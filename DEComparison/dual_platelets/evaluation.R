# DUAL PLATELETS EVALUATION
library(CATALYST)

sce_dual <- readRDS("DataGeneration/platelets_dual/sce.rds")
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# BOXPLOT

features <- NULL
color_by = "platelets"
ncol <- 4
shape_by <- NULL
fun <- "median"
assay <- "exprs"
x <- sce_dual[CATALYST:::.get_features(sce_dual, features), ]
by <- "sample_id"
ms <- CATALYST:::.agg(x, by, fun, assay)
df <- melt(ms, varnames = c("antigen", by[length(by)]))
i <- match(df$sample_id, x$sample_id)
j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
df <- cbind(df, colData(x)[i, j])
ncs <- table(as.list(colData(x)[by]))
ncs <- rep(c(t(ncs)), each = nrow(x))
df <- df[ncs > 0, , drop = FALSE]
df$antigen <- as.character(df$antigen)
df[df == "PAC1"] <- "GPIIbIIIa"
df$antigen <- as.factor(df$antigen)

mean_data <- aggregate(value ~ platelets, data = df, mean)

boxplot <- ggplot(df, aes(x=platelets, y = value, fill = platelets))+
  geom_boxplot(width = 0.8, outlier.color = NA) +
  geom_point(alpha = 0.8, show.legend = FALSE, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) +
  geom_text(aes(label = sprintf("%.2f", mean(value)), y = mean(value)), 
            stat = "summary", fun.y = "mean", position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3, color = "black") +
  facet_wrap( ~ antigen, scales = "free", ncol=ncol) +
  scale_fill_manual(values =colorBlindBlack8[c(3,8)], name="Condition", labels=c("Activated", "Baseline")) + 
  xlab("") +
  ylab("Median Expression") +
  labs(color="Platelets") +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=16),
    axis.text = element_text(color = "black", size=14), 
    axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16))
  
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
  ylab("Normalized Density") +
  xlab("Expression") +
  theme_classic() + theme(
    panel.grid = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size=16),
    axis.text = element_text(color = "black", size=14), 
    axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16)) +  scale_color_manual(values =colorBlindBlack8[c(3,8)], name="Condition", labels=c("Activated", "Baseline"))


library("cowplot")
ggdraw() + 
  draw_plot(g, 0, 0.5, 1, 0.5) + 
  draw_plot(spike, 0, 0, 0.4, 0.5) +
  draw_plot(boxplot, 0.43, 0.25, 0.57, 0.25) + 
  draw_plot(exprs, 0.43, 0, 0.57, 0.25) + 
  draw_plot_label(c("A", "B", "C", "D"), c(0,0,0.4,0.4), c(1,0.5,0.5,0.25), size=20)



# HEATMAP
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

path_no_r <- "DEComparison/dual_platelets_no_random/sce_dual_res_timed.rds"
path_with_r <- "DEComparison/dual_platelets/sce_dual_res_timed.rds"   #with random effects
results_r <- readRDS(path_with_r)[["results"]]
results_no_r <- readRDS(path_no_r)[["results"]]
times_r <- readRDS(path_with_r)[["times"]]
times_no_r <- readRDS(path_no_r)[["times"]]
eff_r <- readRDS(path_with_r)[["eff"]]


tmp <- data.frame(method = results_r$method, marker_id = results_r$marker_id, p_adj = results_r$p_adj)
tmp$significant <- results_r$p_adj < 0.05
tmp$p_adj <- NULL
tmp <- as.data.table(tmp)
tmp$significant <- as.factor(tmp$significant)
tmp$class <- tmp$marker_id %in% c("CD62P", "CD63", "CD154", "CD107a")
tmp$class[tmp$class == TRUE] <- "State"
tmp$class[tmp$class == FALSE] <- "Type"
tmp$marker_id <- as.character(tmp$marker_id)
tmp[tmp == "PAC1"] <- "GPIIbIIIa"
tmp$marker_id <- as.factor(tmp$marker_id)

# add marker_id column and class column
marker_classes <- tmp[,c("marker_id", "class")]
marker_classes <- unique(marker_classes)
eff_r$marker_id <- sapply(strsplit(eff_r$group2,'::'), "[", 1)
eff_r[eff_r == "PAC1"] <- "GPIIb/IIIa"

eff_r <- merge(eff_r, marker_classes, by =c("marker_id"))

# rename and order methods
tmp[tmp == "t_test"] <- "t-test"
tmp[tmp == "wilcoxon_median"] <- "Wilcoxon test"
tmp[tmp == "kruskal_median"] <- "Kruskal-Wallis test"
tmp$method[tmp$method == "CyEMD"] <- "CyEMD"
tmp$method <- factor(tmp$method, levels=rev(c("diffcyt-DS-limma", "diffcyt-DS-LMM", "t-test", "Wilcoxon test","Kruskal-Wallis test", "CytoGLM","CytoGLMM", "logRegression", "ZAGA", "BEZI", "CyEMD")))

tmp$significant <- as.character(tmp$significant)
tmp$significant[tmp$significant == TRUE] <- "Yes"
tmp$significant[tmp$significant == FALSE] <- "No"

eff_r$magnitude <- as.character(eff_r$magnitude)
eff_r$magnitude[eff_r$magnitude == "small"] <- "Small"
eff_r$magnitude[eff_r$magnitude == "negligible"] <- "Negligible"
eff_r$magnitude[eff_r$magnitude == "large"] <- "Large"
eff_r$magnitude[eff_r$magnitude == "moderate"] <- "Moderate"
eff_r$magnitude <- factor(eff_r$magnitude, levels=c("Negligible", "Small", "Moderate", "Large"))

with_random_effect <- ggplot(tmp, aes(marker_id, method)) + 
  geom_tile(aes(fill=significant),color="white", size=1) + 
  ggtitle("") + xlab(label="Marker") + ylab("Method") +
  facet_wrap(~class, scales = "free_x") + 
  theme(text = element_text(size = 18),  axis.text.x = element_text(angle = 90, vjust=0.5))+
  scale_fill_manual(values = colorBlindBlack8[c(7,3)], na.value="transparent", name="Significant") + 
  ggside::geom_xsidetile(data=eff_r, aes(y=overall_group, xfill=magnitude), color="white", size=0.2) + 
  ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,4,6)], name='Cohen’s d\nEffect size\nMagnitude', na.value="transparent")

# plot heatmap without random
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

tmp$marker_id <- as.character(tmp$marker_id)
tmp[tmp == "PAC1"] <- "GPIIbIIIa"
tmp$marker_id <- as.factor(tmp$marker_id)

# rename and order methods
tmp[tmp == "t_test"] <- "t-test"
tmp[tmp == "wilcoxon_median"] <- "Wilcoxon test"
tmp[tmp == "kruskal_median"] <- "Kruskal-Wallis test"

tmp$method[tmp$method == "CyEMD"] <- "CyEMD"
tmp$method <- factor(tmp$method, levels=rev(c("diffcyt-DS-limma", "diffcyt-DS-LMM", "t-test", "Wilcoxon test","Kruskal-Wallis test", "CytoGLM","CytoGLMM", "logRegression", "ZAGA", "BEZI", "CyEMD")))

tmp$significant <- as.character(tmp$significant)
tmp$significant[tmp$significant == TRUE] <- "Yes"
tmp$significant[tmp$significant == FALSE] <- "No"

no_random_effect <- ggplot(tmp, aes(marker_id, method)) + 
  geom_tile(aes(fill=significant), color="white", size=1) + 
  ggtitle("") + xlab(label="Marker") + ylab("Method") +
  facet_wrap(~class, scales = "free_x") + 
  theme(text = element_text(size = 18),  axis.text.x = element_text(angle = 90, vjust=0.5))+
  scale_fill_manual(values = colorBlindBlack8[c(7,3)], na.value="transparent", name="Significant") + 
  ggside::geom_xsidetile(data=eff_r, aes(y=overall_group, xfill=magnitude), color="white", size=0.2) + 
  ggside::scale_xfill_manual(values=colorBlindBlack8[c(8,5,4,6)], name='Cohen’s d\nEffect size\nMagnitude', na.value="transparent")


library(ggpubr)
ggarrange(with_random_effect, no_random_effect, nrow=2,
          labels=c(' A Paired Analysis', 'B Unpaired Analysis'),
          font.label=list(size=18),
          legend = "right",
          common.legend = T)




# plot times
times_r$random_effect <- "Yes"
times_no_r$random_effect <- "No"

times <- rbind(times_r, times_no_r)

ggplot(times, aes(x = reorder(method,elapsed), y = elapsed, fill = random_effect))+ #+ facet_wrap(~random_effect) + 
  geom_bar(size = 3, stat = "identity", position="dodge") + labs(fill="Random Effect") + ylab("elapsed (in sec.)") + scale_fill_manual(values =colorBlindBlack8[c(2,4)]) +
  theme_bw() + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust=1)) + xlab(label="method")

#plot negligibe marker expressions
exprsDual <- assay(sce_dual, "exprs")
exprsDual <- exprsDual[c("CD61", "CD47", "CD9", "CD31"), ]
colnames(exprsDual) <- colData(sce_dual)$sample_id
exprsDualDT <- as.data.table(exprsDual)
exprsDualDT$marker <- rownames(exprsDual)
exprsDualDT <- melt(exprsDualDT, id.vars = "marker", variable.name = "patient_id", value.name = "expression")
exprsDualDT$platelets <- substr(exprsDualDT$patient_id, 8,8)
exprsDualDT$patient_id <- substr(exprsDualDT$patient_id, 1,7)


ggplot(exprsDualDT, aes(x = expression, color = platelets))+
  geom_density()+
  facet_grid(patient_id ~ marker)+
  scale_color_manual(values = colorBlindBlack8[c(3,8)], name = "Condition", labels=c("Activated", "Baseline"))+
  theme_bw() + 
  labs(x="Expression", y = "Density") +
  theme(
    #strip.background = element_blank(),
    strip.text = element_text(size=16),
    axis.text = element_text(color = "black", size=14), 
    axis.title = element_text(color = "black", size=16), legend.text = element_text(size=14), legend.title=element_text(size=16))


  


