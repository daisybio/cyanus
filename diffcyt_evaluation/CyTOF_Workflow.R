library(CATALYST)
RNGversion("3.5.3")

library(readxl)
url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
md <- "PBMC8_metadata.xlsx"
download.file(file.path(url, md), destfile = md, mode = "wb")
md <- read_excel(md)
head(data.frame(md))

library(HDCytoData)
fs <- Bodenmiller_BCR_XL_flowSet()

panel <- "PBMC8_panel_v3.xlsx"
download.file(file.path(url, panel), destfile = panel, mode = "wb")
panel <- read_excel(panel)
head(data.frame(panel))

all(panel$fcs_colname %in% colnames(fs))

# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

# construct SingleCellExperiment
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6
p

plotCounts(sce, group_by = "sample_id", color_by = "condition")

pbMDS(sce, color_by = "condition", label_by = "sample_id")

plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(10, "YlGnBu")))

plotNRS(sce, features = "type", color_by = "condition")

set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta20", 
                bars = TRUE, perc = TRUE)

plotClusterExprs(sce, k = "meta20", features = "type")

plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "pS6", k = "meta20", 
                 row_anno = FALSE, bars = TRUE, perc = TRUE)

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", cells = 500, features = "type")
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")
plotDR(sce, "UMAP", color_by = "CD4")
library(ggplot2)
p1 <- plotDR(sce, "TSNE", color_by = "meta20") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta20")
library(cowplot)
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))

# facet by sample
plotDR(sce, "UMAP", color_by = "meta20", facet_by = "sample_id")

# facet by condition
plotDR(sce, "UMAP", color_by = "meta20", facet_by = "condition")

plotCodes(sce, k = "meta20")

plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "pS6", k = "som100", m = "meta20", 
                 row_anno = FALSE, col_anno = FALSE, bars = TRUE, perc = TRUE)

####MANUAL CLUSTER MERGING####''
merging_table1 <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, merging_table1), 
              destfile = merging_table1, mode = "wb")
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))

# convert to factor with merged clusters in desired order
merging_table1$new_cluster <- factor(merging_table1$new_cluster, 
                                     levels = c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells",
                                                "CD8 T-cells", "DC", "NK cells", "monocytes", "surface-"))

# apply manual merging
sce <- mergeClusters(sce, k = "meta20", 
                     table = merging_table1, id = "merging1")

plotDR(sce, "UMAP", color_by = "merging1")



merging_table2 <- "PBMC8_cluster_merging2_v3.xlsx"
download.file(file.path(url, merging_table2), 
              destfile = merging_table2, mode = "wb")
merging_table2 <- read_excel(merging_table2)
data.frame(merging_table2)

# convert to factor with merged clusters in desired order
merging_table2$new_cluster <- factor(
  merging_table2$new_cluster, 
  levels = levels(merging_table1$new_cluster))

# apply manual merging
sce <- mergeClusters(sce, k = "meta12", 
                     table = merging_table2, id = "merging2")

# tabular comparison of algorithmic & manual merging
table(manual = cluster_codes(sce)[cluster_ids(sce), "merging2"],
      algorithm = cluster_codes(sce)[cluster_ids(sce), "meta8"] )

plot_grid(labels = c("A", "B"),
          plotDR(sce, "UMAP", color_by = "merging2"),
          plotDR(sce, "UMAP", color_by = "meta8"))

plotAbundances(sce, k = "merging1", by = "sample_id")
plotAbundances(sce, k = "merging1", by = "cluster_id", shape_by = "patient_id")

##### DIFFERENTIAL ANALYSIS #######
library(diffcyt)
ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id"))
(da_formula2 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = c("sample_id", "patient_id")))
contrast <- createContrast(c(0, 1))
da_res1 <- diffcyt(sce, 
                   formula = da_formula1, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "merging1", verbose = FALSE)
da_res2 <- diffcyt(sce, 
                   formula = da_formula2, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "merging1", verbose = FALSE)
topTable(da_res2, show_props = TRUE, format_vals = TRUE, digits = 2)

plotDiffHeatmap(sce, rowData(da_res2$res), all = TRUE, fdr = 0.05)
ds_formula1 <- createFormula(ei, cols_fixed = "condition")
ds_formula2 <- createFormula(ei, 
                             cols_fixed = "condition", cols_random = "patient_id")
ds_res1 <- diffcyt(sce, 
                   formula = ds_formula1, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)
ds_res2 <- diffcyt(sce, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)

plotDiffHeatmap(sce, rowData(ds_res2$res), top_n = 50, fdr = 0.05)

sce <- mergeClusters(sce, k = "meta20", id = "merging_all",
                     table = data.frame(old_cluster = seq_len(20), new_cluster = "all"))
saveRDS(sce, file = "data/cytof_workflow_SCE.rds")
