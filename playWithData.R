library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)

exp1 <- list.files("/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/extdata/MouseData/fcs/", pattern = "*.fcs", full.names = T)
sce <- prepData(exp1[[1]], transform = T)
print(rowData(sce))
chs <- c("DNA1", "DNA2")
plotScatter(sce, chs)

chs <- grep("^CD", rownames(sce), value = TRUE)
plotScatter(sce, chs)

sce <- prepData(exp1, transform = T)
plotCounts(sce, group_by = "sample_id", color_by = NULL)

#add information to rowData (cluster_id from FlowSOM), colData (marker_class type or state), 
#metadata (SOM_codes, cluster_codes, delta_area)
sce <- cluster(sce, feature = "state")    

#represents the amount of extra cluster stability gained when clustering into k groups as compared to k-1 groups.
delta_area(sce)

#expression values not comparable but marker expressions comparable
plotExprHeatmap(sce, by = "sample_id", scale = "first")
#each marker's value between 0 and 1. Comparability between markers lost but differences across samples visible
plotExprHeatmap(sce, by = "sample_id", scale = "last")
plotExprHeatmap(sce, by = "cluster_id", scale = "last")
plotExprHeatmap(sce, features = "DNA1",by = "both", scale = "last")

plotPbExprs(sce, facet_by = "cluster_id", color_by = "sample_id")

#do Dimensionsionality Reduction
set.seed(1601)
sce <- runDR(sce, dr = "UMAP", features = "state", cells = 500)
sce <- runDR(sce, dr = "TSNE", features = "state", cells = 500)
sce <- runDR(sce, dr = "PCA", features = "state", cells = 500)
sce <- runDR(sce, dr = "MDS", features = "state", cells = 500)
sce <- runDR(sce, dr = "DiffusionMap", features = "state", cells = 500)
reducedDimNames(sce)

plotDR(sce, dr = "UMAP", color_by = c("BC1", "BC2"))
plotDR(sce, dr = "TSNE", color_by = c("BC1", "BC2"))
plotDR(sce, dr = "PCA", color_by = c("BC1", "BC2"))
plotDR(sce, dr = "MDS", color_by = c("BC1", "BC2"))
plotDR(sce, dr = "DiffusionMap", color_by = c("BC1", "BC2"))


###other data
exp2 <- list.files("~/Downloads/SystemsBioMedicine/systems_biomedicine_spill_applied_fcs_files/", pattern = "*.fcs", full.names = T)
names <- list.files("~/Downloads/SystemsBioMedicine/systems_biomedicine_spill_applied_fcs_files/", pattern = "*.fcs")
exp2 <- exp2[c(1:6)]
names <- names[c(1:6)]
info <- sapply(names, function(x){strsplit(x, split = "RPS....")[[1]][2]})
a_b <- sapply(info, function(x){factor(strsplit(x, split = "_")[[1]][1])})
dual_triple <- sapply(info, function(x){
  factor(
    strsplit(strsplit(x, split = "_")[[1]][2], split = "\\.")[[1]][1]
    )
  })
library(data.table)
metadata <- data.table("file_name" = names, 
                       "sample_id" = names, "ab" = a_b, "dualTriple" = dual_triple)

sce2 <- prepData(exp2, md = metadata, md_cols = list(file = "file_name", 
                                                     id = "sample_id", 
                                                     factors = c("ab", "dualTriple")))

#do Dimensionsionality Reduction
set.seed(1601)
sce2 <- runDR(sce2, dr = "UMAP", features = "state", cells = 1000)
sce2 <- runDR(sce2, dr = "TSNE", features = "state", cells = 1000)
chs <- grep("^CD", rownames(sce2), value = TRUE)[c(1:5)]

plotDR(sce2, dr = "UMAP", color_by = "ab", facet_by = "dualTriple")
plotDR(sce2, dr = "UMAP", color_by = chs, facet_by = "dualTriple")
plotDR(sce2, dr = "TSNE", color_by = "ab", facet_by = "dualTriple")
plotDR(sce2, dr = "TSNE", color_by = chs, facet_by = "dualTriple")

library(CytobankAPI)
#flowsom, viSNE
library(vegan)
#isomap
dis <- vegdist(assays(sce2)$exprs)
tr <- spantree(dis)
pl <- ordiplot(cmdscale(dis), main="cmdscale")
lines(tr, pl, col="red")
ord <- isomap(dis, epsilon = 0.9, ndim=1)

plotExprs(sce2, color_by = "dualTriple")

