library(CATALYST)
library(data.table)
library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(magrittr)
sapply(list.files("functions", full.names = TRUE), source)

## 1. Read in SCE objects 

pathToData = "./data/diffcyt"
sceDiffcytMain <- readRDS(paste0(pathToData, "/sce_full.rds"))
sceDiffcyt0pc <- readRDS(paste0(pathToData, "/sce_0.rds"))
sceDiffcyt25pc <- readRDS(paste0(pathToData, "/sce_25.rds"))
sceDiffcyt50pc <- readRDS(paste0(pathToData, "/sce_50.rds"))
sceDiffcyt75pc <- readRDS(paste0(pathToData, "/sce_75.rds"))
sceDiffcyt100pc <- readRDS(paste0(pathToData, "/sce_100.rds"))

exprsDT <- as.data.table(t(assays(sceDiffcytMain)$counts))
coldataDT <- as.data.table(colData(sceDiffcytMain))[, c("patient_id", "condition", "cluster_id")]
exprsDT <- cbind(exprsDT[, c("pBtk", "pNFkB", "pSlp76", "pS6")], coldataDT)
exprsDT <- exprsDT[cluster_id == "B-cells IgM+" | cluster_id == "B-cells IgM-"]

pBtk_means <- exprsDT[, mean(pBtk), by = .(patient_id, condition, cluster_id)]
colnames(pBtk_means) <- c("patient_id", "condition", "cluster_id", "mean_full")
pNFkB_means <- exprsDT[, mean(pNFkB), by = .(patient_id, condition, cluster_id)]
colnames(pNFkB_means) <- c("patient_id", "condition", "cluster_id", "mean_full")
pSlp76_means <- exprsDT[, mean(pSlp76), by = .(patient_id, condition, cluster_id)]
colnames(pSlp76_means) <- c("patient_id", "condition", "cluster_id", "mean_full")
pS6_means <- exprsDT[, mean(pS6), by = .(patient_id, condition, cluster_id)]
colnames(pS6_means) <- c("patient_id", "condition", "cluster_id", "mean_full")
pBtk_means[, sd_full := exprsDT[, sd(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, sd_full := exprsDT[, sd(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, sd_full := exprsDT[, sd(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, sd_full := exprsDT[, sd(pS6), by = .(patient_id, condition, cluster_id)]$V1]

make_exprsDT <- function(sce){
  exprsDT <- as.data.table(t(assays(sce)$counts))
  coldataDT <- as.data.table(colData(sce))[, c("patient_id", "condition", "cluster_id")]
  exprsDT <- cbind(exprsDT[, c("pBtk", "pNFkB", "pSlp76", "pS6")], coldataDT)
  exprsDT <- exprsDT[cluster_id == "B-cells IgM+" | cluster_id == "B-cells IgM-"]
  return(exprsDT)
}

exprsDT <- make_exprsDT(sceDiffcyt0pc)
pBtk_means[, mean_0 := exprsDT[, mean(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, mean_0 := exprsDT[, mean(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, mean_0 := exprsDT[, mean(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, mean_0 := exprsDT[, mean(pS6), by = .(patient_id, condition, cluster_id)]$V1]
pBtk_means[, sd_0 := exprsDT[, sd(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, sd_0 := exprsDT[, sd(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, sd_0 := exprsDT[, sd(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, sd_0 := exprsDT[, sd(pS6), by = .(patient_id, condition, cluster_id)]$V1]

exprsDT <- make_exprsDT(sceDiffcyt25pc)
pBtk_means[, mean_25 := exprsDT[, mean(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, mean_25 := exprsDT[, mean(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, mean_25 := exprsDT[, mean(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, mean_25 := exprsDT[, mean(pS6), by = .(patient_id, condition, cluster_id)]$V1]
pBtk_means[, sd_25 := exprsDT[, sd(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, sd_25 := exprsDT[, sd(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, sd_25 := exprsDT[, sd(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, sd_25 := exprsDT[, sd(pS6), by = .(patient_id, condition, cluster_id)]$V1]

exprsDT <- make_exprsDT(sceDiffcyt50pc)
pBtk_means[, mean_50 := exprsDT[, mean(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, mean_50 := exprsDT[, mean(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, mean_50 := exprsDT[, mean(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, mean_50 := exprsDT[, mean(pS6), by = .(patient_id, condition, cluster_id)]$V1]
pBtk_means[, sd_50 := exprsDT[, sd(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, sd_50 := exprsDT[, sd(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, sd_50 := exprsDT[, sd(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, sd_50 := exprsDT[, sd(pS6), by = .(patient_id, condition, cluster_id)]$V1]

exprsDT <- make_exprsDT(sceDiffcyt75pc)
pBtk_means[, mean_75 := exprsDT[, mean(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, mean_75 := exprsDT[, mean(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, mean_75 := exprsDT[, mean(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, mean_75 := exprsDT[, mean(pS6), by = .(patient_id, condition, cluster_id)]$V1]
pBtk_means[, sd_75 := exprsDT[, sd(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, sd_75 := exprsDT[, sd(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, sd_75 := exprsDT[, sd(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, sd_75 := exprsDT[, sd(pS6), by = .(patient_id, condition, cluster_id)]$V1]

exprsDT <- make_exprsDT(sceDiffcyt100pc)
pBtk_means[, mean_100 := exprsDT[, mean(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, mean_100 := exprsDT[, mean(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, mean_100 := exprsDT[, mean(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, mean_100 := exprsDT[, mean(pS6), by = .(patient_id, condition, cluster_id)]$V1]
pBtk_means[, sd_100 := exprsDT[, sd(pBtk), by = .(patient_id, condition, cluster_id)]$V1]
pNFkB_means[, sd_100 := exprsDT[, sd(pNFkB), by = .(patient_id, condition, cluster_id)]$V1]
pSlp76_means[, sd_100 := exprsDT[, sd(pSlp76), by = .(patient_id, condition, cluster_id)]$V1]
pS6_means[, sd_100 := exprsDT[, sd(pS6), by = .(patient_id, condition, cluster_id)]$V1]

all_stats <- rbindlist(list("pBtk" = pBtk_means, 
                            "pNFkB" = pNFkB_means, 
                            "pSlp76" = pSlp76_means, 
                            "pS6" = pS6_means), idcol = "marker")

all_stats <- melt(all_stats, id.vars = c("marker", "patient_id", "condition", "cluster_id"), 
                  value.name = "value", variable.name="variable")

all_stats[, model := tstrsplit(variable, "_", keep=2)]
all_stats[, model := factor(model, levels = c("full", "0", "25", "50", "75", "100"))]
all_stats[, mean_sd := tstrsplit(variable, "_", keep=1)]
all_stats[, mean_sd := as.factor(mean_sd)]
all_stats[, variable := NULL]

all_stats <- tidyr::unite(all_stats, "joint_condition", condition:cluster_id, sep = ":")
all_stats[, joint_condition := as.factor(joint_condition)]
library(ggplot2)

ggplot(all_stats, aes(x = model, y = value, fill = joint_condition))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(mean_sd ~ marker, scales = "free")+
  scale_fill_manual(values = c("indianred1", "indianred3", "lightblue1", "lightblue3"))
  


