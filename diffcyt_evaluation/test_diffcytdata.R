library(CATALYST)
library(data.table)
library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(magrittr)
sapply(list.files("functions", full.names = TRUE), source)

## 1. Read in SCE objects 

pathToData = "../data/diffcyt"
sceDiffcytMain <- readRDS(paste0(pathToData, "/main/sce.rds"))
sceDiffcyt75pc <- readRDS(paste0(pathToData, "/less_75pc/sce.rds"))  
sceDiffcyt50pc <- readRDS(paste0(pathToData, "/less_50pc/sce.rds"))

## 2. Start Analysis

### 2.1. Main
#prepare parameters
designMain <- createDesignMatrix(
  ei(sceDiffcytMain),
  cols_design = c("condition", "patient_id")
)
contrastMainLimma <- createContrast(
  c(0, 1, rep(0, 7))
)
contrastMainLMM <- createContrast(
  c(0, 1)
)
formulaMain <- createFormula(
  ei(sceDiffcytMain),
  cols_fixed = "condition",
  cols_random = "patient_id"
)
markersToTest <- getMarkersToTest(sceDiffcytMain, "limma", "all")

