library(CATALYST)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)
library(data.table)
library(readxl)
library(stringr)

checkNullTable <- function(toCheck) {
  if (is.null(toCheck))
    return(data.frame("Nothing" = ""))
  else
    return(toCheck)
}


exp <-
  list.files(
    "/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/extdata/HumanPlatelets",
    pattern = "\\.fcs$",
    full.names = T
  )

names <-
  list.files(
    "/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/extdata/HumanPlatelets/",
    pattern = "*.fcs"
  )
samples <-
  sapply(names, function(x) {
    strsplit(x, split = ".fcs")[[1]]
  })
patients <-
  sapply(names, function(x) {
    str_sub(strsplit(x, split = "_")[[1]][1], end = -2)
  })
a_b <-
  sapply(names, function(x) {
    str_sub(strsplit(x, split = "_")[[1]][1], start = 8)
  })
dual_triple <-
  sapply(names, function(x) {
    strsplit(strsplit(x, split = "_")[[1]][2], split = ".fcs")[[1]][1]
  })
metadata <- data.table(
  "file_name" = names,
  "sample_id" = samples,
  "patient_id" = patients,
  "platelets" = a_b,
  "therapy" = dual_triple
)
m <- match("condition", names(metadata))



storage <-
  sprintf("%.2f MB" , file.size(
    list.files(
      "/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/extdata/HumanPlatelets",
      pattern = "\\.fcs$",
      full.names = T
    )
  ) / 1000000)
files <-
  list.files(
    "/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/extdata/HumanPlatelets",
    pattern = "\\.fcs$"
  )

fcs_data <- data.frame(name = files, size = storage)

panel <-
  read_excel(
    "/Users/lisiarend/Desktop/Uni/Master/SysBioMed/CyTOF/extdata/HumanPlatelets/panel_platelets.xlsx"
  )


sce <-
  prepData(
    exp,
    panel,
    metadata,
    transform = FALSE,
    md_cols = list(
      file = "file_name",
      id = "sample_id",
      factors = c("platelets", "therapy","patient_id")
    )
  )
saveRDS(sce, file = "sce_transformed.rds")
saveRDS(fcs_data, file = "fcs.rds")
saveRDS(metadata, file = "md.rds")
saveRDS(panel, file = "panel.rds")


# Info for PBMC data
exampleData <- "PBMC"
infoText <-
  "Data consisting of paired samples of healthy peripheral blood mononuclear cells (PBMCs), where one sample from each pair was stimulated with B cell receptor / Fc receptor cross-linker (BCR-XL). The dataset contains 16 samples (8 paired samples); a total of 172,791 cells; and a total of 24 protein markers."
paper <-
  "Bodenmiller B, Zunder ER, Finck R, et al. Multiplexed mass cytometry profiling of cellular states perturbed by small-molecule regulators. Nat Biotechnol. 2012;30(9):858-867. doi:10.1038/nbt.2317"
doiLink <- "https://doi.org/10.1038/nbt.2317"

value <- list(strong(sprintf("Info about %s Data :", exampleData)),
              div(infoText),
              a(href = doiLink, paper))
value
saveRDS(value, "help.rds")



# Info for Mouse Data
exampleData <- "Mouse"
infoText <-
  "Data from 10 replicates of mice bone marrow cells. A total of about 840 000 cells were measured using a CyTOF panel of N = 39 markers. "
paper <-
  "Samusik N, Good Z, Spitzer MH, Davis KL, Nolan GP. Automated mapping of phenotype space with single-cell data. Nat Methods. 2016;13(6):493-496. doi:10.1038/nmeth.3863"
doiLink <- "https://doi.org/10.1038/nmeth.3863"
value <- list(strong(sprintf("Info about %s Data :", exampleData)),
              div(infoText),
              a(href = doiLink, paper))
value
saveRDS(value, "help.rds")

# Info for Platelet Data
exampleData <- "Platelets"
infoText <-
  "Data from human platelets of patients with chronic coronary syndrome undergoing different therapy: dual antiplatelet therapy versus triple antiplatelet therapy, before and after platelet activation with 10µm TRAP. There are files of 7 patients with triple therapy and 12 patients with dual therapy (each in two condtions)."
paper <- ""
doiLink <- ""
value <- list(strong(sprintf("Info about %s Data :", exampleData)),
              div(infoText),
              a(href = doiLink, paper))
value
saveRDS(value, "help.rds")
