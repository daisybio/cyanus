checkNullTable <- function(toCheck) {
  if (is.null(toCheck))
    return(data.frame("Nothing" = ""))
  else
    return(toCheck)
}

observeEvent(input$fcsFiles, {
  fileTable <- input$fcsFiles
  fileTable <- fileTable[, c("name", "size")]
  fileTable$size <- sprintf("%.2f MB", fileTable$size / 1000000)
  reactiveVals$fcs <- fileTable
})

observeEvent(input$metaFile, {
  if(endsWith(input$metaFile$datapath, ".csv")){
    reactiveVals$md <- read.table(input$metaFile$datapath, header = T, sep = ",")
  }else if(endsWith(input$metaFile$datapath, ".xls") | 
           endsWith(input$metaFile$datapath, ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    reactiveVals$md <- read.xlsx2(input$metaFile$datapath, 1)
  }
})

observeEvent(input$panelFile, {
  if(endsWith(input$panelFile$datapath, ".csv")){
    reactiveVals$panel <-
      read.table(input$panelFile$datapath, header = T, sep = ",")
  }else if(endsWith(input$panelFile$datapath, ".xls") | 
            endsWith(input$panelFile$datapath, ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    reactiveVals$panel <- read.xlsx2(input$panelFile$datapath, 1)
  }
})

observeEvent(input$exampleData, {
  reactiveVals$fcs <- readRDS(file.path(input$exampleData, "fcs.rds"))
  if (input$exampleData=="data/pbmc"){
    reactiveVals$panel <- readRDS(file.path(input$exampleData, "panel.rds"))
    reactiveVals$md <- readRDS(file.path(input$exampleData, "md.rds"))
  } else if (input$exampleData=="data/platelets"){
    reactiveVals$panel <- readRDS(file.path(input$exampleData, "panel.rds"))
    reactiveVals$md <- NULL
  } else if (input$exampleData=="data/mousedata"){
    reactiveVals$md <- NULL
    reactiveVals$panel <- NULL
  }
}, ignoreInit = TRUE)

observeEvent(input$loadData, {
  updateButton(session, "loadData", label = " Loading...", disabled = TRUE)
  library(CATALYST)
  if (input$chooseDataTab == "dataUpload") {
    reactiveVals$sce <- CATALYST::prepData(
      dirname(input$fcsFiles$datapath)[1],
      reactiveVals$panel,
      reactiveVals$md,
      transform = FALSE
      #TODO: check if we have other columns
      #panel_cols = names(reactiveVals$panel),
      #md_cols = names(reactiveVals$md)
    )
  } else if (input$chooseDataTab == "dataExample") {
    reactiveVals$sce <-
      readRDS(file.path(input$exampleData, "sce.rds"))
  } else
    stop("Which tab is selected?")
  updateButton(session, "loadData", label = " Load Data", disabled = FALSE)
  updateButton(session, "continue", label = " Preprocessing")
  shinyjs::show("continue")
  runjs("document.getElementById('continue').scrollIntoView();")
})

output$currentData <- renderInfoBox({
  status <- "warning"
  value <-
    list(
      renderTable(
        checkNullTable(reactiveVals$fcs),
        caption = "FCS Data",
        caption.placement = "top"
      ),
      renderTable(
        checkNullTable(reactiveVals$panel),
        caption = "Panel Data",
        caption.placement = "top"
      ),
      renderTable(
        checkNullTable(reactiveVals$md),
        caption = "Metadata",
        caption.placement = "top"
      )
    )
  
  if (input$chooseDataTab == "dataUpload" &
      !is.null(input$fcsFiles)) {
    status <- "success"
  }
  else if (input$chooseDataTab == "dataExample" &
           input$exampleData != "") {
    status <- "success"
    
    if (input$exampleData=="data/pbmc"){
      exampleData <- "PBMC"
      infoText <- "Data with PBMCs samples from 6 patients.For each sample, the expression of 10 cell surface and 14 signaling markers was measured before (REF) and upon BCR/FcR-XL stimulation (BCRXL) with B cell receptor/ Fc receptor crosslinking for 30', resulting in a total of 12 samples."
      paper <- "Bodenmiller B, Zunder ER, Finck R, et al. Multiplexed mass cytometry profiling of cellular states perturbed by small-molecule regulators. Nat Biotechnol. 2012;30(9):858-867. doi:10.1038/nbt.2317"
      doiLink <- "https://doi.org/10.1038/nbt.2317"
    } else if (input$exampleData=="data/mousedata"){
      exampleData <- "Mouse"
      infoText <- "Data from 10 replicates of mice bone marrow cells. A total of about 840 000 cells were measured using a CyTOF panel of N = 39 markers. "
      paper <- "Samusik N, Good Z, Spitzer MH, Davis KL, Nolan GP. Automated mapping of phenotype space with single-cell data. Nat Methods. 2016;13(6):493-496. doi:10.1038/nmeth.3863"
      doiLink <- "https://doi.org/10.1038/nmeth.3863"
    } else if (input$exampleData=="data/platelets"){
      exampleData <- "Platelets"
      infoText <- "Data from human platelets of patients with chronic coronary syndrome undergoing different therapy: dual antiplatelet therapy versus triple antiplatelet therapy, before and after platelet activation with 10µm TRAP. There are files of 7 patients with triple therapy and 12 patients with dual therapy (each in two condtions)."
      paper <- ""
      doiLink <- ""
    }
    
    value <-
      c(list(
        #div(
        #sprintf(
         # "Found %s: %s",
         # input$exampleData,
         # file.exists(input$exampleData)
        #)
      #),
      strong(sprintf("Info about %s Data :", exampleData)),
      div(infoText),
      a(href=doiLink,paper)
    
      ),
      value)
    
    
    
  }
  
  if (status == "success") {
    value <- c(list(
      bsButton(
        "loadData",
        "Load Data",
        icon("database"),
        style = "success",
        block = TRUE
      )
    ),
    value)
  }
  
  shinydashboard::box(value, title = "Selected Data", status = status)
})
