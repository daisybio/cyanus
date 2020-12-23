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
  reactiveVals$md <- read.table(input$metaFile$datapath, header = T)
})

observeEvent(input$panelFile, {
  reactiveVals$panel <-
    read.table(input$panelFile$datapath, header = T)
})

observeEvent(input$exampleData, {
  reactiveVals$fcs <- readRDS(file.path(input$exampleData, "fcs.rds"))
  # reactiveVals$panel <- file.path(input$exampleData, "panel.rds")
  # reactiveVals$md <- file.path(input$exampleData, "md.rds")
}, ignoreInit = TRUE)

observeEvent(input$loadData, {
  updateButton(session, "loadData", label = " Loading...", disabled = TRUE)
  library(CATALYST)
  #TODO: check if data was uploaded or example selected
  if (input$chooseDataTab == "dataUpload") {
    CATALYST::prepData(
      input$fcsFiles$datapath[1],
      reactiveVals$panel,
      reactiveVals$md,
      transform = FALSE,
      panel_cols = names(reactiveVals$panel),
      md_cols = names(reactiveVals$md)
    )
  } else if (input$chooseDataTab == "dataExample") {
    reactiveVals$sce <-
      readRDS(file.path(input$exampleData, "sce.rds"))
  } else
    stop("Which tab is selected?")
  updateButton(session, "loadData", label = " Load Data", disabled = FALSE)
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
        checkNullTable(reactiveVals$md),
        caption = "Metadata",
        caption.placement = "top"
      ),
      renderTable(
        checkNullTable(reactiveVals$panel),
        caption = "Panel Data",
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
    value <-
      c(list(div(
        sprintf(
          "Found %s: %s",
          input$exampleData,
          file.exists(input$exampleData)
        )
      ),
      div(
        sprintf("info about %s", input$exampleData)
      )),
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
