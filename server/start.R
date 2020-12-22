reactiveVals$fcs <- data.frame("Nothing" = "")
reactiveVals$panel <- data.frame("Nothing" = "")
reactiveVals$md <- data.frame("Nothing" = "")

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
  reactiveVals$panel <- read.table(input$panelFile$datapath, header = T)
})

observeEvent(input$exampleData, {
  reactiveVals$fcs <- readRDS(file.path(input$exampleData, "fcs.rds"))
  # reactiveVals$panel <- file.path(input$exampleData, "panel.rds")
  # reactiveVals$md <- file.path(input$exampleData, "md.rds")
}, ignoreInit = TRUE)

#TODO add load data

output$currentData <- renderInfoBox({
  status <- "warning"
  value <- list(renderTable(reactiveVals$fcs, caption = "FCS Data", caption.placement = "top"),
                  renderTable(reactiveVals$md, caption = "Metadata", caption.placement = "top"),
                  renderTable(reactiveVals$panel, caption = "Panel Data", caption.placement = "top"))
  
  if (input$selectedData == "dataUpload" & !is.null(input$fcsFiles)) {
      status <- "success"
  }
  else if (input$selectedData == "dataExample" & input$exampleData != "") {
      status <- "success"
      value <-
        c(div(sprintf(
          "Found %s: %s",
          input$exampleData,
          file.exists(input$exampleData)
        )),
        div(sprintf("info about %s", input$exampleData)),
        value)
      # TODO: remove this as soon as preprocessing is done
      library(CATALYST)
      reactiveVals$sce <- readRDS(file.path(input$exampleData, "sce.rds"))
  } 
  if (status == "success") {
    updateActionButton(session, "continue", label = "Preprocessing")
    shinyjs::show("continue")
  } else {
    shinyjs::hide("continue")
  }
  
  shinydashboard::box(value, title = "Selected Data", status = status)
})