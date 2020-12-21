output$fcsInfo <- renderTable({
  if (is.null(input$fcsFiles))
    return(data.frame("Nothing" = ""))
  fileTable <- input$fcsFiles
  fileTable <- fileTable[, c("name", "size")]
  fileTable$size <- sprintf("%.2f MB", fileTable$size / 1000000)
  fileTable
}, caption="FCS Data", caption.placement = "top")

output$metaInfo <- renderTable({
  if (is.null(input$metaFile))
    return(data.frame("Nothing" = ""))
  reactiveVals$md <- read.table(input$metaFile$datapath, header = T)
  reactiveVals$md
}, caption="Metadata", caption.placement = "top")

output$panelInfo <- renderTable({
  if (is.null(input$panelFile))
    return(data.frame("Nothing" = ""))
  reactiveVals$panel <- read.table(input$panelFile$datapath, header = T)
  reactiveVals$panel
}, caption="Panel Data", caption.placement = "top")

output$currentData <- renderInfoBox({
  status <- "warning"
  
  if (input$selectedData == "dataUpload") {
    if (!is.null(input$fcsFiles))
      status <- "success"
      library(CATALYST)
      value <- list(tableOutput("fcsInfo"),
                    tableOutput("metaInfo"),
                    tableOutput("panelInfo"))
  }
  else if (input$selectedData == "dataExample") {
    if (input$exampleData == "")
      value <- div("No Example Data selected...")
    else {
      status <- "success"
      value <- sprintf("info about %s", input$exampleData)
    }
  } else
    stop("This is unexpected. Which Tab is selected?")
  
  if (status == "success") {
    updateActionButton(session, "continue", label = "Preprocessing")
    shinyjs::show("continue")
  } else {
    shinyjs::hide("continue")
  }
  
  shinydashboard::box(value, title = "Selected Data", status = status)
})