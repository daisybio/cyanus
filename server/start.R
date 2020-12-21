fileOverview <- function(fileTable) {
  fileTable <- fileTable[, c("name", "size")]
  fileTable$size <- sprintf("%.2f MB", fileTable$size/1000000)
  return(fileTable)
}


output$fcsInfo <- renderTable(fileOverview())


output$currentData <- renderInfoBox({
  status <- "success"
  # if (!reactiveVals$useExampleData & !reactiveVals$useUploadedData) {
  #   value <- "You have currently no data selected."
  #   status <- "warning"
  # }
  # else 
  if (input$selectedData == "dataUpload") {
    library(CATALYST)
    # tmp <- CATALYST::prepData(, transform = FALSE)
    fileTable <- input$fcsFiles
    fileTable <- fileTable[, c("name", "size")]
    fileTable
    fileTable$size <- sprintf("%.2f MB", fileTable$size/1000000)
    value <- list(div("The following file(s) have been uploaded:"),
                  renderTable(fileTable)
    )
  }
  else if (reactiveVals$useExampleData) {
    value <- sprintf("info about %s", input$exampleData)
  } else
    stop("This is unexpected. Is Data selected or not?")
  if (status == "success"){
    updateActionButton(session, "continue", label = "Preprocessing")
    shinyjs::show("continue")
  }
  shinydashboard::box(value, title = "Selected Data", status = status)
})