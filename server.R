server <- function(input, output, session) {
  ### menu ----
  observeEvent(input$continue, {
    reactiveVals$current_tab <- reactiveVals$current_tab + 1
    shinyjs::hide("continue")
  })
  
  # grow menu depending on current tab
  output$sidebar <- renderUI({
    curr_menu <- sidebarMenu(id = "tabs",
                             tabs[1:reactiveVals$current_tab])
    updateTabItems(session, "tabs", tab_ids[reactiveVals$current_tab])
    return(curr_menu)
  })
  
  
  ### start ----
  observeEvent(input$exampleData, {
    if (input$exampleData != ""){
      reactiveVals$useExampleData <- T
      reactiveVals$useUploadedData <- F
    }
  })
  
  observeEvent(input$fcsFiles, {
    if (!is.null(input$fcsFiles)){
      reactiveVals$useExampleData <- F
      reactiveVals$useUploadedData <- T
    }
  })
  
  
  #TODO: make tabs
  #TODO: make optional metadata input box
  #TODO: make panel information 
  
  output$currentData <- renderInfoBox({
    status <- "success"
    if (!reactiveVals$useExampleData & !reactiveVals$useUploadedData) {
      value <- "You have currently no data selected."
      status <- "warning"
    }
    else if (reactiveVals$useUploadedData) {
      library(CATALYST)
      # tmp <- CATALYST::prepData(, transform = FALSE)
      fileTable <- input$fcsFiles
      fileTable <- fileTable[, c("name", "size")]
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
  
}
