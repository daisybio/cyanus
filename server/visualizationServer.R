output$currentDataVis <- renderInfoBox({
  status <- "error"
  req(reactiveVals$useExampleDataVis)
  if (reactiveVals$useExampleDataVis) {
    value <- sprintf("info about %s", input$exampleDataVis)
    status <- "success"
  } else
    stop("This is unexpected. Is Data selected or not?")
  if (status == "success"){
    updateActionButton(session, "continue", label = "Clustering")
    shinyjs::show("continue")
  }
  shinydashboard::box(value, title = "Selected Data", status = status)
  
  req(reactiveVals$visMethod)
  if(reactiveVals$visMethod){
    value <- sprintf("you chose %s", input$selectedVisMethod)
  }
  shinydashboard::box(value, title = "Selected Vis")
})

output$plot1 <- renderPlot({
  req(reactiveVals$useExampleDataVis)
  #if(reactiveVals$useExampleDataVis){
  #  library(CATALYST)
  #  plotCounts(readRDS(input$exampleDataVis), group_by = "sample_id", color_by = NULL)
  #}
})