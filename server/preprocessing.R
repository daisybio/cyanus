# Preprocessing Server

observeEvent(input$prepData,{
  # construct SingleCellExperiment
  fcs <- dirname(input$fcsFiles$datapath)[1]
  md <- dirname(input$metaFile$datapath)[1]
  panel <- dirname(input$panelFile$datapath)[1]
  
  if (input$transformation == "no"){
    sce <- CATALYST::prepData(filepath)
    
  } else if (input$transformation == "log"){
    sce <- CATALYST:prepData(filepath)
  } else {
    cofactor <- as.numeric(input$cofactor)
    print(cofactor)
    sce <- CATALYST::prepData(filepath, transform = TRUE, cofactor = cofactor)
  }
})

output$plotCounts <- renderPlot({
  plotCounts(sce, group_by="sample_id", color_by=NULL)
})