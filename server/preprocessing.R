# Preprocessing Server

options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

# function calls the prepData function of the CATALYST package
preprocessingData <- function(fcs, md, panel, transform, cofactor){
  if (transform=="arcsinh"){
    sce <- CATALYST::prepData(fcs, transform = TRUE, cofactor = cofactor)
  } else if (transform=="log"){
    sce <- CATALYST::prepData(fcs)
  } else {
    sce <- CATALYST::prepData(fcs)
  }
  return (sce)
}

observeEvent(input$prepButton, {
  # hide start preprocessing button
  shinyjs::hide("prepButton")
  
  # show box of plots
  show(id="plots", anim=TRUE)
  
  # data preparation
  fcs <- dirname(input$fcsFiles$datapath)[1]
  #md <- dirname(input$metaFile$datapath)[1]
  #panel <- dirname(input$panelFile$datapath)[1]
  transform <- input$transformation
  cofactor <- as.numeric(input$cofactor)
  
  # data preprocessing
  sce <- preprocessingData(fcs=fcs,transform=transform, cofactor=cofactor) # add md and panel
  reactiveVals$sce <- sce
  
  # show markers and samples selection
  show(id="markers", anim=TRUE)
  show(id="samples", anim=TRUE)
  updateSelectizeInput(session, 'markerSelection', choices = names(channels(sce)), server = TRUE)
  updateSelectizeInput(session, 'sampleSelection', choices = unique(sample_ids(sce)), server = TRUE)
})  

output$countsPlot <- renderPlot({
  CATALYST::plotCounts(reactiveVals$sce, group_by="sample_id", color_by="sample_id")
})

output$mdsPlot <- renderPlot({
  CATALYST::pbMDS(reactiveVals$sce,color_by="sample_id")
  
})

output$nrsPlot <- renderPlot({
  CATALYST::plotNRS(reactiveVals$sce,features="state", color_by="sample_id")
})

output$exprsPlot <- renderPlot({
  CATALYST::plotExprHeatmap(reactiveVals$sce)
})

observeEvent(input$markerSelection,{
  if (length(input$markerSelection)!=0){
    updateActionButton(session, "continue", label = "Visualization")
    shinyjs::show("continue")
  }
})





