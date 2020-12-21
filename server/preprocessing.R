# Preprocessing Server


preprocessingData <- reactive({
  fcs <- dirname(input$fcsFiles$datapath)[1]
  #md <- dirname(input$metaFile$datapath)[1]
  #panel <- dirname(input$panelFile$datapath)[1]
  transform <- input$transformation
  cofactor <- as.numeric(input$cofactor)
  
  if (transform=="arcsinh"){
    sce <- CATALYST::prepData(fcs, transform = TRUE, cofactor = cofactor)
  } else if (transform=="log"){
    sce <- CATALYST::prepData(fcs)
  } else {
    sce <- CATALYST::prepData(fcs)
  }
  return (sce)
})


ntext <- eventReactive(input$prepButton, { preprocessingData()})
output$nText <- eventReactive(input$prepButton, {
  withProgress(message = 'Progress indicator', {
      ntext()
    })
})

