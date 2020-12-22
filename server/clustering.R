observeEvent(input$startClustering, {
  if (input$clusteringMethod == "flowSOM") {
    reactiveVals$sce <- cluster(reactiveVals$sce,
                                input$features,
                                input$xdim,
                                input$ydim,
                                input$k)
  }
})

output$k <- renderUI({
  if (input$clusteringMethod == "clusterX")
    return(NULL)
  else if (input$clusteringMethod == "flowSOM") {
    label = "Maximum Number of Clusters to Evaluate in the Metaclustering"
    value = 20
  } else if (input$clusteringMethod == "rphenoGraph") {
    label = "Number of Nearest Neighbours"
    value = 30
  } else
    stop("which clustermethod was chosen?")
  sliderInput(
    inputId = "k",
    label = label,
    value = value,
    min = 2,
    max = 100
  )
})

output$ydim <- renderUI({
  if (input$clusteringMethod %in% c("clusterX", "rphenoGraph"))
    return(NULL)
  else if (input$clusteringMethod == "flowSOM") {
    return(
      sliderInput(
        inputId = "ydim",
        label = "ydim of the grid size of the self-orginizing map",
        value = 10,
        min = 2,
        max = 100
      )
    )
  } else
    stop("which clustermethod was chosen?")
})

output$xdim <- renderUI({
  if (input$clusteringMethod %in% c("clusterX", "rphenoGraph"))
    return(NULL)
  else if (input$clusteringMethod == "flowSOM") {
    return(
      sliderInput(
        inputId = "xdim",
        label = "xdim of the grid size of the self-orginizing map",
        value = 10,
        min = 2,
        max = 100
      )
    )
  } else
    stop("which clustermethod was chosen?")
})

output$useFeatures <- renderUI({
  if (input$clusteringMethod %in% c("clusterX", "rphenoGraph"))
    return(NULL)
  else if (input$clusteringMethod == "flowSOM") {
    return(selectInput(
      "useFeatures",
      label = "Parameters to choose from",
      choices = c("Marker by Class",
                  "Marker by Name",
                  "Choose All"),
      selected = "Marker by Class"
    ))
  } else
    stop("which clustermethod was chosen?")
})

output$features <- renderUI({
  if (input$clusteringMethod %in% c("clusterX", "rphenoGraph"))
    return(NULL)
  else if (input$clusteringMethod == "flowSOM") {
    if (is.null(input$useFeatures) || input$useFeatures == "Choose All") {
      return(NULL)
    } else if (input$useFeatures == "Marker by Class") {
      choices = unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
    } else if (input$useFeatures == "Marker by Name") {
      choices = rownames(reactiveVals$sce)
    }
    return(
      selectInput(
        inputId = "features",
        label = "features to use for clustering",
        choices = choices,
        selected = choices[1],
        multiple = TRUE
      )
    )
    
  } else
    stop("which clustermethod was chosen?")
})

output$dimReduction <- renderUI({
  if (input$clusteringMethod %in% c("flowSOM", "rphenoGraph"))
    return(NULL)
  else if (input$clusteringMethod == "clusterX") {
    return(
      selectInput(
        inputId = "dimReduction",
        label = "Dimensionality reduction method.",
        choices = reducedDimNames(reactiveVals$sce)
      )
    )
  } else
    stop("which clustermethod was chosen?")
})

output$clusterPlot <- renderUI({
  
})

output$parameters <- renderUI({
  (
    shinydashboard::box(
      column(uiOutput("k"),
             uiOutput("xdim"),
             uiOutput("ydim"), width = 6),
      column(
        uiOutput("useFeatures"),
        uiOutput("features"),
        uiOutput("dimReduction"),
        width = 6
      ),
      title = "Choose Clustering Parameters",
      width = 9
    )
  )
})