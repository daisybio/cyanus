observeEvent(input$startClustering, {
  updateButton(session,
               "startClustering",
               label = " Clustering...",
               disabled = TRUE)
  shinyjs::disable("clusteringInputs")
  shinyjs::disable("clusteringOutput")
  
  
  # TODO: remove this
  reactiveVals$sce <-
    prepData(
      "extdata/PBMC/PBMC8_fcs_files/",
      readRDS("data/pbmc/panel.rds"),
      readRDS("data/pbmc/md.rds")
    )
  
  if (input$clusteringMethod == "flowSOM") {
    if (input$useFeatures == "Choose All")
      features <- NULL
    else
      features <- isolate(input$featuresIn)
    reactiveVals$sce <- CATALYST::cluster(reactiveVals$sce,
                                          features,
                                          input$xdim,
                                          input$ydim,
                                          input$k)
  }
  shinyjs::enable("clusteringInputs")
  shinyjs::enable("clusteringOutput")
  updateButton(session,
               "startClustering",
               label = " Start Clustering",
               disabled = FALSE)
})

observeEvent(input$plotClusterExpr, {
  #TODO: plot exproutput
  browser()
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
                  "Choose All")
    ))
  } else
    stop("which clustermethod was chosen?")
})

output$featuresOut <- renderUI({
  if (input$clusteringMethod %in% c("clusterX", "rphenoGraph"))
    return(NULL)
  else if (input$clusteringMethod == "flowSOM") {
    if (is.null(input$useFeatures) ||
        input$useFeatures == "Choose All") {
      return(NULL)
    } else if (input$useFeatures == "Marker by Class") {
      choices = unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
    } else if (input$useFeatures == "Marker by Name") {
      choices = rownames(reactiveVals$sce)
    }
    return(
      selectInput(
        inputId = "featuresIn",
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

output$delta_area <- renderPlot({
  if (!"delta_area" %in% names(metadata(reactiveVals$sce)))
    return(NULL)
  
  CATALYST::delta_area(reactiveVals$sce)
})

output$clusterSizes <- renderTable({
  res <-
    as.data.frame(table(cluster_ids(reactiveVals$sce, input$clusterCode)))
  names(res) <- c("ClusterID", "Size")
  res
})

output$clusterExprsPlot <- renderPlot({
  CATALYST::plotClusterExprs(reactiveVals$sce,
                             k = input$clusterCode)#, features = TODO: implement)
})

output$clusterHeatmapPlot <- renderPlot({
  CATALYST::plotClusterHeatmap(reactiveVals$sce,
                               k = input$clusterCode) # TODO: implement other parameters
})

output$clusterExprs <- renderUI({
  fluidRow(
    selectInput("clusterCode",
                "Clusters",
                rev(names(
                  cluster_codes(reactiveVals$sce)
                )),
                selected = input$k),
    uiOutput("clusterSizes"),
    bsButton("plotClusterExpr",
             "Plot Cluster Expression"),
    tabsetPanel(
      tabPanel("Densities", plotOutput("clusterExprsPlot")),
      tabPanel("Heatmap", plotOutput("clusterHeatmapPlot"))
    )
  )
})


output$clusteringOutput <- renderUI({
  if (!"cluster_id" %in% names(colData(reactiveVals$sce)))
    return(NULL)
  
  shinydashboard::box(fluidRow(uiOutput("clusterExprs")),
                      fluidRow(plotOutput("delta_area")),
                      title = "Clustering Visualization",
                      width = 12)
})

output$parameters <- renderUI({
  (
    shinydashboard::box(
      column(uiOutput("k"),
             uiOutput("xdim"),
             uiOutput("ydim"), width = 6),
      column(
        uiOutput("useFeatures"),
        uiOutput("featuresOut"),
        uiOutput("dimReduction"),
        width = 6
      ),
      title = "Choose Clustering Parameters",
      width = 9
    )
  )
})
