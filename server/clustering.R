### Helpers ----
source("server/clusterFun.R", local = TRUE)

### Observer ----
observeEvent(input$startClustering, {
  updateButton(session,
               "startClustering",
               label = " Clustering...",
               disabled = TRUE)
  shinyjs::disable("clusteringInputs")
  shinyjs::disable("clusteringOutput")
  shinyjs::disable("clusteringVisualizationSelection")
  shinyjs::disable("sidebar")
  
  reactiveVals$sce <-
    clusterSCE(
      reactiveVals$sce,
      input$clusteringMethod,
      input$assayTypeIn,
      input$featuresIn,
      input$xdim,
      input$ydim,
      input$k
    )
  
  shinyjs::enable("clusteringInputs")
  shinyjs::enable("clusteringOutput")
  shinyjs::enable("clusteringVisualizationSelection")
  shinyjs::enable("sidebar")
  updateButton(session, "startClustering", label = " Start Clustering", disabled = FALSE)
  updateButton(session, "continue", label = " Differential Expression Analysis")
  shinyjs::show("continue")
})

observe({
  if (is.null(input$featuresIn))
    disable("startClustering")
  else
    enable("startClustering")
})

### Renderer ----
output$parameters <- renderUI({
  shinydashboard::box(
    column(
      uiOutput("useFeaturesOut"),
      uiOutput("featuresOut"),
      uiOutput("assayTypeOut"),
      width = 6
    ),
    column(uiOutput("k"),
           uiOutput("xdim"),
           uiOutput("ydim"), width = 6),
    title = "Choose Clustering Parameters",
    width = 10
  )
})

output$useFeaturesOut <- renderUI({
  selectInput(
    "useFeaturesIn",
    label = "Features to choose from",
    choices = c("Marker by Class",
                "Marker by Name")
  )
})

output$featuresOut <- renderUI({
  req(input$useFeaturesIn)
  if (input$useFeaturesIn == "Marker by Class") {
    choices <-
      levels(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
  } else if (input$useFeaturesIn == "Marker by Name") {
    choices <- rownames(reactiveVals$sce)
    names(choices) <- sprintf("%s (%s)", choices, as.character(rowData(reactiveVals$sce)$marker_class))
  } else
    stop("by name or by class?")
  shinyWidgets::pickerInput(
    inputId = "featuresIn",
    label = "Features to use for clustering",
    choices = choices,
    selected = choices[1],
    multiple = TRUE,
    options = list(
      `actions-box` = TRUE,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    )
  )
})

output$assayTypeOut <- renderUI({
  choices <- c("Transformed" = "exprs", "Raw" = "counts")
  choices <- choices[choices %in% assayNames(reactiveVals$sce)]
  selectizeInput(
    "assayTypeIn",
    "Expression Type",
    choices = choices
  )
})

output$k <- renderUI({
  req(input$clusteringMethod != "clusterX", input$featuresIn)
  
  if (input$clusteringMethod == "flowSOM") {
    label = "Maximum Number of Clusters to Evaluate in the Metaclustering"
    value = 20
    maxK = 100
  } else if (input$clusteringMethod == "rphenoGraph") {
    label = "Number of Nearest Neighbours"
    value = 30
    maxK = length(.get_features(reactiveVals$sce, input$featuresIn)) - 1
  } else
    stop("which clustermethod was chosen?")
  sliderInput(
    inputId = "k",
    label = label,
    value = min(value, maxK),
    min = 1,
    max = maxK
  )
})

output$ydim <- renderUI({
  req(input$clusteringMethod == "flowSOM")
  sliderInput(
    inputId = "ydim",
    label = "ydim of the grid size of the self-organizing map",
    value = 10,
    min = 2,
    max = 100
  )
  
})

output$xdim <- renderUI({
  req(input$clusteringMethod == "flowSOM")
  sliderInput(
    inputId = "xdim",
    label = "xdim of the grid size of the self-organizing map",
    value = 10,
    min = 2,
    max = 100
  )
})

output$clusteringVisualizationSelection <- renderUI({
  req(metadata(reactiveVals$sce)$cluster_run)
  
  
  shinydashboard::box(
    selectizeInput(
      "clusteringRuns",
      "Successfull Run",
      names(metadata(reactiveVals$sce)$cluster_run)
    ),
    div(
      bsButton(
        "visualizeClustering",
        "Visualize Clustering",
        icon = icon("palette"),
        style = "success"
      ),
      style = "float: right;"
    ),
    column(tableOutput("clusterRunParams"),
           width = 5),
    column(
      uiOutput("selectClusterCode"),
      withSpinner(uiOutput("clusterSizes")),
      width = 5,
      style = "overflow-x: scroll;"
    ),
    fluidRow(withSpinner(uiOutput("delta_area"))),
    title = "Visualize Clustering Results",
    width = 12
  )
})

output$clusterRunParams <- renderTable({
  
  runjs("document.getElementById('clusteringVisualizationSelection').scrollIntoView();")
  runParams <-
    metadata(reactiveVals$sce)$cluster_run[[input$clusteringRuns]]
  runParams$features <- paste(runParams$features, collapse = ",")
  data.frame(runParams)
},
caption = "Run Parameters",
caption.placement = "top")

output$selectClusterCode <- renderUI({
  req(input$clusteringRuns == "flowSOM")
  
  
  selectInput("clusterCode",
              "Clusters",
              rev(names(cluster_codes(reactiveVals$sce))))
})

output$clusterSizes <- renderTable({
  if (input$clusteringRuns == "flowSOM") {
    req(input$clusterCode)
    res <-
      as.data.frame(table(cluster_ids(reactiveVals$sce, input$clusterCode), useNA =
                            "ifany"))
    
  }
  else
    res <-
      as.data.frame(table(colData(reactiveVals$sce)[[sprintf("%s_id", input$clusteringRuns)]], useNA =
                            "ifany"))
  
  names(res) <- c("Cluster", "Size")
  t(res)
},
caption = "Cluster Sizes",
caption.placement = "top", rownames = TRUE, colnames = FALSE)

output$delta_area <- renderUI({
  req(input$clusteringRuns == "flowSOM")
  
  
  runjs("document.getElementById('clusteringVisualizationSelection').scrollIntoView();")
  shinydashboard::box(renderPlot(CATALYST::delta_area(reactiveVals$sce)),
                      title = "Delta Area",
                      width = 12)
})

output$clusteringOutput <- renderUI({
  req(input$visualizeClustering)
  
  shinydashboard::box(fluidRow(
    shinydashboard::tabBox(
      tabPanel("Densities", withSpinner(plotOutput(
        "clusterExprsPlot"
      ))),
      tabPanel("Heatmap", withSpinner(plotOutput(
        "clusterHeatmapPlot"
      ))),
      title = "Cluster Visualization",
      width = 12
    )
  ),
  title = "Clustering Results",
  width = 12)
})

output$clusterExprsPlot <- renderPlot({
  CATALYST::plotClusterExprs(reactiveVals$sce,
                             k = input$clusterCode)#, features = TODO: implement)
})

output$clusterHeatmapPlot <- renderPlot({
  CATALYST::plotFreqHeatmap(reactiveVals$sce,
                            k = input$clusterCode) # TODO: implement other parameters
})
