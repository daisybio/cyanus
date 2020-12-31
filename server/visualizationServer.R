shinyjs::hide("visPlotBox")
#reactiveVals$sce <- readRDS("data/plateletsSCE1-4.Rds")
observeEvent(input$markersState, {
  shinyjs::hide("visPlotBox")
  reactiveVals$useClassesRun <- F
  reactiveVals$useMarkersRun <- T
  req(input$selectedRunMethod)
  if(input$markersState != "" & input$selectedRunMethod != ""){
    updateButton(session, "runDRButton", disabled = FALSE)
  }
})

observeEvent(input$markersType, {
    shinyjs::hide("visPlotBox")
    reactiveVals$useClassesRun <- F
    reactiveVals$useMarkersRun <- T
    req(input$selectedRunMethod)
    if(input$markersType != "" & input$selectedRunMethod != ""){
      updateButton(session, "runDRButton", disabled = FALSE)
    }
  })

observeEvent(input$markersNone, {
  shinyjs::hide("visPlotBox")
  reactiveVals$useClassesRun <- F
  reactiveVals$useMarkersRun <- T
  req(input$selectedRunMethod)
  if(input$markersNone != "" & input$selectedRunMethod != ""){
    updateButton(session, "runDRButton", disabled = FALSE)
  }
})

observeEvent(input$classes, {
  shinyjs::hide("visPlotBox")
  reactiveVals$useClassesRun <- T
  reactiveVals$useMarkersRun <- F
  req(input$selectedRunMethod)
  if(input$classes != "" & input$selectedRunMethod != ""){
    updateButton(session, "runDRButton", disabled = FALSE)
  }
})

observeEvent(input$selectedRunMethod, {
  if(!is.null(reactiveVals$useClassesRun)){
    if(input$selectedRunMethod != "" & reactiveVals$useClassesRun){
      updateButton(session, "runDRButton", disabled = FALSE)
    }
  }
  if(!is.null(reactiveVals$useMarkersRun)){
    if(input$selectedRunMethod != "" & reactiveVals$useMarkersRun){
      updateButton(session, "runDRButton", disabled = FALSE)
    }
  }
  if(input$selectedRunMethod == "Isomap"){
    shinyjs::show("isobox")
  }else{
    shinyjs::hide("isobox")
  }
})

observeEvent(input$runDRButton, {
  disable("continue")
  disable("runDRButton")
  disable("visBox")
  if(reactiveVals$useClassesRun){
    if(input$classes != "All features"){
      reactiveVals$classes <- input$classes
    }else{
      reactiveVals$classes <- NULL
    }
    reactiveVals$featuresDR <- reactiveVals$classes
    reactiveVals$markers <- c()
  }else if(reactiveVals$useMarkersRun){
    reactiveVals$markers <- c(input$markersState, input$markersType, input$markersNone)
    reactiveVals$featuresDR <- reactiveVals$markers
    reactiveVals$classes <- c()
  }
  
  if(input$selectedRunMethod == "UMAP"){
    library(uwot)
    reactiveVals$umapDF <- data.frame(method = "UMAP",
                         features = toString(reactiveVals$featuresDR),
                         counts = input$assayRunSelected,
                         cells = input$nrCellsRun,
                         scale = input$scaleRun)
    runCatalystDR("UMAP", input$nrCellsRun, reactiveVals$featuresDR, 
                  input$assayRunSelected, input$scaleRun, NULL)
  }else if(input$selectedRunMethod == "T-SNE"){
    reactiveVals$tsneDF <- data.frame(method = "TSNE",
                                      features = toString(reactiveVals$featuresDR),
                                      counts = input$assayRunSelected,
                                      cells = input$nrCellsRun,
                                      scale = input$scaleRun)
    runCatalystDR("TSNE", input$nrCellsRun, reactiveVals$featuresDR, 
                  input$assayRunSelected, input$scaleRun, NULL)
  }else if(input$selectedRunMethod == "PCA"){
    reactiveVals$pcaDF <- data.frame(method = "PCA",
                                      features = toString(reactiveVals$featuresDR),
                                      counts = input$assayRunSelected,
                                      cells = input$nrCellsRun,
                                      scale = input$scaleRun)
    runCatalystDR("PCA", input$nrCellsRun, reactiveVals$featuresDR, 
                  input$assayRunSelected, input$scaleRun, NULL)
  }else if(input$selectedRunMethod == "MDS"){
    reactiveVals$mdsDF <- data.frame(method = "MDS",
                                     features = toString(reactiveVals$featuresDR),
                                     counts = input$assayRunSelected,
                                     cells = input$nrCellsRun,
                                     scale = input$scaleRun)
    runCatalystDR("MDS", input$nrCellsRun, reactiveVals$featuresDR, 
                  input$assayRunSelected, input$scaleRun, NULL)
  }else if(input$selectedRunMethod == "DiffusionMap"){
    reactiveVals$diffMapDF <- data.frame(method = "Diffusion Map",
                                     features = toString(reactiveVals$featuresDR),
                                     counts = input$assayRunSelected,
                                     cells = input$nrCellsRun,
                                     scale = input$scaleRun)
    runCatalystDR("DiffusionMap", input$nrCellsRun, reactiveVals$featuresDR, 
                  input$assayRunSelected, input$scaleRun, NULL)
  }else if(input$selectedRunMethod == "Isomap"){
    reactiveVals$isomapDF <- data.frame(method = "Isomap",
                                        features = toString(reactiveVals$featuresDR),
                                        counts = input$assayRunSelected,
                                        cells = input$nrCellsRun,
                                        scale = input$scaleRun, 
                                        k = input$valueGraph)
    runCatalystDR("Isomap", input$nrCellsRun, reactiveVals$featuresDR, 
                  input$assayRunSelected, input$scaleRun, input$valueGraph)
  }
  reactiveVals$useMarkersRun <- c()
  reactiveVals$useClassesRun <- c()
  enable("visBox")
  disable("startDimRed")
  enable("continue")
  enable("runDRButton")
})

observeEvent(input$selectedVisMethod, {
  shinyjs::hide("visPlotBox")
})

observeEvent(reactiveVals$featuresDR, {
  shinyjs::hide("visPlotBox")
})

observeEvent(input$assayVisSelected, {
  shinyjs::hide("visPlotBox")
})

observeEvent(input$scaleVis, {
  shinyjs::hide("visPlotBox")
})

observeEvent(input$plt_color_by, {
  shinyjs::hide("visPlotBox")
  if(input$plt_color_by != "" & input$selectedVisMethod != ""){
    enable("startDimRed")
  }
})


observeEvent(input$plt_facet_by, {
  shinyjs::hide("visPlotBox")
})

observeEvent(input$nrCells, {
  shinyjs::hide("visPlotBox")
})

observeEvent(input$startDimRed, {
  library(ggplot2)
  library(CATALYST)
  output$visPlot <- renderPlotly({
    color <- isolate(input$plt_color_by)
    facet <- isolate(input$plt_facet_by)
    method <- isolate(input$selectedVisMethod)
    assay <- isolate(input$assayVisSelected)
    scale <- isolate(input$scaleVis)
    sceObj <- isolate(reactiveVals$sce)
    plotData(sceObj, method, color, facet, assay, scale)
  })
  output$plotInfo <- renderUI({
    method <- isolate(input$selectedVisMethod)
    if(method == "UMAP"){
      value <- renderTable(
        checkNullTable(reactiveVals$umapDF),
        caption = "UMAP Run Features",
        caption.placement = "top"
      )
    }else if(method == "TSNE"){
      value <- renderTable(
        checkNullTable(reactiveVals$tsneDF),
        caption = "T-SNE Run Features",
        caption.placement = "top"
      )
    }else if(method == "PCA"){
      value <- renderTable(
        checkNullTable(reactiveVals$pcaDF),
        caption = "PCA Run Features",
        caption.placement = "top"
      )
    }else if(method == "MDS"){
      value <- renderTable(
        checkNullTable(reactiveVals$mdsDF),
        caption = "MDS Run Features",
        caption.placement = "top"
      )
    }else if(method == "DiffusionMap"){
      value <- renderTable(
        checkNullTable(reactiveVals$diffMapDF),
        caption = "Diffusion Map Run Features",
        caption.placement = "top"
      )
    }else if(method == "Isomap"){
      value <- renderTable(
        checkNullTable(reactiveVals$isomapDF),
        caption = "Isomap Run Features",
        caption.placement = "top"
      )
    }else{
      value <- "failed"
    }
    shinydashboard::box(value, title = "Info", width = 4)
  })
  shinyjs::show("visPlotBox")
  updateActionButton(session, "continue", label = "Clustering")
  shinyjs::show("continue")
})

plotData <- function(sceObj, method, color, facet, assay, scale){
  disable("startDimRed")
  disable("continue")
  disable("runDRButton")
  
  g <- makeDR(sceObj, method, color, facet, assay, scale)
  g <- g + theme(plot.margin = unit(c(1, 1, 1, 2), "cm"))
  g <- ggplotly(g)
  enable("startDimRed")
  enable("continue")
  enable("runDRButton")
  return(g)
}

observeEvent(input$radioButtonsColor, {
  if(input$radioButtonsColor == "expression"){
    shinyjs::show("assayVis")
    shinyjs::show("scaleVis")
  }else{
    shinyjs::hide("assayVis")
    shinyjs::hide("scaleVis")
  }
})

output$markerClassVis <- renderUI({
  classes <- unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
  choices <- c(as.vector(classes), "All features")
  fluidRow(
    shinydashboard::box(
      selectizeInput(
        inputId = "classes",
        label = "Available marker classes",
        choices = choices,
        options = list(
          placeholder = "Select your marker class",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = FALSE
      ), 
      id = "classBoxVis",
      title = span("Select Class", icon("question-circle"), id = "classQ"), 
      width = 12), 
    bsPopover(
      id = "classQ",
      title = "Compute your dimensionality reduction using all data or just a subclass of markers",
      content = "The dimensionality reduction can be computed using all the data (By Marker Class -> All features), using just a certain marker class like \"type\" or using a subset of markers (By Marker). In order to obtain a meaningful visualization, you should choose the marker class that captures best how dissimilar your cells are."
    )
  )
})

output$expressionVis <- renderUI({
  fluidRow(
    shinydashboard::box(
    pickerInput(
      inputId = "markersState",
      label = "Markers from class state",
      choices = rownames(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == "state"],
      options = list(
        onInitialize = I("function() { this.setValue(''); }"), 
        `actions-box` = TRUE
      ),
      multiple = TRUE
    ), 
    pickerInput(
      inputId = "markersType",
      label = "Markers from class type",
      choices = rownames(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == "type"],
      options = list(
        onInitialize = I("function() { this.setValue(''); }"), 
        `actions-box` = TRUE
      ),
      multiple = TRUE
    ), 
    pickerInput(
      inputId = "markersNone",
      label = "Markers from class none",
      choices = rownames(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == "none"],
      options = list(
        onInitialize = I("function() { this.setValue(''); }"), 
        `actions-box` = TRUE
      ),
      multiple = TRUE
    ), 
    id = "markerBoxVis",
    title = span("Select Markers", icon("question-circle"), id = "markersQ"),
    width = 12
  ), 
  bsPopover(
    id = "markersQ",
    title = "Compute your dimensionality reduction using all data or just a subclass of markers",
    content = "Choose the markers that capture best how dissimilar your cells are. The dimensionality reduction can be computed using all the data (By Marker Class -> All features), using just a certain marker class like \"type\" or using a subset of markers (By Marker)."
  )
  )
})


output$runDRparBox <- renderUI({
  vis_methods <- c("UMAP", "T-SNE", "PCA", "MDS", "DiffusionMap", "Isomap")
  choices <- assayNames(reactiveVals$sce)
  
  if(all(choices == c("counts", "exprs"))){
    choices <- c("Raw" = "counts", "Normalized" = "exprs")
  }
  returnbox <- shinydashboard::box(
    selectizeInput(
      "selectedRunMethod",
      "Dimensionality Reduction Method",
      vis_methods,
      options = list(
        placeholder = "Select your method",
        onInitialize = I("function() { this.setValue(''); }")
      )
    ), 
    selectizeInput(
      "assayRunSelected",
      label = "Do you want to use raw counts or the normalization?",
      choices = choices,
      multiple = FALSE
    ),
    numericInput(
      "nrCellsRun",
      label = span(sprintf("From how many cells do you want to sample (all: %s)?", dim(reactiveVals$sce)[2]), icon("question-circle"), id = "cellQ"),
      value = 100,
      min = 0,
      max = dim(reactiveVals$sce)[2],
      step = 100
    ), 
    bsPopover(
      id = "cellQ",
      title = "Speed up your computations",
      content = "Specify the maximal number of cells per sample to use for dimension reduction, e.g. 100. Then, from each sample 100 cells are randomly taken in order to compute the dimensionality reduction. This is considerably faster than computing the dimensionality reduction for all cells available. If you specify high values that exceed the number of cells per sample, all cells will be taken."
    ),
    radioButtons(
      inputId = "scaleRun",
      label = span("Scale ",  icon("question-circle"), id = "scaleQ"),
      choices = c("yes", "no"), 
      inline = TRUE
    ), 
    bsPopover(
      id = "scaleQ",
      title = "Should the expression values be standardized?",
      content = "Run the dimensionality reduction either with standardized or unchanged counts / expression values"
    ),
    div(id = "isobox",
           fluidRow(numericInput(
             "valueGraph",
             label = span("Choose k for the KNN construction", icon("question-circle"), id = "kQ"),
             min = 0,
             value = 5,
             max = 200
           )),
           bsPopover(
             id = "kQ", 
             title = "Special parameter for Isomap",
             content = "Isomap approximates a manifold using geodesic distances on a k nearest neighbor graph. You can specify k here, the number of nearest neighbors in the graph. "
           )),
    div(
      bsButton(
        "runDRButton",
        "Run with these parameters",
        icon = icon("mouse-pointer"),
        style = "success", 
        disabled = TRUE
      ),
      style = "float: right;"
    ),
    width = 12
  )
  return(returnbox)
})

output$methodsVis <- renderUI({
  vis_methods <- reactiveVals$availableDRs
  div(selectizeInput(
    inputId = "selectedVisMethod",
    label = span("Dimensionality Reduction Method", icon("question-circle"), id = "drVisQ"),
    vis_methods,
    options = list(
      placeholder = "Select how you want to visualize your data",
      onInitialize = I("function() { this.setValue(''); }")
    )
  ),
  bsPopover(
    id = "drVisQ", 
    title = "Available Dimensionality Reductions",
    content = "Here, the dimensionality reductions you already ran are displayed. If you run the same dimensionality reduction twice with different parameters, the old one will be overwritten."
  )
  )
})



output$assayVis <- renderUI({
  choices <- assayNames(reactiveVals$sce)
  if(all(choices == c("counts", "exprs"))){
    choices <- c("Raw" = "counts", "Normalized" = "exprs")
  }
  return(div(selectizeInput(
    "assayVisSelected",
    label = span("Do you want to use raw counts or the normalization?", icon("question-circle"), id = "assayQ"),
    choices = choices,
    multiple = FALSE
  ), 
  bsPopover(
    id = "assayQ", 
    title = "For coloring by expression",
    content = "When you color by marker expression, the colors will either correspond to the raw counts or the normalized counts."
  )
  )
  )
})

output$color_by <- renderUI({
  shinydashboard::box(
  radioButtons(
    inputId = "radioButtonsColor",
    label = "Color by: ",
    choices = c("expression", "condition / sample"), 
    inline = TRUE
  ), 
  uiOutput("selectColorBy")
    )
})

output$selectColorBy <- renderUI({
  if(input$radioButtonsColor == "expression"){
    choices = c(rownames(reactiveVals$sce))
    return(pickerInput(
      inputId = "plt_color_by",
      label = "Color by: ",
      choices = choices,
      options = list(
        onInitialize = I("function() { this.setValue(''); }"), 
        `actions-box` = TRUE
      ),
      multiple = TRUE
    ))
    }else{
      choices = names(colData(reactiveVals$sce))
      return(selectizeInput(
        inputId = "plt_color_by",
        label = "Color by: ",
        choices = choices,
        options = list(
          placeholder = "Select your coloring",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = FALSE
      ))
    }
})

output$facet_by <- renderUI({
  choices = names(colData(reactiveVals$sce))
  shinydashboard::box(
  selectizeInput(
    inputId = "plt_facet_by",
    label = "Facet by: ",
    choices = choices,
    options = list(
      placeholder = "Select your faceting or nothing",
      onInitialize = I("function() { this.setValue(''); }")
    ),
    multiple = FALSE
  )
  )
})

runCatalystDR <- function(dr_chosen, cells_chosen, feature_chosen, assay_chosen, scale, k){
  if(scale == "yes"){
    scale <- T
  }else{
    scale <- F
  }
  if(dr_chosen == "Isomap"){
    reactiveVals$sce <- runIsomap(
      reactiveVals$sce, 
      cells = cells_chosen, 
      features = feature_chosen, 
      assay = assay_chosen, 
      scale = scale, 
      k = k)
    
  }else{
  reactiveVals$sce <- runDR(
    reactiveVals$sce, 
    dr = dr_chosen, 
               cells = cells_chosen, 
               features = feature_chosen, 
               assay = assay_chosen, 
               scale = scale)
  }
  reactiveVals$availableDRs <- reducedDimNames(reactiveVals$sce)
}

makeDR <- function(sce, dr_chosen, color_chosen, facet_chosen, assay_chosen, scale) {
  if(color_chosen == ""){
    color_chosen <- NULL
  }
  if(facet_chosen == ""){
    facet_chosen <- NULL
  } 
  if(assay_chosen == ""){
    assay_chosen <- NULL
  }
  if(scale == "yes"){
    scale <- T
  }else{
    scale <- F
  }
  g <- plotDR(sce, dr = dr_chosen, color_by = color_chosen, facet_by = facet_chosen, assay = assay_chosen, scale = scale)
  return(g)
}
