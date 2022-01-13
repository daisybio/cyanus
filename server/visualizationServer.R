shinyjs::hide("visPlotBox")

# ---------------------------------------------------------------
# For resetting after loading a new object 
# ---------------------------------------------------------------

resetVisualization <- function(){
  shinyjs::hide("visPlotBox")
  reactiveVals$availableDRs <- NULL
  reactiveVals$lastMethod <- NULL
  session$sendCustomMessage(type = "resetValue", message = "selectedVisMethod")
  reactiveVals$lastPlot <- NULL
  reactiveVals$visInfo <- NULL
}

# ---------------------------------------------------------------
# Observers
# ---------------------------------------------------------------

observeEvent(input$visTabs, {
  if (input$visTabs == "expressionTab") {
    disable("runDRButton")
  }
})

observeEvent(input$useFeaturesInVis, {
  req(input$useFeaturesInVis)
  if (input$useFeaturesInVis == "Marker by Class") {
    reactiveVals$useClassesRun <- T
    reactiveVals$useMarkersRun <- F
  } else{
    reactiveVals$useClassesRun <- F
    reactiveVals$useMarkersRun <- T
  }
  req(input$selectedRunMethod)
  req(input$pickedFeaturesVis)
  if (input$pickedFeaturesVis != "" &
      input$selectedRunMethod != "") {
    enable("runDRButton")
  }
})


observeEvent(input$selectedRunMethod, {
  if (!is.null(reactiveVals$useClassesRun)) {
    if (input$selectedRunMethod != "" & reactiveVals$useClassesRun) {
      enable("runDRButton")
    }
  }
  if (!is.null(reactiveVals$useMarkersRun)) {
    if (input$selectedRunMethod != "" & reactiveVals$useMarkersRun) {
      enable("runDRButton")
    }
  }
  if (input$selectedRunMethod == "Isomap") {
    shinyjs::show("isobox")
  } else{
    shinyjs::hide("isobox")
  }
})

observeEvent(input$runDRButton, {
  waiter_show(id = "app",html = tagList(spinner$logo, 
                                        HTML("<br>Computing Dimensionality Reduction...")), 
              color=spinner$color)

  reactiveVals$stopVis <- F
  if (reactiveVals$useClassesRun) {
    if (all(
      unique(
        SummarizedExperiment::rowData(reactiveVals$sce)$marker_class %in% input$pickedFeaturesVis
      )
    )) {
      reactiveVals$classes <- NULL
    } else if (length(input$pickedFeaturesVis) > 1) {
      reactiveVals$classes <-
        unname(unlist(sapply(input$pickedFeaturesVis, function(x) {
          rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) == x]
        })))
    } else{
      reactiveVals$classes <- input$pickedFeaturesVis
    }
    reactiveVals$featuresDR <- reactiveVals$classes
    reactiveVals$markers <- c()
  } else if (reactiveVals$useMarkersRun) {
    reactiveVals$markers <- input$pickedFeaturesVis
    reactiveVals$featuresDR <- reactiveVals$markers
    reactiveVals$classes <- c()
  }
  method <- isolate(input$selectedRunMethod)
  
  if(is.null(reactiveVals$visInfo)){
    reactiveVals$visInfo <- list()
  }
  #error handling
  if(method == "TSNE" & input$nrDimensions > 3){
    showNotification("For t-SNE, not more than 3 dimensions are possible. ", type = "error")
    reactiveVals$stopVis <- T
  }else{
    reactiveVals$stopVis <- F
  }
  if(method == "UMAP"){
    library(uwot)
  }
  if(method == "Isomap"){
    k <- isolate(input$valueGraph)
  }else{
    k <- NULL
  }
  
  nrCells <- isolate(input$nrCellsRun)
  seed <- isolate(input$seedSubsamplingDR)
  features <- isolate(reactiveVals$featuresDR)
  assay <- isolate(input$assayRunSelected)
  scale <- isolate(input$scaleRun)
  nrDims <- isolate(input$nrDimensions)

  reactiveVals$visInfo[[method]] <- data.frame(
    method = method,
    features = toString(features),
    counts = assay,
    cells = nrCells,
    seed = seed,
    scale = scale,
    k = toString(k),
    dimensions = nrDims
  )
  if(!reactiveVals$stopVis){
    set.seed(seed)
    reactiveVals$sce <- runCatalystDR(
      sce = reactiveVals$sce,
      dr_chosen = method,
      cells_chosen = nrCells,
      seed = seed,
      feature_chosen = features,
      assay_chosen = assay,
      scale = scale,
      k = k,
      dimensions = nrDims
    )
  }

  reactiveVals$availableDRs <- reducedDimNames(reactiveVals$sce)
  waiter_hide(id="app")
  if(!reactiveVals$stopVis){
    shinyjs::show("visPlotBox")
    if (length(reactiveVals$availableDRs) == 1) {
      output$visPlot <- renderPlot({
        ggplotObject <- makeDR(sce = isolate(reactiveVals$sce), 
                               dr_chosen = method, 
                               color_chosen = renameColorColumn(names(colData(reactiveVals$sce)), T)[[1]], 
                               facet_chosen = "", 
                               assay_chosen = assay, 
                               scale =T, 
                               dims = c(1, 2))
        return(ggplotObject)
      })
    }
  }
})

observeEvent(input$selectedVisMethod, {
  req(input$selectedVisMethod)
  req(input$plt_color_by)
  if (input$plt_color_by != "" & input$selectedVisMethod != "") {
    enable("startDimRed")
  } else{
    disable("startDimRed")
  }
})

observeEvent(input$plt_color_by, {
  if (input$plt_color_by != "" & input$selectedVisMethod != "") {
    enable("startDimRed")
  } else{
    disable("startDimRed")
  }
})



observeEvent(input$nrCellsRun, {
  if (input$nrCellsRun > min(metadata(reactiveVals$sce)$experiment_info$n_cells)) {
    showNotification(
      "This number is higher than the smallest sample count. Therefore, a different number of cells will be sampled from each sample.",
      type = "warning"
    )
  }
})

observeEvent(input$startDimRed, {
  library(ggplot2)
  library(CATALYST)
  output$visPlot <- renderPlot({
    color <- isolate(input$plt_color_by)
    facet <- isolate(input$plt_facet_by)
    method <- isolate(input$selectedVisMethod)
    reactiveVals$lastMethod <- method
    assay <- isolate(input$assayVisSelected)
    scale <- isolate(input$scaleVis)
    sceObj <- isolate(reactiveVals$sce)
    dim1 <- as.numeric(isolate(input$dimension1))
    dim2 <- as.numeric(isolate(input$dimension2))
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Visualizing Dimensionality Reduction...")), 
                color=spinner$color)
    ggplotObject <-
      makeDR(sceObj, method, color, facet, assay, scale, c(dim1, dim2))
    waiter_hide(id="app")
    reactiveVals$lastPlot <- ggplotObject
    return(ggplotObject)
  })
  output$plotInfo <- renderUI({
    method <- isolate(input$selectedVisMethod)
    value <- renderTable(
      checkNullTable(reactiveVals$visInfo[[method]]), 
      caption = paste(method, "Run Features"),
      caption.placement = "top"
    )
    dropdownButton(
      shinydashboard::box(value, title = "Info", width = 12),
      icon = icon("info-circle"),
      status = "info",
      right = TRUE
    )
  })
})

observeEvent(input$radioButtonsColor, {
  if (input$radioButtonsColor == "expression") {
    shinyjs::show("scaleVis")
  } else{
    shinyjs::hide("scaleVis")
  }
})

# ---------------------------------------------------------------
# Render UI 
# ---------------------------------------------------------------

output$useFeaturesInVis <- renderUI({
  selectInput(
    "useFeaturesInVis",
    label = "Features to choose from",
    choices = c("Marker by Class",
                "Marker by Name")
  )
})

output$markersVis <- renderUI({
  req(input$useFeaturesInVis)
  if (input$useFeaturesInVis == "Marker by Class") {
    classes <-
      unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
    choices <- as.vector(classes)
    selected <- choices[1]
    if ("type" %in% choices) {
      selected <- "type"
    }
  } else if (input$useFeaturesInVis == "Marker by Name") {
    choices <- rownames(reactiveVals$sce)
    names(choices) <-
      sprintf("%s (%s)", choices, as.character(marker_classes(reactiveVals$sce)))
    selected <-
      rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) == "type"]
    choices <-
      sortMarkerNames(choices, as.character(marker_classes(reactiveVals$sce)), first = "type")
  } else
    stop("by name or by class?")
  div(
    shinyWidgets::pickerInput(
      inputId = "pickedFeaturesVis",
      label = span(
        "Features to use for the dimensionality reduction",
        icon("question-circle"),
        id = "pickedFeaturesVisQ"
      ),
      choices = choices,
      selected = selected,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      )
    ),
    bsPopover(
      id = "pickedFeaturesVisQ",
      title = "Compute your dimensionality reduction using all data or just a subclass of markers",
      content = "The dimensionality reduction can be computed using all the data (Marker by Class -> type + state), using just a certain marker class like \"type\" or using a subset of markers (Marker by Name). In order to obtain a meaningful visualization, you should choose the marker class that captures best how dissimilar your cells are."
    )
  )
})

output$dimensionsVis <- renderUI({
  div(
    numericInput(
    "nrDimensions",
    label = span("How many dimensions should be returned?",
      icon("question-circle"),
      id = "dimensionsQ"
    ),
    value = 2,
    min = 2,
    max = ifelse(input$selectedRunMethod == "TSNE", 3, 10),
    step = 1
    ),
    bsPopover(
    id = "dimensionsQ",
    title = "Number of dimensions",
    content = "Dimensionality reduction methods reduce the high dimensional input space (e.g. 8 million cells x 20 markers) to a low dimensional representation. The first component that is returned contains the most information, the second one the second most information and so on. Sometimes it might be useful not just to plot component 1 vs. 2 but also component 2 vs. 3 because one component might e.g. represent batch effects. For t-SNE, max. 3 dimensions are possible."
    )
  )
})



output$runDRparBox <- renderUI({
  vis_methods <-
    c("UMAP", "TSNE", "PCA", "MDS", "DiffusionMap", "Isomap")
  choices <- assayNames(reactiveVals$sce)
  selected <- "counts"
  if (all(choices == c("counts", "exprs"))) {
    choices <- c("Raw" = "counts", "Normalized" = "exprs")
    selected <- "exprs"
  }
  shinydashboard::box(
    column(
      width = 6,
      uiOutput("useFeaturesInVis"),
      uiOutput("markersVis"),
      uiOutput("dimensionsVis")
    ),
    column(
      width = 6,
      selectizeInput(
        "selectedRunMethod",
        "Dimensionality Reduction Method",
        vis_methods,
        selected = "UMAP"
      ),
      selectizeInput(
        "assayRunSelected",
        label = "Do you want to use raw counts or the normalization?",
        choices = choices,
        multiple = FALSE,
        selected = selected
      ),
      numericInput(
        "nrCellsRun",
        label = span(
          sprintf(
            "From how many cells do you want to sample (minimal number of cells in a sample: %s)?",
            min(metadata(reactiveVals$sce)$experiment_info$n_cells)
          ),
          icon("question-circle"),
          id = "cellQ"
        ),
        value = 100,
        min = 0,
        max = 10000,
        step = 100
      ),
      bsPopover(
        id = "cellQ",
        title = "Speed up your computations",
        content = "Specify the maximal number of cells per sample to use for dimension reduction, e.g. 100. Then, from each sample 100 cells are randomly taken in order to compute the dimensionality reduction. If you specify values higher than the number of cells in one sample, all cells from these samples will be taken. For the bigger samples, the number you specified will be taken."
      ),
      numericInput(
        "seedSubsamplingDR",
        label = span(
            "Set a seed for your subsampling so that your results are reproducible",
          icon("question-circle"),
          id = "seedsubsamplingDRQ"
        ),
        value = 42,
        min = 1,
        max = 10000,
        step = 1
      ),
      bsPopover(
        id = "seedsubsamplingDRQ",
        title = "Seed for Subsampling",
        content = "Set a seed for your subsampling so that your results are reproducible. If you run different methods with the same seed, the same cells will be used for the dimensionality reduction."
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
      hidden(div(
        id = "isobox",
        numericInput(
            "valueGraph",
            label = span(
              "Choose k for the KNN construction",
              icon("question-circle"),
              id = "kQ"
            ),
            min = 0,
            value = 5,
            max = 200
        ),
        bsPopover(
          id = "kQ",
          title = "Special parameter for Isomap",
          content = "Isomap approximates a manifold using geodesic distances on a k nearest neighbor graph. You can specify k here, the number of nearest neighbors in the graph. "
        )
      )),
      div(
        bsButton(
          "runDRButton",
          "Run with these parameters",
          icon = icon("mouse-pointer"),
          style = "success"
        ),
        style = "float: right;"
      )
    ),
    width = 12,
    title = "Choose Method and Parameters"
  )
})

output$radioButtonsColorVis <- renderUI({
  radioButtons(
    inputId = "radioButtonsColor",
    label = "Color by: ",
    choices = c("expression", "condition / sample"),
    inline = TRUE,
    selected = "condition / sample"
  )
})

output$radioButtonsScale <- renderUI({
  div(
  radioButtons(
    inputId = "scaleVis",
    label = span("Scale ",  icon("question-circle"), id = "scaleVisQ"),
    choices = c("yes", "no"),
    inline = TRUE
  ),
  bsPopover(
    id = "scaleVisQ",
    title = "For coloring by expression",
    content = HTML("Should the counts / the expression data be scaled between 0 and 1 using lower (1%) and upper (99%) expression quantiles?")
  )
  )
})

output$visPlotBox <- renderUI({
  shinydashboard::box(
    fluidRow(
      column(
        2,
        div(
          uiOutput("methodsVis"),
          uiOutput("radioButtonsColorVis"),
          uiOutput("selectColorBy"),
          uiOutput("assayVis"),
          uiOutput("radioButtonsScale"),
          uiOutput("facet_by"),
          uiOutput("plotDimensions"),
          div(
            bsButton(
              "startDimRed",
              "Start Visualization",
              icon = icon("palette"),
              style = "success", 
              disabled = TRUE
            ),
            style = "float: right;",  
            id = "divStartDimRed"
          ),
          style = "position: relative; height: 500px;"
        ),
      ),
      column(
        9,
        shinycssloaders::withSpinner(plotOutput("visPlot", width = "90%")), 
      ),
      column(
        1,
        uiOutput("plotInfo")
      ), 
      div(
        downloadButton("downloadPlot", "Download Plot"),
        style = "position: absolute; bottom: 10px;right:10px"
      )
    ),
    title = "Dimensionality Reduction",
    width = 12
  )
})


output$methodsVis <- renderUI({
  vis_methods <- reactiveVals$availableDRs
  div(
    selectizeInput(
      inputId = "selectedVisMethod",
      label = span(
        "Dimensionality Reduction Method",
        icon("question-circle"),
        id = "drVisQ"
      ),
      vis_methods
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
  selected <- "counts"
  if (all(choices == c("counts", "exprs"))) {
    choices <- c("Raw" = "counts", "Normalized" = "exprs")
    selected <- "exprs"
  }
  return(div(
    selectizeInput(
      "assayVisSelected",
      label = span(
        "Do you want to use raw counts or the normalization?",
        icon("question-circle"),
        id = "assayQ"
      ),
      choices = choices,
      multiple = FALSE,
      selected = selected
    ),
    bsPopover(
      id = "assayQ",
      title = "For coloring by expression",
      content = "When you color by marker expression, the colors will either correspond to the raw counts or the normalized counts."
    )
  ))
})


output$selectColorBy <- renderUI({
  req(input$radioButtonsColor)
  if (input$radioButtonsColor == "expression") {
    choices = c(rownames(reactiveVals$sce))
    return(
      pickerInput(
        inputId = "plt_color_by",
        label = "Color by: ",
        choices = choices,
        selected = choices[1],
        multiple = TRUE
      )
    )
  } else{
    choices = renameColorColumn(names(colData(reactiveVals$sce)), T)
    return(
      selectizeInput(
        inputId = "plt_color_by",
        label = "Color by: ",
        choices = choices,
        multiple = FALSE
      )
    )
  }
})

output$facet_by <- renderUI({
  choices = renameColorColumn(names(colData(reactiveVals$sce)), F)
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
})

output$plotDimensions <- renderUI({
  req(input$selectedVisMethod)
  req(reactiveVals$availableDRs)
  nrDimensionsPossible <- ncol(reducedDim(reactiveVals$sce, input$selectedVisMethod))
  options <- seq(nrDimensionsPossible)
  div(
    selectInput(
      "dimension1",
      "Dimension 1",
      choices = options,
      selected = options[1],
      multiple = F
    ),
    selectInput(
      "dimension2",
      "Dimension 2",
      choices = options,
      selected = options[2],
      multiple = F
    )
  )
})

output$downloadPlot <- downloadHandler(
  filename = function(){
    paste0(reactiveVals$lastMethod, ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$lastPlot, width=14, height=11)
    waiter_hide(id="app")
  }
)


