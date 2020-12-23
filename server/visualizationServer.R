observeEvent(input$markersState, {
  shinyjs::hide("visBox")
  shinyjs::hide("visPlotBox")
  if(input$markersState != ""){
    updateButton(session, "selectMarkersVis", disabled = FALSE)
  }
})

observeEvent(input$markersType, {
    shinyjs::hide("visBox")
    shinyjs::hide("visPlotBox")
    if(input$markersType != ""){
      updateButton(session, "selectMarkersVis", disabled = FALSE)
    }
  })

observeEvent(input$markersNone, {
  shinyjs::hide("visBox")
  shinyjs::hide("visPlotBox")
  if(input$markersNone != ""){
    updateButton(session, "selectMarkersVis", disabled = FALSE)
  }
})

observeEvent(input$classes, {
  shinyjs::hide("visBox")
  shinyjs::hide("visPlotBox")
  if(input$classes != ""){
    updateButton(session, "selectClassesVis", disabled = FALSE)
  }
})

observeEvent(input$selectMarkersVis, {
  reactiveVals$markers <- c(input$markersState, input$markersType, input$markersNone)
  reactiveVals$featuresDR <- reactiveVals$markers
  reactiveVals$classes <- c()
  shinyjs::show("visBox")
})

observeEvent(input$selectClassesVis, {
  reactiveVals$classes <- input$classes
  reactiveVals$featuresDR <- reactiveVals$classes
  reactiveVals$markers <- c()
  shinyjs::show("visBox")
})

observeEvent(input$selectedVisMethod, {
  shinyjs::hide("visPlotBox")
})

observeEvent(reactiveVals$featuresDR, {
  shinyjs::hide("visPlotBox")
})

observeEvent(input$plt_color_by, {
  shinyjs::hide("visPlotBox")
  if(input$plt_color_by != ""){
    updateButton(session, "startDimRed", disabled = FALSE)
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
  output$visPlot <- renderPlot({
    plotData()
  })
  shinyjs::show("visPlotBox")
  updateActionButton(session, "continue", label = "Clustering")
  shinyjs::show("continue")
})

plotData <- eventReactive(input$startDimRed, {
  disable("startDimRed")
  disable("continue")
  disable("selectClassesVis")
  disable("selectMarkersVis")
  nr_cells <- isolate(input$nrCells)
  feat <- isolate(reactiveVals$featuresDR)
  color <- isolate(input$plt_color_by)
  facet <- isolate(input$plt_facet_by)
  method <- isolate(input$selectedVisMethod)
  sceObj <- isolate(reactiveVals$sce)
  
  if(method == "UMAP"){
    g <- makeDR(sceObj, "UMAP", feat, nr_cells , color, facet)
  }else if(method == "T-SNE"){
    g <- makeDR(sceObj, "TSNE", feat, nr_cells , color, facet)
  }else if(method == "PCA"){
    g <- makeDR(sceObj, "PCA", feat, nr_cells , color, facet)
  }else{
    g <- NULL
  }
  enable("startDimRed")
  enable("continue")
  enable("selectClassesVis")
  enable("selectMarkersVis")
  return(g)
})

output$parametersVis <- renderUI({
  shinydashboard::tabBox(
    tabPanel(uiOutput("markerClassVis"), 
             value = "classTab",
             title = "By Marker Class"),
    tabPanel(uiOutput("expressionVis"),
             value = "expressionTab",
             title = "By Marker"),
    id = "visTabs",
    title = "Choose which antigens to use",
    width = 6
  )
})

output$markerClassVis <- renderUI({
  fluidRow(
    shinydashboard::box(
      selectizeInput(
        inputId = "classes",
        label = "Available marker classes",
        unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class),
        options = list(
          placeholder = "Select your marker class",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = FALSE
      ),
      div(
        bsButton(
          "selectClassesVis",
          "Select this class",
          icon = icon("mouse-pointer"),
          style = "success", 
          disabled = TRUE
        ),
        style = "float: right;"
      ), 
      id = "classBoxVis",
      title = "Select Class",
      height = "12em")
  )
})

output$expressionVis <- renderUI({
      fluidRow(shinydashboard::box(
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
        div(
          bsButton(
            "selectMarkersVis",
            "Select these markers",
            icon = icon("mouse-pointer"),
            style = "success", 
            disabled = TRUE
          ),
          style = "float: right;"
        ),
        id = "markerBoxVis",
        title = "Select Markers",
        height = "12em"
      ))
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
    if(!is.null(reactiveVals$markers)){
      choices = c(rownames(reactiveVals$sce)[rownames(reactiveVals$sce) %in% reactiveVals$markers])
    }else if(!is.null(reactiveVals$classes)){
      choices = c(rownames(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == reactiveVals$classes])
    }else{
      choices = c("Something went wrong")
    }
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
      if(!is.null(reactiveVals$markers)){
        choices = names(colData(reactiveVals$sce))
      }else if(!is.null(reactiveVals$classes)){
        choices = names(colData(reactiveVals$sce))
      }else{
        choices = c("Something went wrong")
      }
      return(selectizeInput(
        inputId = "plt_color_by",
        label = "Color by: ",
        choices = choices,
        options = list(
          placeholder = "Select your coloring or nothing",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = FALSE
      ))
    }
})

output$facet_by <- renderUI({
  if(!is.null(reactiveVals$markers)){
    choices = names(colData(reactiveVals$sce))
  }else if(!is.null(reactiveVals$classes)){
    choices = names(colData(reactiveVals$sce))
  }else{
    choices = c("Something went wrong")
  }
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

makeDR <- function(sce, dr_chosen, feature_chosen, cell_nr, color_chosen, facet_chosen) {
  if(color_chosen == ""){
    color_chosen <- NULL
  }
  if(facet_chosen == ""){
    facet_chosen <- NULL
  } 
  sce <- runDR(sce, dr = dr_chosen, features = feature_chosen, cells = cell_nr)
  g <- plotDR(sce, dr = dr_chosen, color_by = color_chosen, facet_by = facet_chosen)
  return(g)
}
