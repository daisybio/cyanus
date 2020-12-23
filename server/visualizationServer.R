
observeEvent(input$selectMarkersVis, {
    if (input$selectMarkersVis){
      reactiveVals$sceTMP <- reactiveVals$sce[rownames(reactiveVals$sce) %in% input$markers]
      reactiveVals$featuresDR <- input$markers
      shinyjs::show("visBox")
      shinyjs::hide("visPlotBox")
    }
})

observeEvent(input$selectClassesVis, {
  if(input$selectClassesVis){
    reactiveVals$sceTMP <- reactiveVals$sce[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == input$classes]
    reactiveVals$featuresDR <- input$classes
    shinyjs::show("visBox")
    shinyjs::hide("visPlotBox")
  }
})

observeEvent(input$startDimRed, {
  output$visPlot <- renderPlot({
    library(ggplot2)
    library(CATALYST)
    nr_cells <- input$nrCells
    feat <- reactiveVals$featuresDR
    color <- input$plt_color_by
    if(color == "expression"){
      color <- input$markers
    }
    facet <- input$plt_facet_by
    if(input$selectedVisMethod == "UMAP"){
      makeDR(reactiveVals$sceTMP, "UMAP", feat, nr_cells , input$plt_color_by, input$plt_facet_by)
    }else if(input$selectedVisMethod == "T-SNE"){
      makeDR(reactiveVals$sceTMP, "TSNE", feat, nr_cells , input$plt_color_by, input$plt_facet_by)
    }else if(input$selectedVisMethod == "PCA"){
      makeDR(reactiveVals$sceTMP, "PCA", feat, nr_cells , input$plt_color_by, input$plt_facet_by)
    }
    
  })
  shinyjs::show("visPlotBox")
})

output$currentDataVis <- renderUI({
  status <- "error"
  reactiveVals$sce <- readRDS("data/plateletsSCE1-4.Rds")
  if(input$selectedVisMethod != ""){
    status <- "success"
  }
  if (status == "success"){
    updateActionButton(session, "continue", label = "Clustering")
    shinyjs::show("continue")
  }
})

output$parametersVis <- renderUI({
  shinyjs::hide("visBox")
  shinydashboard::tabBox(
    tabPanel(uiOutput("markerClassVis"), 
             value = "classTab",
             title = "By Marker Class"),
    tabPanel(uiOutput("expressionVis"),
             value = "expressionTab",
             title = "By Marker"),
    id = "visTabs",
    title = "Choose",
    width = 12
  )
})


output$expressionVis <- renderUI({
  shinyjs::hide("visBox")
      fluidRow(shinydashboard::box(
        selectizeInput(
          inputId = "markers",
          label = "markers to use for clustering",
          choices = rownames(reactiveVals$sce),
          multiple = TRUE
        ), 
        div(
          bsButton(
            "selectMarkersVis",
            "Select these markers",
            icon = icon("mouse-pointer"),
            style = "success"
          ),
          style = "float: right;"
        ),
        id = "markerBoxVis",
        title = "Select Markers",
        height = "12em",
        width = 12
      ))
})

output$markerClassVis <- renderUI({
  shinyjs::hide("visBox")
  fluidRow(
    shinydashboard::box(
      selectizeInput(
        inputId = "classes",
        label = "marker classes to use for clustering",
        choices = unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class),
        multiple = FALSE
      ),
      div(
        bsButton(
          "selectClassesVis",
          "Select this class",
          icon = icon("mouse-pointer"),
          style = "success"
        ),
        style = "float: right;"
      ),
      id = "classBoxVis",
      title = "Select Class",
      height = "12em",
      width = 12
    )
  )
})

output$color_by <- renderUI({
  colorVector <- c(names(colData(reactiveVals$sceTMP)), rownames(reactiveVals$sceTMP))
  if(input$selectMarkersVis){
    colorVector <- c("expression", colorVector)
  }
  selectizeInput(
    inputId = "plt_color_by",
    label = "Color by: ",
    colorVector,
    options = list(
      placeholder = "Select your coloring or nothing",
      onInitialize = I("function() { this.setValue(''); }")
    ),
    multiple = FALSE
  )
})

output$facet_by <- renderUI({
  selectizeInput(
    inputId = "plt_facet_by",
    label = "Facet by: ",
    names(colData(reactiveVals$sceTMP)),
    options = list(
      placeholder = "Select your faceting or nothing",
      onInitialize = I("function() { this.setValue(''); }")
    ),
    multiple = FALSE
  )
})

makeDR <- function(dr_chosen, feature_chosen, cell_nr, color_chosen, facet_chosen) {
  sce <- runDR(sce, dr = dr_chosen, features = feature_chosen, cells = cell_nr)
  g <- plotDR(sce, dr = dr_chosen, color_by = color_chosen, facet_by = facet_chosen)
  return(g)
}
