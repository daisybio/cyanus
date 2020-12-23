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
  #reactiveVals$sceTMP <- reactiveVals$sce[rownames(reactiveVals$sce) %in% reactiveVals$markers]
  reactiveVals$featuresDR <- reactiveVals$markers
  input$classes <- ""
  shinyjs::show("visBox")
})

observeEvent(input$selectClassesVis, {
  #reactiveVals$sceTMP <- reactiveVals$sce[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == input$classes]
  reactiveVals$featuresDR <- input$classes
  reactiveVals$markers <- ""
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
    nr_cells <- isolate(input$nrCells)
    feat <- isolate(reactiveVals$featuresDR)
    color <- isolate(input$plt_color_by)
    facet <- isolate(input$plt_facet_by)
    method <- isolate(input$selectedVisMethod)
    #sceObj <- isolate(reactiveVals$sceTMP)
    sceObj <- isolate(reactiveVals$sce)
    
    if(method == "UMAP"){
      makeDR(sceObj, "UMAP", feat, nr_cells , color, facet)
    }else if(method == "T-SNE"){
      makeDR(sceObj, "TSNE", feat, nr_cells , color, facet)
    }else if(method == "PCA"){
      makeDR(sceObj, "PCA", feat, nr_cells , color, facet)
    }
  })
  shinyjs::show("visPlotBox")
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
    title = "Choose",
    width = 12
  )
})

output$markerClassVis <- renderUI({
  fluidRow(
    shinydashboard::box(
      selectizeInput(
        inputId = "classes",
        label = "marker classes to use for clustering",
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
      height = "12em",
      width = 12)
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
        height = "12em",
        width = 12
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
    if(reactiveVals$markers != ""){
      choices = c(rownames(reactiveVals$sce)[rownames(reactiveVals$sce) %in% reactiveVals$markers])
    }else if(input$classes != ""){
      choices = c(rownames(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == input$classes])
    }else{
      choices = c("Something went wrong")
    }
    return(pickerInput(
      inputId = "plt_color_by",
      label = "Color by: ",
      #choices = c(rownames(reactiveVals$sceTMP)),
      choices = choices,
      options = list(
        onInitialize = I("function() { this.setValue(''); }"), 
        `actions-box` = TRUE
      ),
      multiple = TRUE
    ))
    }else{
      if(reactiveVals$markers != ""){
        choices = c(colData(reactiveVals$sce)[rownames(reactiveVals$sce) %in% reactiveVals$markers])
      }else if(input$classes != ""){
        choices = c(colData(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == input$classes])
      }else{
        choices = c("Something went wrong")
      }
      return(selectizeInput(
        inputId = "plt_color_by",
        label = "Color by: ",
        #names(colData(reactiveVals$sceTMP)),
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
  if(reactiveVals$markers != ""){
    choices = c(colData(reactiveVals$sce)[rownames(reactiveVals$sce) %in% reactiveVals$markers])
  }else if(input$classes != ""){
    choices = c(colData(reactiveVals$sce)[SummarizedExperiment::rowData(reactiveVals$sce)$marker_class == input$classes])
  }else{
    choices = c("Something went wrong")
  }
  shinydashboard::box(
  selectizeInput(
    inputId = "plt_facet_by",
    label = "Facet by: ",
    #names(colData(reactiveVals$sceTMP)),
    choices = choices,
    options = list(
      placeholder = "Select your faceting or nothing",
      onInitialize = I("function() { this.setValue(''); }")
    ),
    multiple = FALSE
  )
  )
})

output$currentDataVis <- renderUI({
  status <- "error"
  reactiveVals$sce <- readRDS("data/plateletsSCE1-4.Rds")
  if(input$startDimRed){
    status <- "success"
  }
  if (status == "success"){
    updateActionButton(session, "continue", label = "Clustering")
    shinyjs::show("continue")
  }
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
