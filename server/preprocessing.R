# Preprocessing Server

# transform SingleCellExperiment
transformData <-
  function (sce,
            transform,
            cf = 5,
            ain = "counts",
            aout = "exprs") {
    
    if (transform != "no") {
      y <- assay(sce, ain)
      if (transform == "arcsinh") {
        chs <- channels(sce)
        stopifnot(is.numeric(cf), cf > 0)
        if (length(cf) == 1) {
          int_metadata(sce)$cofactor <- cf
          cf <- rep(cf, nrow(sce))
        }
        else {
          stopifnot(!is.null(names(cf)), chs %in% names(cf))
          cf <- cf[match(chs, names(cf))]
          int_metadata(sce)$cofactor <- cf
        }
        fun <- asinh
        op <- "/"
        y <- fun(sweep(y, 1, cf, op))
      } else if (transform == "log") {
        y <- log(y + 1)
      }
      assay(sce, aout, FALSE) <- y
    } else {
      y <- assay(sce, ain)
      assay(sce,aout,FALSE) <- y
    }
    sce
  }

# render markers box
output$markersBox <- renderUI({
  pickerInput(
    inputId = "markerSelection",
    label = "Markers",
    choices = names(channels(reactiveVals$sce)),
    selected = names(channels(reactiveVals$sce)),
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  )
})

# render samples box
output$samplesBox <- renderUI({
  pickerInput(
    "sampleSelection",
    choices = unique(colData(reactiveVals$sce)$sample_id),
    selected = unique(colData(reactiveVals$sce)$sample_id),
    label = "Samples",
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  )
})

# render counts plot
output$countsPlot <- renderPlot({
  CATALYST::plotCounts(
    reactiveVals$sce,
    group_by = input$countsGroupBy,
    color_by = input$countsColorBy,
    prop = as.logical(input$countsProp)
  )
})

# render mds plot
output$mdsPlot <- renderPlot({
  feature <- input$mdsFeatures
  if (feature=="all"){
    feature <- NULL
  }
  CATALYST::pbMDS(
    reactiveVals$sce,
    label_by = input$mdsLabelBy,
    color_by = input$mdsColorBy,
    features = feature,
    assay = input$mdsAssay,
  )
})

# render nrs plot
output$nrsPlot <- renderPlot({
  feature <- input$nrsFeatures
  if (feature=="all"){
    feature <- NULL
  }
  CATALYST::plotNRS(
    reactiveVals$sce,
    color_by = input$nrsColorBy,
    features = feature,
    assay = input$nrsAssay
  )
})

# render exprs plot
output$exprsPlot <- renderPlot({
  feature <- input$exprsFeatures
  if (feature=="all"){
    feature <- NULL
  }
  CATALYST::plotExprs(
    reactiveVals$sce,
    color_by = input$exprsColorBy,
    features = feature,
    assay = input$exprsAssay
  )
})

# render exprs heatmap plot
output$exprsHeatmapPlot <- renderPlot({
  feature <- input$exprsHeatmapFeatures
  if (feature=="all"){
    feature <- NULL
  }
  CATALYST::plotExprHeatmap(
    reactiveVals$sce,
    scale = input$exprsHeatmapScale,
    features = feature,
    assay = input$exprsHeatmapAssay
  )
})

# marker selection -> if markers are selected -> visualization can start
observeEvent({
  input$markerSelection
  input$sampleSelection
}, {
  if ((length(input$markerSelection) != 0) &&
      (length(input$sampleSelection) != 0)) {
    updateActionButton(session, "continue", label = "Visualization")
    shinyjs::show("continue")
  }
})

# if start transformation button is clicked
observeEvent(input$prepButton, {
  # data transformation
  reactiveVals$sce <-transformData(sce = reactiveVals$sce, transform = input$transformation, cf = as.numeric(input$cofactor))
  print(assayNames(reactiveVals$sce))
})

output$designCounts <- renderUI({
  fluidRow(
    column(4,
    dropdownButton(
      tags$h3("Plot Options"),
      selectizeInput("countsGroupBy",
                     "Group by:",
                     names(colData(reactiveVals$sce)), multiple = F),
      selectizeInput("countsColorBy",
                     "Color by:",
                     names(colData(reactiveVals$sce)), multiple = F),
      selectizeInput("countsProp",
                     "Stacked or dodged:",
                     c( "dodged (total cell counts)"=FALSE,"stacked (relative abundance)"=TRUE), multiple = F),
      circle = TRUE,
      status = "info",
      icon = icon("gear"),
      width = "100%",
      tooltip = tooltipOptions(title="Click to see plot options")
    )),
    column(8, shinycssloaders::withSpinner(plotOutput(
      "countsPlot", width = "100%", height = "400px"
    )))
  )
})

output$designMDS <- renderUI({
  fluidRow(
    column(4,
           dropdownButton(
             tags$h3("Plot Options"),
             selectizeInput("mdsLabelBy",
                            "Label by:",
                            names(colData(reactiveVals$sce)), multiple = F),
             selectizeInput("mdsColorBy",
                            "Color by:",
                            names(colData(reactiveVals$sce)), multiple = F),
             selectizeInput("mdsAssay",
                            "Assay:",
                            assayNames(reactiveVals$sce), multiple = F),
             selectizeInput("mdsFeatures",
                            "Features:",
                            c("all",as.character(unique(rowData(reactiveVals$sce)$marker_class))), multiple = F),
             circle = TRUE,
             status = "info",
             icon = icon("gear"),
             width = "100%",
             tooltip = tooltipOptions(title="Click to see plot options")
           )),
    column(8, shinycssloaders::withSpinner(plotOutput(
      "mdsPlot", width = "100%", height = "400px"
    )))
  )
})

output$designNRS <- renderUI({
  fluidRow(
    column(4,
           dropdownButton(
             tags$h3("Plot Options"),
             selectizeInput("nrsColorBy",
                            "Color by:",
                            names(colData(reactiveVals$sce)), multiple = F),
             selectizeInput("nrsAssay",
                            "Assay:",
                            assayNames(reactiveVals$sce), multiple = F),
             selectizeInput("nrsFeatures",
                            "Features:",
                            c("all",as.character(unique(rowData(reactiveVals$sce)$marker_class))), multiple = F),
             circle = TRUE,
             status = "info",
             icon = icon("gear"),
             width = "100%",
             tooltip = tooltipOptions(title="Click to see plot options")
           )),
    column(8, shinycssloaders::withSpinner(plotOutput(
      "nrsPlot", width = "100%", height = "400px"
    )))
  )
})

output$designExprs <- renderUI({
  fluidRow(
    column(4,
           dropdownButton(
             tags$h3("Plot Options"),
             selectizeInput("exprsColorBy",
                            "Color by:",
                            names(colData(reactiveVals$sce)), multiple = F),
             selectizeInput("exprsAssay",
                            "Assay:",
                            assayNames(reactiveVals$sce), multiple = F),
             selectizeInput("exprsFeatures",
                            "Features:",
                            c("all",as.character(unique(rowData(reactiveVals$sce)$marker_class))), multiple = F),
             circle = TRUE,
             status = "info",
             icon = icon("gear"),
             width = "100%",
             tooltip = tooltipOptions(title="Click to see plot options")
           )),
    column(8, shinycssloaders::withSpinner(plotOutput(
      "exprsPlot", width = "100%", height = "400px"
    )))
  )
})

output$designExprsHeatmap <- renderUI({
  fluidRow(
    column(4,
           dropdownButton(
             tags$h3("Plot Options"),
             selectizeInput("exprsHeatmapScale",
                            "Scale:",
                            c("never","first","last"), multiple = F),
             selectizeInput("exprsHeatmapAssay",
                            "Assay:",
                            assayNames(reactiveVals$sce), multiple = F),
             selectizeInput("exprsHeatmapFeatures",
                            "Features:",
                            c("all",as.character(unique(rowData(reactiveVals$sce)$marker_class))), multiple = F),
             circle = TRUE,
             status = "info",
             icon = icon("gear"),
             width = "100%",
             tooltip = tooltipOptions(title="Click to see plot options")
           )),
    column(8, shinycssloaders::withSpinner(plotOutput(
      "exprsHeatmapPlot", width = "100%", height = "400px"
    )))
  )
})

