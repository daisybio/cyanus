# Preprocessing Server

# transform SingleCellExperiment
transformData <-
  function (sce,
            cf = 5,
            ain = "counts",
            aout = "exprs") {
    y <- assay(sce, ain)
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
    assay(sce, aout, FALSE) <- y
    sce
  }

# marker, sample, patient selection -> if markers, patients, samples are selected -> prepValButton can be clicked
observeEvent({
  input$markerSelection
  input$sampleSelection
  input$patientSelection
}, {
  if ((length(input$markerSelection) != 0) &&
      (length(input$sampleSelection) != 0) &&
      (length(input$patientSelection) != 0)) {
    updateActionButton(session, "continue", label = "Visualization")
    shinyjs::show("continue")
    shinyjs::enable("prepVisButton")
  } else {
    shinyjs::disable("prepVisButton")
  }
})

observe({
  if (reactiveVals$current_tab == 3){
    plotPreprocessing(reactiveVals$sce)
  }
})

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
    choices = as.character(unique(colData(reactiveVals$sce)$sample_id)),
    selected = as.character(unique(colData(reactiveVals$sce)$sample_id)),
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

# render samples box
output$patientsBox <- renderUI({
  pickerInput(
    "patientSelection",
    choices = as.character(unique(colData(reactiveVals$sce)$patient_id)),
    selected = as.character(unique(colData(reactiveVals$sce)$patient_id)),
    label = "Patients",
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  )
})

# if start transformation button is clicked
observeEvent(input$prepButton, {
  # data transformation
  reactiveVals$sce <-
    transformData(
      sce = reactiveVals$sce,
      cf = as.numeric(input$cofactor)
    )
})

# if start visualization button is clicked
observeEvent(input$prepVisButton, {
  markers <- isolate(input$markerSelection)
  samples <- isolate(input$sampleSelection)
  patients <- isolate(input$patientSelection)
  sce <- reactiveVals$sce[, reactiveVals$sce$sample_id %in% samples]
  sce <- sce[, sce$patient_id %in% patients]
  sce <- sce[rownames(sce) %in% markers,]
  plotPreprocessing(sce)
})

# method for plotting all kinds of preprocessing plots
plotPreprocessing <- function(sce) {
  
  groupColorLabelBy <- names(colData(sce))
  possAssays <- assayNames(sce)
  if(all(possAssays == c("counts", "exprs"))){
    possAssays <- c("Raw" = "counts", "Normalized" = "exprs")
  }
  features <- c("all", as.character(unique(rowData(sce)$marker_class)))

  output$designCounts <- renderUI({
    fluidRow(column(
      4,
      dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("countsGroupBy",
                       "Group by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput("countsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "countsProp",
          "Stacked or dodged:",
          c(
            "dodged (total cell counts)" = FALSE,
            "stacked (relative abundance)" = TRUE
          ),
          multiple = F
        ),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "100%",
        tooltip = tooltipOptions(title = "Click to see plot options")
      )
    ),
    column(8, shinycssloaders::withSpinner(
      plotlyOutput("countsPlot", width = "100%", height = "400px")
    )))
  })
  
  output$designMDS <- renderUI({
    fluidRow(column(
      4,
      dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("mdsLabelBy",
                       "Label by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput("mdsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "mdsAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput(
          "mdsFeatures",
          "Features:",
          features,
          multiple = F
        ),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "100%",
        tooltip = tooltipOptions(title = "Click to see plot options")
      )
    ),
    column(8, shinycssloaders::withSpinner(
      plotlyOutput("mdsPlot", width = "100%", height = "400px")
    )))
  })
  
  output$designNRS <- renderUI({
    fluidRow(column(
      4,
      dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("nrsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "nrsAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput(
          "nrsFeatures",
          "Features:",
          features,
          multiple = F
        ),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "100%",
        tooltip = tooltipOptions(title = "Click to see plot options")
      )
    ),
    column(8, shinycssloaders::withSpinner(
      plotlyOutput("nrsPlot", width = "100%", height = "400px")
    )))
  })
  
  output$designExprs <- renderUI({
    fluidRow(column(
      4,
      dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("exprsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "exprsAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput(
          "exprsFeatures",
          "Features:",
          features,
          multiple = F
        ),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "100%",
        tooltip = tooltipOptions(title = "Click to see plot options")
      )
    ),
    column(8, shinycssloaders::withSpinner(
      plotlyOutput("exprsPlot", width = "100%", height = "400px")
    )))
  })
  
  output$designExprsHeatmap <- renderUI({
    fluidRow(column(
      4,
      dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput(
          "exprsHeatmapScale",
          "Scale:",
          c("never", "first", "last"),
          multiple = F
        ),
        selectizeInput(
          "exprsHeatmapAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput(
          "exprsHeatmapFeatures",
          "Features:",
          features,
          multiple = F
        ),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "100%",
        tooltip = tooltipOptions(title = "Click to see plot options")
      )
    ),
    column(8, shinycssloaders::withSpinner(
      plotOutput("exprsHeatmapPlot", width = "100%", height = "400px")
    )))
  })
  
  # render counts plot
  output$countsPlot <- renderPlotly({
    ggplotly(
      CATALYST::plotCounts(
        sce,
        group_by = input$countsGroupBy,
        color_by = input$countsColorBy,
        prop = as.logical(input$countsProp)
      )
    )
  })
  
  # render mds plot
  output$mdsPlot <- renderPlotly({
    feature <- input$mdsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    g <- CATALYST::pbMDS(
      sce,
      label_by = input$mdsLabelBy,
      color_by = input$mdsColorBy,
      features = feature,
      assay = input$mdsAssay,
    )
    ggplotly(g)
  })
  
  # render nrs plot
  output$nrsPlot <- renderPlotly({
    feature <- input$nrsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    g <- CATALYST::plotNRS(
      sce,
      color_by = input$nrsColorBy,
      features = feature,
      assay = input$nrsAssay
    )
    ggplotly(g)
  })
  
  # render exprs plot
  output$exprsPlot <- renderPlotly({
    feature <- input$exprsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    g <- CATALYST::plotExprs(
      sce,
      color_by = input$exprsColorBy,
      features = feature,
      assay = input$exprsAssay
    )
    ggplotly(g)
  })
  
  # render exprs heatmap plot
  output$exprsHeatmapPlot <- renderPlot({
    feature <- input$exprsHeatmapFeatures
    if (feature == "all") {
      feature <- NULL
    }
    CATALYST::plotExprHeatmap(
      sce,
      scale = input$exprsHeatmapScale,
      features = feature,
      assay = input$exprsHeatmapAssay
    )
  })
  
  
}