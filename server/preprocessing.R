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
}, {
  if (length(input$sampleSelection) != 0){
    updateActionButton(session, "continue", label = "Visualization")
    shinyjs::show("continue")
    shinyjs::enable("prepSelectionButton")
    shinyjs::enable("filterSelectionButton")
  } else {
    shinyjs::disable("prepSelectionButton")
    shinyjs::disable("filterSelectionButton")
  }
}, ignoreNULL = FALSE)

# check current tab
observe({
  if (reactiveVals$current_tab == 3) {
    plotPreprocessing(reactiveVals$sce)
    if (!("patient_id" %in% colnames(colData(reactiveVals$sce)))){
      shinyjs::hide("patientsBox")
    }
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
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("filterSelectionButton")
  shinyjs::disable("continue")
  # data transformation
  reactiveVals$sce <-
    transformData(sce = reactiveVals$sce,
                  cf = as.numeric(input$cofactor))
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("filterSelectionButton")
  shinyjs::enable("continue")
})

# if visualize selection button is clicked
observeEvent(input$prepSelectionButton, {
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("continue")
  shinyjs::disable("filterSelectionButton")
  allpatients <- length(as.character(unique(colData(reactiveVals$sce)$patient_id)))
  allsamples <- length(as.character(unique(colData(reactiveVals$sce)$sample_id)))
  if ((length(input$patientSelection) != allpatients) || (length(input$sampleSelection) != allsamples)){
    showNotification(HTML(
      "<b>Attention!</b><br>
      The unselected samples and patients are <b>deleted</b> from the data when pressing the <b>Confirm Selection</b> button. Further analysis is being performed only on the selected patients and samples!"
    ),
    duration = 10,
    type = "warning")
  }
  markers <- isolate(input$markerSelection)
  samples <- isolate(input$sampleSelection)
  patients <- isolate(input$patientSelection)
  sce <- filterSCE(reactiveVals$sce,sample_id %in% samples)
  if (("patient_id" %in% colnames(colData(reactiveVals$sce)))){
    sce <- filterSCE(sce,patient_id %in% patients)
  }
  
  sce <- sce[rownames(sce) %in% markers, ]
  plotPreprocessing(sce)
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("continue")
  shinyjs::enable("filterSelectionButton")
})

# if filtering button is clicked -> selection is applied to sce
observeEvent(input$filterSelectionButton,{
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("filterSelectionButton")
  shinyjs::disable("continue")
  allpatients <- length(as.character(unique(colData(reactiveVals$sce)$patient_id)))
  allsamples <- length(as.character(unique(colData(reactiveVals$sce)$sample_id)))
  
  if (length(input$sampleSelection) != allsamples){
    reactiveVals$sce <- filterSCE(reactiveVals$sce,sample_id %in% input$sampleSelection)
  }
  if (("patient_id" %in% colnames(colData(reactiveVals$sce)))){
    if (length(input$patientSelection) != allpatients){
      reactiveVals$sce <- filterSCE(reactiveVals$sce,patient_id %in% input$patientSelection)
    }
  }
  
  markers <- isolate(input$markerSelection)
  sce <- reactiveVals$sce[rownames(reactiveVals$sce) %in% markers, ]
  plotPreprocessing(sce)
  
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("filterSelectionButton")
  shinyjs::enable("continue")
})

# method for plotting all kinds of preprocessing plots
plotPreprocessing <- function(sce) {
  groupColorLabelBy <- names(colData(sce))
  possAssays <- assayNames(sce)
  if (all(possAssays == c("counts", "exprs"))) {
    possAssays <- c("Raw" = "counts", "Normalized" = "exprs")
  }
  features <-
    c("all", as.character(unique(rowData(sce)$marker_class)))
  
  ## COUNTS
  
  # ui for counts
  output$designCounts <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
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
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),

      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("countsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("countsPlotDownload"),
      style = "position: absolute; bottom: 10px; right:10px;"
    ))
  })
  
  # render counts plot
  output$countsPlot <- renderPlot({
    reactiveVals$countsPlot <-  CATALYST::plotCounts(
      sce,
      group_by = input$countsGroupBy,
      color_by = input$countsColorBy,
      prop = as.logical(input$countsProp)
    )
    reactiveVals$countsPlot
  })
  
  # ui for download button
  output$countsPlotDownload <- renderUI({
    req(reactiveVals$countsPlot)
    downloadButton("downloadPlotCounts", "Download Plot")
  })
  
  # function for downloading count plot
  output$downloadPlotCounts <- downloadHandler(
    filename = function(){
      paste0("Counts_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$countsPlot, width=12, height=6)
    }
  )
  
  ## MDS 
  
  # ui for MDS
  output$designMDS <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
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
        selectizeInput("mdsFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      ),
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("mdsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("mdsPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px"
    ),)
  })
  
  # render mds plot
  output$mdsPlot <- renderPlot({
    feature <- input$mdsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$mdsPlot <- CATALYST::pbMDS(
      sce,
      label_by = input$mdsLabelBy,
      color_by = input$mdsColorBy,
      features = feature,
      assay = input$mdsAssay,
    )
    reactiveVals$mdsPlot
    
  })
  
  # ui for download button
  output$mdsPlotDownload <- renderUI({
    req(reactiveVals$mdsPlot)
    downloadButton("downloadPlotMDS", "Download Plot")
  })
  
  # function for downloading MDS plot
  output$downloadPlotMDS <- downloadHandler(
    filename = function(){
      paste0("MDS_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$mdsPlot, width=16, height=11)
    }
  )
  
  ## NRS
  
  # ui for NRS
  output$designNRS <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
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
        selectizeInput("nrsFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("nrsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("nrsPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px;"
    ))
  })
  
  # render nrs plot
  output$nrsPlot <- renderPlot({
    feature <- input$nrsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$nrsPlot <- CATALYST::plotNRS(
      sce,
      color_by = input$nrsColorBy,
      features = feature,
      assay = input$nrsAssay
    )
    reactiveVals$nrsPlot
  })
  
  # ui for download button
  output$nrsPlotDownload <- renderUI({
    req(reactiveVals$nrsPlot)
    downloadButton("downloadPlotNRS", "Download Plot")
  })
  
  # function for downloading NRS plot
  output$downloadPlotNRS <- downloadHandler(
    filename = function(){
      paste0("NRS_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$nrsPlot, width=12, height=6)
    }
  )
  
  ## Exprs

  # ui for expr
  output$designExprs <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
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
        selectizeInput("exprsFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("exprsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("exprsPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px;"
    ),)
  })
  
  # render exprs plot
  output$exprsPlot <- renderPlot({
    feature <- input$exprsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$exprsPlot <- CATALYST::plotExprs(
      sce,
      color_by = input$exprsColorBy,
      features = feature,
      assay = input$exprsAssay
    )
    reactiveVals$exprsPlot
  })
  
  # ui for download button
  output$exprsPlotDownload <- renderUI({
    req(reactiveVals$exprsPlot)
    downloadButton("downloadPlotExprs", "Download Plot")
  })
  
  # function for downloading exprs plot
  output$downloadPlotExprs <- downloadHandler(
    filename = function(){
      paste0("Expr_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$exprsPlot, width=14, height=9)
    }
  )
  

  ## Exprs Heatmap
  
  # ui for exprs heatmap
  output$designExprsHeatmap <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
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
        selectizeInput("exprsHeatmapFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("exprsHeatmapPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("exprsHeatmapPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px;"
    ),)
  })
  
  # render exprs heatmap plot
  output$exprsHeatmapPlot <- renderPlot({
    feature <- input$exprsHeatmapFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$exprsPlotHeatmap <- CATALYST::plotExprHeatmap(
      sce,
      scale = input$exprsHeatmapScale,
      features = feature,
      assay = input$exprsHeatmapAssay
    )
    reactiveVals$exprsPlotHeatmap
  })
  
  # ui for download button
  output$exprsHeatmapPlotDownload <- renderUI({
    req(reactiveVals$exprsPlotHeatmap)
    library(ComplexHeatmap)
    downloadButton("downloadPlotExprsHeatmap", "Download Plot")
  })
  
  # function for downloading exprs heatmap
  output$downloadPlotExprsHeatmap <- downloadHandler(
    filename = "Expression_Heatmap.pdf", 
    content = function(file){
      pdf(file, width = 12, height = 8)
      draw(reactiveVals$exprsPlotHeatmap)
      dev.off()
    }
  )
  
}