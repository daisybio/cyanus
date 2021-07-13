# Preprocessing Server

resetPreprocessing <- function(){
  reactiveVals$preprocessingShowed <- FALSE
}

# marker, sample, patient selection -> if markers, patients, samples are selected -> prepValButton can be clicked
observeEvent({
  input$patientSelection
  input$sampleSelection
}, {
  if (length(input$sampleSelection) == 0 || length(input$patientSelection)==0){
    shinyjs::disable("prepSelectionButton")
    shinyjs::disable("filterSelectionButton")
  } else {
    shinyjs::enable("prepSelectionButton")
    shinyjs::enable("filterSelectionButton")
  }
}, ignoreNULL = FALSE)

# check current tab
observe({
  if (reactiveVals$current_tab == 3) {
    if (!reactiveVals$preprocessingShowed){
      plotPreprocessing(reactiveVals$sce)
      reactiveVals$preprocessingShowed <- TRUE
    }
    if (!("patient_id" %in% colnames(colData(reactiveVals$sce)))) {
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

# render sortable lists for all conditions
output$reorderingTabs <- renderUI({
  library(sortable)
  conditions <- names(metadata(reactiveVals$sce)$experiment_info)
  conditions <- conditions[!conditions %in% c("sample_id", "patient_id", "n_cells")]
  lapply(conditions, function(condition){
    print(condition)
    rank_list(
      text = paste0("Reorder the condition: ",condition),
      labels = levels(metadata(reactiveVals$sce)$experiment_info[[condition]]),
      input_id = condition
    )
  })
  
})

# if conditions are ordered
observeEvent(input$reorderButton, {
  waiter_show(id = "app",html = tagList(spinner$logo, 
                                        HTML("<br>Reordering Data...<br>Please be patient")), 
              color=spinner$color)
  
  # reorder levels 
  conditions <- names(metadata(reactiveVals$sce)$experiment_info)
  conditions <- conditions[!conditions %in% c("sample_id", "patient_id", "n_cells")]
  lapply(conditions, function(condition){
    ordered <- input[[condition]]
    reactiveVals$sce[[condition]] <- factor(reactiveVals$sce[[condition]], levels=ordered)
    metadata(reactiveVals$sce)$experiment_info[[condition]] <- factor(metadata(reactiveVals$sce)$experiment_info[[condition]], levels=ordered)
  })
  plotPreprocessing(reactiveVals$sce)
  waiter_hide(id = "app")

})

# if start transformation button is clicked
observeEvent(input$prepButton, {
  waiter_show(id = "app",html = tagList(spinner$logo, 
                                        HTML("<br>Transforming Data...<br>Please be patient")), 
              color=spinner$color)
  # data transformation
  reactiveVals$sce <-
    transformData(sce = reactiveVals$sce,
                  cf = as.numeric(input$cofactor))
  plotPreprocessing(reactiveVals$sce)
  waiter_hide(id = "app")
  runjs("document.getElementById('nextTab').scrollIntoView();")
})

# if visualize selection button is clicked
observeEvent(input$prepSelectionButton, {
  waiter_show(id = "app",html = tagList(spinner$logo, 
                                        HTML("<br>Applying Selection...<br>Please be patient")), 
              color=spinner$color)
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
  waiter_hide(id = "app")
})

# if filtering button is clicked -> selection is applied to sce
observeEvent(input$filterSelectionButton,{
  waiter_show(id = "app",html = tagList(spinner$logo, 
                                        HTML("<br>Applying Filtering...<br>Please be patient")), 
              color=spinner$color)
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
  
  waiter_hide(id = "app")
})

observeEvent(reactiveVals$sce, {
  if ("exprs" %in% names(assays(reactiveVals$sce)))
    shinyjs::hide("noTransformationWarning")
  else 
    shinyjs::show("noTransformationWarning")
})

# method for plotting all kinds of preprocessing plots
plotPreprocessing <- function(sce) {
  groupColorLabelBy <- names(colData(sce))
  possAssays <- assayNames(sce)
  if (all(possAssays == c("counts", "exprs"))) {
    possAssays <- c("Normalized" = "exprs", "Raw" = "counts")
  }
  if (all(possAssays == c("counts"))){
    possAssays <- c("Raw" = "counts")
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
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      ggsave(file, plot = reactiveVals$countsPlot, width=12, height=6)
      waiter_hide(id="app")
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
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      ggsave(file, plot = reactiveVals$mdsPlot, width=16, height=11)
      waiter_hide(id="app")
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
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      ggsave(file, plot = reactiveVals$nrsPlot, width=12, height=6)
      waiter_hide(id="app")
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
    column(11, div(HTML("Plot of the smoothed expression densities. <b style='color:#FF3358';>Attention:</b> This may take some time")),
      shinycssloaders::withSpinner(
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
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      ggsave(file, plot = reactiveVals$exprsPlot, width=14, height=9)
      waiter_hide(id="app")
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
    reactiveVals$exprsPlotHeatmap <- plotExprHeatmapCustom(
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
    downloadButton("downloadPlotExprsHeatmap", "Download Plot")
  })
  
  # function for downloading exprs heatmap
  output$downloadPlotExprsHeatmap <- downloadHandler(
    filename = "Expression_Heatmap.pdf", 
    content = function(file){
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      pdf(file, width = 12, height = 8)
      ComplexHeatmap::draw(reactiveVals$exprsPlotHeatmap)
      dev.off()
      waiter_hide(id="app")
    }
  )
  
}