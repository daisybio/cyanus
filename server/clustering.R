### Helpers ----
source("server/clusterFun.R", local = TRUE)

### Observer ----
observeEvent(input$startClustering, {
  updateButton(session,
               "startClustering",
               label = " Clustering...",
               disabled = TRUE)
  toggle_inputs()
  
  showNotification(ui =
                     HTML(
                       sprintf(
                         "<div id='clusteringProgress'><b>Clustering Progress with %s:</b><div>",
                         input$clusteringMethod
                       )
                     ),
                   duration = NULL,
                   id = "clusteringProgressNote")
  
  withCallingHandlers({
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
  },
  message = function(m) {
    shinyjs::html(id = "clusteringProgress",
                  html = sprintf("<br>%s", HTML(m$message)),
                  add = TRUE)
  })
  
  toggle_inputs(enable_inputs = TRUE)
  updateButton(session,
               "startClustering",
               label = " Start Clustering",
               disabled = FALSE)
  updateButton(session, "continue", label = " Differential Expression Analysis")
  shinyjs::show("continue")
  removeNotification("clusteringProgressNote")
  showNotification(HTML(
    sprintf(
      "<b>Finished clustering with %s.</b><br>
      You can visualize the newly assigned clusters in the Visualization tab.",
      input$clusteringMethod
    )
  ),
  duration = 10,
  type = "message")
})

observeEvent(input$visualizeClustering, {
  method <- isolate(input$clusteringRuns)
  sce <- isolate(reactiveVals$sce)
  clusterCode <- isolate(input$clusterCode)
  assay <- isolate(input$assayTypeVisIn)
  by <- isolate(input$abundanceBy)
  group_by <- isolate(input$abundanceGroup)
  shape_by <- isolate(input$abundanceShape)
  if (is.null(shape_by) | shape_by == "")
    shape_by <- NULL
  if (is.null(clusterCode) | method != "flowSOM")
    clusterCode <- NULL
  
  reactiveVals[[method]]$clusterAbundancePlot <-
    plotAbundancesCustom(
      sce,
      method = method,
      k = clusterCode,
      by = by,
      group_by = group_by,
      shape_by = shape_by
    )
  # TODO: implement for rphenograph
  
  reactiveVals[[method]]$clusterHeatmapPlot <-
    plotFreqHeatmapCustom(isolate(reactiveVals$sce),
                          method,
                          clusterCode)
  
  reactiveVals[[method]]$clusterExprsPlot <-
    plotClusterExprsCustom(
      sce,
      method = method,
      k = clusterCode,
      features = metadata(sce)$cluster_run[[method]]$features,
      assay
    )
  
  
  
})

observeEvent(input$featuresIn, {
  if (input$useFeaturesIn == "Marker by Class")
    reactiveVals$featureNames <-
      rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) %in% input$featuresIn]
  else
    reactiveVals$featureNames <- input$featuresIn
})

observe({
  if (length(reactiveVals$featureNames) < 2)
    disable("startClustering")
  else
    enable("startClustering")
})

observeEvent(input$mergeClusteringButton, {
  req(input$mergeClusters_cell_edit)
  reactiveVals$mergingFrame <-
    editData(isolate(reactiveVals$mergingFrame),
             input$mergeClusters_cell_edit)
  reactiveVals$sce <-
    mergeClusters(
      isolate(reactiveVals$sce),
      isolate(input$clusterCode),
      isolate(reactiveVals$mergingFrame),
      id = sprintf("merged_%s", isolate(input$clusterCode)),
      overwrite = TRUE
    )
  runjs("('#mergeClusteringButton').closest('.box-header').find('[data-widget=collapse]').click();")
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
           uiOutput("ydim"),
           width = 6),
    title = "Choose Clustering Parameters",
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE
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
    selected <- "type"
  } else if (input$useFeaturesIn == "Marker by Name") {
    choices <- rownames(reactiveVals$sce)
    names(choices) <-
      sprintf("%s (%s)", choices, as.character(marker_classes(reactiveVals$sce)))
    selected <-
      rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) == "type"]
  } else
    stop("by name or by class?")
  shinyWidgets::pickerInput(
    inputId = "featuresIn",
    label = "Features to use for clustering",
    choices = choices,
    selected = selected,
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
  selectizeInput("assayTypeIn",
                 "Expression Type",
                 choices = choices)
})

output$k <- renderUI({
  req(input$clusteringMethod != "clusterX",
      length(reactiveVals$featureNames) > 1)
  
  if (input$clusteringMethod == "flowSOM") {
    label = "Maximum Number of Clusters to Evaluate in the Metaclustering"
    value = 20
    maxK = 100
  } else if (input$clusteringMethod == "rphenoGraph") {
    label = "Number of Nearest Neighbours"
    value = 30
    maxK = length(reactiveVals$featureNames) - 1
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
  
  abundanceByChoices <-
    c("Samples" = "sample_id", "Clusters" = "cluster_id")
  abundanceChoices <- names(colData(reactiveVals$sce))
  names(abundanceChoices) <- names(abundanceChoices)
  
  shinydashboard::box(
    fluidRow(
      column(
        selectizeInput(
          "clusteringRuns",
          "Successfull Run",
          names(metadata(reactiveVals$sce)$cluster_run)
        ),
        width = 6
      ),
      column(uiOutput("assayTypeVisOut"),
             width = 6)
    ),
    column(tableOutput("clusterRunParams"),
           width = 8),
    column(uiOutput("selectClusterCode"),
           width = 4),
    column(withSpinner(uiOutput("clusterSizes")),
           width = 12,
           style = "overflow-x: scroll;"),
    fluidRow(withSpinner(uiOutput("delta_area"))),
    fluidRow(withSpinner(uiOutput("clusterMergingBox")),
    ),
    fluidRow(withSpinner(uiOutput(
      "clusteringOutput"
    ))),
    fluidRow(
      column(
        width = 4,
        selectizeInput("abundanceBy", "By Samples or Clusters?", abundanceByChoices)
      ),
      column(
        width = 4,
        selectizeInput("abundanceGroup", "Group By", abundanceChoices[!abundanceChoices %in% abundanceByChoices])
      ),
      column(width = 4,
             selectizeInput(
               "abundanceShape",
               "Shape by",
               c("Nothing" = "", abundanceChoices[!abundanceChoices %in% abundanceByChoices]),
             ))
    ),
    div(
      downloadButton("downloadClusters", "Download Cluster Assignments"),
      style = "float: left;"
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
    title = "Visualize Clustering Results",
    width = 12
  )
})

output$clusterRunParams <- renderTable({
  runjs(
    "document.getElementById('clusteringVisualizationSelection').scrollIntoView();"
  )
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

output$assayTypeVisOut <- renderUI({
  choices <- c("Transformed" = "exprs", "Raw" = "counts")
  choices <- choices[choices %in% assayNames(reactiveVals$sce)]
  selectizeInput("assayTypeVisIn",
                 "Expression Type",
                 choices = choices)
})

output$clusterSizes <- renderTable({
  res <-
    as.data.frame(table(
      cluster_ids(reactiveVals$sce, input$clusterCode, input$clusteringRuns),
      useNA =
        "ifany"
    ))
  names(res) <- c("Cluster", "Size")
  t(res)
},
caption = "Cluster Sizes",
caption.placement = "top", rownames = TRUE, colnames = FALSE)

output$mergeClusters <- renderDT({
  req(input$clusteringRuns)
  curr_cluster_ids <-
    levels(
      cluster_ids(
        x = reactiveVals$sce,
        k = input$clusterCode,
        method = input$clusteringRuns
      )
    )
  reactiveVals$mergingFrame <-
    data.frame(
      old_cluster = curr_cluster_ids,
      new_cluster = rep("new_cluster_id", times = length(curr_cluster_ids))
    )
  reactiveVals$mergingFrame
},
editable = list(target = "cell", disable = list(columns = 1)),
options = list(pageLength = nlevels(
  cluster_ids(
    x = reactiveVals$sce,
    k = input$clusterCode,
    method = input$clusteringRuns
  )
)))

output$delta_area <- renderUI({
  req(input$clusteringRuns == "flowSOM")
  
  
  runjs(
    "document.getElementById('clusteringVisualizationSelection').scrollIntoView();"
  )
  shinydashboard::box(
    renderPlotly(ggplotly(
      CATALYST::delta_area(reactiveVals$sce)
    )),
    title = "Delta Area",
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE
  )
})

output$clusterMergingBox <- renderUI({
  req(input$clusteringRuns == "flowSOM")
  
  shinydashboard::box(
    HTML(
      "You can assign new cluster names in the <b>new_cluster</b> column by double-clicking.<br>
                                      <i>Make sure to assign new names to <b>all</b> clusters.</i>"
    ),
    DTOutput("mergeClusters"),
    div(
      bsButton(
        "mergeClusteringButton",
        "Merge Clustering",
        icon = icon("border-all"),
        style = "success"
      ),
      style = "margin-top: 5px; float: right;"
    ),
    title = "Merge Clusters",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
  )})

output$downloadClusters <- downloadHandler(
  filename = function() {
    paste(input$clusteringRuns, ".csv", sep = "")
  },
  content = function(file) {
    to_write <- colData(reactiveVals$sce)
    if (input$clusteringRuns == "flowSOM") {
      add_meta <- list(cluster_ids(reactiveVals$sce, input$clusterCode))
      names(add_meta) <- input$clusterCode
      to_write <- cbind(to_write, add_meta)
    }
    write.csv(to_write, file, row.names = FALSE)
  }
)

output$clusteringOutput <- renderUI({
  req(input$visualizeClustering, reactiveVals[[input$clusteringRuns]])
  
  shinydashboard::box(
    fluidRow(
      shinydashboard::tabBox(
        tabPanel(
          "Cluster Abundances",
          div("Relative population abundances of the specified clustering."),
          fluidRow(
            withSpinner(plotOutput("clusterAbundancePlot",
                                   height = "800px")),
            div(uiOutput("clusterAbundanceDownload"),
                style = "float: right;")
          )
        ),
        tabPanel(
          "Cluster Frequencies",
          div(
            "Heatmap of relative cluster abundances (frequencies) by sample."
          ),
          fluidRow(
            withSpinner(plotOutput("clusterHeatmapPlot",
                                   height = "800px")),
            div(uiOutput("clusterHeatmapDownload"),
                style = "float: right;")
          )
        ),
        tabPanel(
          "Marker Densities",
          div("Smoothed densities of marker intensities by cluster."),
          fluidRow(
            withSpinner(plotOutput("clusterExprsPlot",
                                   height = "800px")),
            div(uiOutput("clusterDensitiyDownload"),
                style = "float: right;")
          )
        ),
        title = "Cluster Visualization",
        width = 12
      )
    ),
    title = "Cluster Visualizations",
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE
  )
})

output$clusterAbundanceDownload <- renderUI({
  req(reactiveVals$abundanceCluster)
  
  downloadButton("downloadPlotAbundance", "Download Plot")
})

output$clusterDensitiyDownload <- renderUI({
  req(reactiveVals$exprsCluster)
  
  downloadButton("downloadPlotDensity", "Download Plot")
})

output$clusterHeatmapDownload <- renderUI({
  req(reactiveVals$heatmapCluster)
  
  downloadButton("downloadPlotFrequency", "Download Plot")
})

output$clusterAbundancePlot <- renderPlot({
  reactiveVals$abundanceCluster <-
    reactiveVals[[input$clusteringRuns]]$clusterAbundancePlot
  reactiveVals$abundanceCluster
})

output$clusterHeatmapPlot <- renderPlot({
  reactiveVals$heatmapCluster <-
    reactiveVals[[input$clusteringRuns]]$clusterHeatmapPlot
  reactiveVals$heatmapCluster
})

output$clusterExprsPlot <- renderPlot({
  reactiveVals$exprsCluster <-
    reactiveVals[[input$clusteringRuns]]$clusterExprsPlot
  reactiveVals$exprsCluster
})

output$downloadPlotAbundance <-
  downloadPlotFunction("Population_Abundances", reactiveVals$abundanceCluster)

output$downloadPlotDensity <-
  downloadPlotFunction(
    "Cluster_Expression",
    reactiveVals$exprsCluster,
    width = 16,
    height = 12
  )

output$downloadPlotFrequency <- downloadHandler(
  filename = "Cluster_Heatmap.pdf",
  content = function(file) {
    pdf(file, width = 12, height = 8)
    draw(reactiveVals$heatmapCluster)
    dev.off()
  }
)
