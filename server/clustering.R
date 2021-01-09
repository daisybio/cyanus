### Helpers ----
source("server/clusterFun.R", local = TRUE)

### Observer ----
observeEvent(input$startClustering, {
  updateButton(session,
               "startClustering",
               label = " Clustering...",
               disabled = TRUE)
  toggle_inputs()
  
  showNotification(
    ui =
      HTML(
        "<div id='clusteringProgress'><b>Clustering Progress:</b><div>"
      ),
    duration = NULL,
    id = "clusteringProgressNote"
  )
  
  withCallingHandlers({
    reactiveVals$sce <-
      clusterSCE(
        reactiveVals$sce,
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
  
  reactiveVals$sce <-
    mergeClusters(
      reactiveVals$sce,
      k = "meta20",
      id = "all",
      table = data.frame(old_cluster = seq_len(20), new_cluster = "all")
    )
  
  assays <- c("exprs" = "Transformed", "counts" = "Raw")
  
  reactiveVals$clusterRun <- list(
    features = CATALYST:::.get_features(reactiveVals$sce, input$featuresIn),
    assayType = assays[input$assayTypeIn],
    xdim = input$xdim,
    ydim = input$ydim,
    maxK = input$k
  )
  
  toggle_inputs(enable_inputs = TRUE)
  updateButton(session,
               "startClustering",
               label = " Start Clustering",
               disabled = FALSE)
  updateButton(session, "continue", label = " Differential Expression Analysis")
  shinyjs::show("continue")
  removeNotification("clusteringProgressNote")
  showNotification(
    HTML(
      "<b>Finished clustering.</b><br>
      You can visualize the newly assigned clusters in the Visualization tab."
    ),
    duration = 10,
    type = "message"
  )
})

observeEvent(input$clusterCode, {
  curr_cluster_ids <-
    levels(cluster_ids(x = reactiveVals$sce,
                       k = input$clusterCode))
  reactiveVals$mergingFrame <-
    data.frame(
      old_cluster = curr_cluster_ids,
      new_cluster = rep("new_cluster_id", times = length(curr_cluster_ids))
    )
})

observeEvent(input$mergeClustersDT_cell_edit, {
  edit <- input$mergeClustersDT_cell_edit
  edit$col <- edit$col + 1
  reactiveVals$mergingFrame <-
    editData(reactiveVals$mergingFrame,
             edit)
  
})

observeEvent(input$mergeClusteringButton, {
  reactiveVals$sce <-
    mergeClusters(
      isolate(reactiveVals$sce),
      isolate(input$clusterCode),
      isolate(reactiveVals$mergingFrame),
      id = sprintf("merged_%s", isolate(input$clusterCode)),
      overwrite = TRUE
    )
  
  # runjs(
  #   "$('#mergeClusteringButton').closest('.box-header').find('[data-widget=collapse]').click();"
  # )
})

### Renderer ----

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

output$clusteringVisualizationSelection <- renderUI({
  req(reactiveVals$clusterRun)
  
  shinydashboard::box(
    column(tableOutput("clusterRunParams"),
           width = 8, 
           style = "overflow-x: scroll;"),
    column(
      uiOutput("selectClusterCode"),
      div(
        downloadButton("downloadClusters", "Download Cluster Assignments"),
        style = "float: right;"
      ),
      width = 4
    ),
    column(withSpinner(uiOutput("clusterSizes")),
           width = 12,
           style = "overflow-x: scroll;"),
    fluidRow(withSpinner(uiOutput("delta_area"))),
    fluidRow(withSpinner(uiOutput(
      "clusteringOutput"
    ))),
    fluidRow(withSpinner(uiOutput(
      "clusterMergingBox"
    ))),
    title = "Clustering Output",
    width = 12
  )
})

output$clusterRunParams <- renderTable({
  runjs(
    "document.getElementById('clusteringVisualizationSelection').scrollIntoView();"
  )
  runParams <-
    reactiveVals$clusterRun
  runParams$features <- paste(runParams$features, collapse = ",")
  data.frame(runParams)
},
caption = "Run Parameters",
caption.placement = "top")

output$selectClusterCode <- renderUI({
  req(reactiveVals$clusterRun)
  
  choicesClusterCode <- names(cluster_codes(reactiveVals$sce))
  choicesClusterCode <-
    c('all', choicesClusterCode[choicesClusterCode != 'all'])
  selectInput("clusterCode",
              "Clusters",
              rev(choicesClusterCode))
})

output$clusterSizes <- renderTable({
  req(input$clusterCode)
  res <-
    as.data.frame(table(cluster_ids(reactiveVals$sce, input$clusterCode),
                        useNA =
                          "ifany"))
  names(res) <- c("Cluster", "Size")
  t(res)
},
caption = "Cluster Sizes",
caption.placement = "top", rownames = TRUE, colnames = FALSE)

output$mergeClustersDT <- renderDT({
  req(reactiveVals$clusterRun)
  
  reactiveVals$mergingFrame
},
rownames = FALSE,
editable = list(target = "cell", disable = list(columns = 0)),
selection = "none",
options = list(pageLength = nlevels(
  cluster_ids(x = reactiveVals$sce,
              k = input$clusterCode)
)))

output$delta_area <- renderUI({
  req(reactiveVals$clusterRun)
  
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
    collapsed = TRUE
  )
})

output$clusterMergingBox <- renderUI({
  req(reactiveVals$clusterRun)
  
  shinydashboard::box(
    HTML(
      "You can assign new cluster names in the <b>new_cluster</b> column by double-clicking.<br>
                                      <i>Make sure to assign new names to <b>all</b> clusters.</i>"
    ),
    withSpinner(DTOutput("mergeClustersDT")),
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
  )
})

output$downloadClusters <- downloadHandler(
  filename = "Clustering_Assignments.csv",
  content = function(file) {
    add_meta <- list(cluster_ids(reactiveVals$sce, input$clusterCode))
    names(add_meta) <- input$clusterCode
    
    write.csv(cbind(colData(reactiveVals$sce), add_meta), file, row.names = FALSE)
  }
)

output$clusteringOutput <- renderUI({
  req(reactiveVals$clusterRun)
  
  abundanceByChoices <-
    c("Samples" = "sample_id", "Clusters" = "cluster_id")
  abundanceChoices <- names(colData(reactiveVals$sce))
  names(abundanceChoices) <- names(abundanceChoices)
  
  densityChoices <- c("Transformed" = "exprs", "Raw" = "counts")
  densityChoices <-
    densityChoices[densityChoices %in% assayNames(reactiveVals$sce)]
  
  shinydashboard::box(
    fluidRow(
      shinydashboard::tabBox(
        tabPanel(
          "Cluster Abundances",
          div("Relative population abundances of the specified clustering."),
          br(),
          div(
            dropdownButton(
              tags$h3("Plot Options"),
              selectizeInput("abundanceBy", "By Samples or Clusters?", abundanceByChoices),
              selectizeInput("abundanceGroup", "Group By", abundanceChoices[!abundanceChoices %in% abundanceByChoices]),
              conditionalPanel(
                "input.abundanceBy == 'cluster_id'",
                selectizeInput(
                  "abundanceShape",
                  "Shape by",
                  c("Nothing" = "", abundanceChoices[!abundanceChoices %in% abundanceByChoices])
                )
              ),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div(uiOutput("clusterAbundanceDownload"),
              style = "float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterAbundancePlot",
                       height = "800px")
          ))
        ),
        tabPanel(
          "Cluster Frequencies",
          div(
            "Heatmap of relative cluster abundances (frequencies) by sample."
          ),
          div(uiOutput("clusterHeatmapDownload"),
              style = "float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterHeatmapPlot",
                       height = "800px")
          ))
        ),
        tabPanel(
          "Marker Densities",
          div("Smoothed densities of marker intensities by cluster."),
          br(),
          div(
            dropdownButton(
              tags$h3("Plot Options"),
              selectizeInput("assayTypeVisIn",
                             "Expression Type",
                             choices = densityChoices),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div(uiOutput("clusterDensitiyDownload"),
              style = "float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterExprsPlot",
                       height = "800px")
          ))
        ),
        title = "Cluster Visualization",
        id = "clusterVisTabBox",
        width = 12
      )
    ),
    title = "Visualize Clustering Results",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
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
  shape_by <- input$abundanceShape
  if (is.null(shape_by) | shape_by == "")
    shape_by <- NULL
  reactiveVals$abundanceCluster <-
    plotAbundancesCustom(
      reactiveVals$sce,
      k = input$clusterCode,
      by = input$abundanceBy,
      group_by = input$abundanceGroup,
      shape_by = shape_by
    )
  reactiveVals$abundanceCluster
})

output$clusterHeatmapPlot <- renderPlot({
  reactiveVals$heatmapCluster <-
    plotFreqHeatmapCustom(reactiveVals$sce,
                          input$clusterCode)
  reactiveVals$heatmapCluster
})

output$clusterExprsPlot <- renderPlot({
  reactiveVals$exprsCluster <-
    plotClusterExprsCustom(
      reactiveVals$sce,
      k = input$clusterCode,
      features = reactiveVals$clusterRun$features,
      input$assayTypeVisIn
    )
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
