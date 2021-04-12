### Helpers ----
source("server/clusterFun.R", local = TRUE)

### Observer ----
observeEvent(input$startClustering, {
  updateButton(session,
               "startClustering",
               label = " Clustering...",
               disabled = TRUE)
  
  waiter_show(html = tagList(spinner$logo, 
                             HTML("<br>Clustering in Progress...<br>Please be patient")), 
              color=spinner$color)
  
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
        reactiveVals$clusterFeatureNames,
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
  
  assays <- c("exprs" = "Transformed", "counts" = "Raw")
  
  reactiveVals$clusterRun <- list(
    features = reactiveVals$clusterFeatureNames,
    assayType = assays[input$assayTypeIn],
    xdim = input$xdim,
    ydim = input$ydim,
    maxK = input$k
  )
  
  waiter_hide()
  
  updateButton(session,
               "startClustering",
               label = " Start Clustering",
               disabled = FALSE)
  reactiveVals$continue <- TRUE
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

observeEvent(input$featuresIn,
             {
               if (input$useFeaturesIn == "Marker by Class")
                 reactiveVals$clusterFeatureNames <-
                   rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) %in% input$featuresIn]
               else
                 reactiveVals$clusterFeatureNames <-
                   input$featuresIn
             },
             ignoreNULL = FALSE,
             ignoreInit = TRUE)

observe({
  if (length(reactiveVals$clusterFeatureNames) == 0)
    disable("startClustering")
  else
    enable("startClustering")
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

observeEvent(input$clusterCode, {
  if (nlevels(cluster_ids(reactiveVals$sce, input$clusterCode)) == 1) {
    hideTab("clusterVisTabBox", "Cluster Frequencies", session = getDefaultReactiveDomain())
    hideTab("clusterVisTabBox", "Marker Densities", session = getDefaultReactiveDomain())
  } else {
    showTab(
      "clusterVisTabBox",
      "Cluster Frequencies",
      select = FALSE,
      session = getDefaultReactiveDomain()
    )
    showTab(
      "clusterVisTabBox",
      "Marker Densities",
      select = FALSE,
      session = getDefaultReactiveDomain()
    )
  }
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
    choices <-
      sortMarkerNames(choices, as.character(marker_classes(reactiveVals$sce)), first = "type")
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
    column(
      tableOutput("clusterRunParams"),
      width = 8,
      style = "overflow-x: scroll;"
    ),
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
              "Meta-Clusters",
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
    div(
      'It is recommended to choose a metacluster where the plateau is reached.
                          "The delta area represents the amount of extra cluster stability gained when clustering into k groups as compared to k-1 groups.
                          It can be expected that high stability of clusters can be reached when clustering into the number of groups that best fits the data.
                          The "natural" number of clusters present in the data should thus corresponds to the value of k where there is no longer a considerable increase in stability (plateau onset)."',
      style = "text-align: center;vertical-align: middle;"
    ),
    renderPlotly(ggplotly(
      CATALYST::delta_area(reactiveVals$sce)
    )),
    title = "1. Delta Area",
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE
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
    title = "3. Merge Clusters",
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
  
  densityAssayChoices <-
    c("Transformed" = "exprs", "Raw" = "counts")
  densityAssayChoices <-
    densityAssayChoices[densityAssayChoices %in% assayNames(reactiveVals$sce)]
  
  # sceEI <- ei(reactiveVals$sce)
  # starMarkerFacets <- names(which(sapply(sceEI, function(feature) nlevels(as.factor(feature)) == 2)))
  # names(starMarkerFacets) <- starMarkerFacets
  # starMarkerFacets <- c("No Facetting" = NA, starMarkerFacets)
  
  shinydashboard::box(
    fluidRow(
      shinydashboard::tabBox(
        tabPanel(
          "Cluster Abundances",
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
          div(
            "Relative population abundances of the specified clustering.",
            style = "text-align: center;vertical-align: middle;"
          ),
          div(uiOutput("clusterAbundanceDownload"),
              style = "position: relative; z-index: 99; float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterAbundancePlot",
                       height = "800px")
          ))
        ),
        tabPanel(
          "Cluster Frequencies",
          div(
            "Heatmap of relative cluster abundances (frequencies) by sample."
            ,
            style = "text-align: center;vertical-align: middle;"
          ),
          div(uiOutput("clusterHeatmapDownload"),
              style = "position: relative; z-index: 99; float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterHeatmapPlot", height = "800px")
          ))
        ),
        tabPanel(
          "Marker Densities",
          div(
            dropdownButton(
              tags$h3("Plot Options"),
              selectizeInput("assayTypeVisIn",
                             "Expression Type",
                             choices = densityAssayChoices),
              selectInput(
                "densityUseFeatureChoicesIn",
                label = "Features",
                choices =
                  c(
                    "all",
                    levels(
                      SummarizedExperiment::rowData(reactiveVals$sce)$marker_class
                    )
                  ),
                selected = "type"
              ),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div("Smoothed densities of marker intensities by cluster.", style = "text-align: center;vertical-align: middle;"),
          div(uiOutput("clusterDensitiyDownload"),
              style = "position: relative; z-index: 99; float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterExprsPlot",
                       height = "800px")
          ))
        ),
        tabPanel(
          "Star Chart (Overall)",
          div(
            "Tree, where each node (cluster) is represented by a star chart indicating median marker values.",
            style = "text-align: center;vertical-align: middle;"
          ),
          div(uiOutput("clusterStarDownload"),
              style = "position: relative; z-index: 99; float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterStarPlot",
                       height = "800px")
          ))
        ),
        tabPanel(
          "Star Chart (Markerwise)",
          div(
            dropdownButton(
              tags$h3("Plot Options"),
              selectInput(
                "plotStarMarkerFeatureIn",
                label = "Marker",
                choices = stats::setNames(
                  rownames(reactiveVals$sce),
                  sprintf(
                    "%s (%s)",
                    rownames(reactiveVals$sce),
                    as.character(marker_classes(reactiveVals$sce))
                  )
                )
              ),
              # selectInput(
              #   "plotStarMarkerFacets",
              #   label = "Facet By",
              #   choices = starMarkerFacets
              # ),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div(
            "Tree, where each node (cluster) is coloured depending on its median value for the given marker.",
            style = "text-align: center;vertical-align: middle;"
          ),
          div(uiOutput("clusterStarMarkerDownload"),
              style = "position: relative; z-index: 99; float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterStarMarkerPlot",
                       height = "800px")
          ))
        ),
        title = "Cluster Visualization",
        id = "clusterVisTabBox",
        width = 12
      )
    ),
    title = "2. Visualize Clustering",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
  )
})

output$clusterStarMarkerDownload <- renderUI({
  req(reactiveVals$starMarkerCluster)
  
  downloadButton("downloadPlotMarkerStar", "Download Plot")
})

output$clusterStarDownload <- renderUI({
  req(reactiveVals$starCluster)
  
  downloadButton("downloadPlotStar", "Download Plot")
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
  library(ComplexHeatmap)
  downloadButton("downloadPlotFrequency", "Download Plot")
})

output$clusterStarPlot <- renderPlot({
  reactiveVals$starCluster <-
    plotStarsCustom(reactiveVals$sce, backgroundValues = cluster_codes(reactiveVals$sce)[[input$clusterCode]])
  reactiveVals$starCluster
})

output$clusterStarMarkerPlot <- renderPlot({
  reactiveVals$starMarkerCluster <-
    plotMarkerCustom(
      reactiveVals$sce,
      input$plotStarMarkerFeatureIn,
      # facet_by = input$plotStarMarkerFacets,
      backgroundValues = cluster_codes(reactiveVals$sce)[[input$clusterCode]]
    )
  reactiveVals$starMarkerCluster
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
  req(nlevels(cluster_ids(reactiveVals$sce, input$clusterCode)) > 1)
  reactiveVals$heatmapCluster <-
    plotFreqHeatmapCustom(reactiveVals$sce,
                          input$clusterCode)
  reactiveVals$heatmapCluster
})

output$clusterExprsPlot <- renderPlot({
  req(nlevels(cluster_ids(reactiveVals$sce, input$clusterCode)) > 1)
  features_tmp <- input$densityUseFeatureChoicesIn
  if (features_tmp == "all")
    features_tmp <- NULL
  reactiveVals$exprsCluster <-
    plotClusterExprsCustom(
      reactiveVals$sce,
      k = input$clusterCode,
      features = features_tmp,
      input$assayTypeVisIn
    )
  reactiveVals$exprsCluster
})

output$downloadPlotStar <- downloadHandler(
  filename = "Star_Charts_overall.pdf",
  content = function(file) {
    pdf(file, width = 12, height = 8)
    print(reactiveVals$starCluster)
    dev.off()
  }
)

output$downloadPlotMarkerStar <- downloadHandler(
  filename = sprintf("Star_Charts_%s.pdf", input$plotStarMarkerFeatureIn),
  content = function(file) {
    pdf(file, width = 12, height = 8)
    print(reactiveVals$starMarkerCluster)
    dev.off()
  }
)

output$downloadPlotAbundance <- downloadHandler(
  filename = function() {
    paste0("Population_Abundances", ".pdf")
  },
  content = function(file) {
    ggsave(
      file,
      plot = reactiveVals$abundanceCluster,
      width = 12,
      height = 6
    )
  }
)

output$downloadPlotDensity <- downloadHandler(
  filename = function() {
    paste0("Cluster_Expression", ".pdf")
  },
  content = function(file) {
    ggsave(
      file,
      plot = reactiveVals$exprsCluster,
      width = 16,
      height = 12
    )
  }
)

output$downloadPlotFrequency <- downloadHandler(
  filename = "Cluster_Heatmap.pdf",
  content = function(file) {
    pdf(file, width = 12, height = 8)
    draw(reactiveVals$heatmapCluster)
    dev.off()
  }
)
