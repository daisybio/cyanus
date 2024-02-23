
resetClustering <- function(){
  reactiveVals$clusterFeatureNames <- NULL
  reactiveVals$mergingFrame <- NULL
  reactiveVals$starMarkerCluster <- NULL
  reactiveVals$starCluster <- NULL
  reactiveVals$abundanceCluster <- NULL
  reactiveVals$exprsCluster <- NULL
  reactiveVals$heatmapCluster <- NULL
}


### Observer ----
observeEvent(input$startClustering, {
  updateButton(session,
               "startClustering",
               label = " Clustering...",
               disabled = TRUE)
  
  waiter_show(id = "app",html = tagList(spinner$logo, 
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
  tryCatch({
    withCallingHandlers({
      reactiveVals$sce <-
        evap(expression(reactiveVals$sce <- clusterSCE(
          reactiveVals$sce,
          assayTypeIn,
          clusterFeatureNames,
          xdim,
          ydim,
          k
        ))[[1]], params = list(assayTypeIn = input$assayTypeIn, clusterFeatureNames = reactiveVals$clusterFeatureNames,
                               xdim = input$xdim, ydim = input$ydim, k = input$k))
    },
    
    error = function(e) {
      shinyjs::html(id = "clusteringProgress",
                    html = sprintf("<br>%s", HTML(e$message)),
                    add = TRUE)
      
      stop(e)
    })
  },
  error = function(e){
    showNotification(HTML(sprintf("Clustering failed with the following message:<br>
                                    <b>%s</b>", e$message)), duration = NULL, type = "error")
    
  })
  
  assays <- c("exprs" = "Transformed", "counts" = "Raw")
  
  metadata(reactiveVals$sce)$clusterRun <- evap(expression(metadata(reactiveVals$sce)$clusterRun <- list(
    features = clusterFeatureNames,
    assayType = assayTypeIn,
    xdim = xdim,
    ydim = ydim,
    maxK = k
  ))[[1]], params = list(clusterFeatureNames = reactiveVals$clusterFeatureNames, assayTypeIn = assays[input$assayTypeIn], 
  xdim = input$xdim, ydim = input$ydim, k = input$k))
  waiter_hide(id = "app")
  
  updateButton(session,
               "startClustering",
               label = " Start Clustering",
               disabled = FALSE)
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


observeEvent({input$featuresIn
  reactiveVals$sce},
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

observeEvent(reactiveVals$clusterFeatureNames, {
  if (length(reactiveVals$clusterFeatureNames) == 0){
    shinyjs::show("invalidClusteringFeaturesWarning")
    disable("startClustering")
  }
  else if (floor(0.8 * (input$xdim * input$ydim)) >= (input$k)){
    shinyjs::hide("invalidClusteringFeaturesWarning")
    enable("startClustering")
  }
  else {
    shinyjs::hide("invalidClusteringFeaturesWarning")
    message("features ok but dimensions?")
  } 
})

observeEvent({
  input$xdim
  input$ydim
  input$k
}, {
  if (floor(0.8 * (input$xdim * input$ydim)) < (input$k)) {
    shinyjs::show("invalidClusteringParamsWarning")
    shinyjs::disable("startClustering")
  } else if (length(reactiveVals$clusterFeatureNames) > 0) {
    shinyjs::hide("invalidClusteringParamsWarning")
    shinyjs::enable("startClustering")
  } else {
    shinyjs::hide("invalidClusteringParamsWarning")
    message("dimensions ok but features?")
  }
})

observeEvent(input$clusterCode, {
  curr_cluster_ids <-
    levels(cluster_ids(x = reactiveVals$sce,
                       k = input$clusterCode))
  reactiveVals$mergingFrame <-
    data.frame(
      old_cluster = curr_cluster_ids,
      new_cluster = curr_cluster_ids
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
  new_metacluster <- isolate(input$mergeClusteringNewName)
  if (is.null(new_metacluster) || is.na(new_metacluster) || new_metacluster == "")
    new_metacluster <- sprintf("merging_%s", isolate(input$clusterCode))
  reactiveVals$sce <-
    evap(expression(reactiveVals$sce <- mergeClusters(
      isolate(reactiveVals$sce),
      clusterCode,
      mergingFrame,
      id = id,
      overwrite = TRUE
    ))[[1]], params = list(clusterCode = isolate(input$clusterCode), mergingFrame = isolate(reactiveVals$mergingFrame), id = new_metacluster))
  
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

observeEvent(reactiveVals$sce, {
  if (!is.null(reactiveVals$sce$cluster_id) &&
      !is.null(cluster_codes(reactiveVals$sce))  &&
      !is.null(metadata(reactiveVals$sce)$clusterRun))
    reactiveVals$continue[which(tab_ids == "clustering")] <- TRUE
  else
    reactiveVals$continue[which(tab_ids == "clustering")] <- FALSE
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
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  shinydashboard::box(
    shinydashboard::box(width=12, 
    div(HTML("Here, you can select the meta-cluster you want to investigate in this part of the analysis.<br>
        The meta-cluster starting with `som` contains the original clusters resulting from the self-organizing map without meta-clustering. You can have a look at it, but you will most likely want to use a meta-cluster identified by ConsensusClusterPlus.<br>
        The meta-clusters starting with `meta` were identified by ConsensusClusterPlus and contain as many clusters as indicated by the number following `meta`. To determine which meta-cluster represents your data best, you can use the delta area plot. (Section 1. Delta Area)<br>
        The meta-cluster `all` contains all cells. It is especially useful for differential expression analysis when comparing conditions overall."))),
    column(
      uiOutput("selectClusterCode"),
      div(
        downloadButton("downloadClusters", "Download Cluster Assignments"),
        style = "float: right;"
      ),
      width = 4
    ),
    column(
      tableOutput("clusterRunParams"),
      width = 8,
      style = "overflow-x: scroll;"
    ),
    column(withSpinner(uiOutput("clusterSizes")),
           width = 12,
           style = "overflow-x: scroll;"),
    column(withSpinner(uiOutput("clusterFrequencies")),
           div(
             downloadButton("downloadFrequencies", "Download Cluster Frequencies per Sample"),
             style = "float: left;"
           ),
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
    metadata(reactiveVals$sce)$clusterRun
  runParams$features <- paste(runParams$features, collapse = ",")
  data.frame(runParams)
},
caption = "Run Parameters",
caption.placement = "top")

output$selectClusterCode <- renderUI({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  choicesClusterCode <- names(cluster_codes(reactiveVals$sce))
  choicesClusterCode <-
    c('all', choicesClusterCode[choicesClusterCode != 'all'])
  selectInput("clusterCode",
              "Select Meta-Cluster",
              rev(choicesClusterCode))
})

output$clusterSizes <- renderTable({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun, input$clusterCode, cluster_ids(reactiveVals$sce, input$clusterCode))
  res <-
    as.data.frame(table(cluster_ids(reactiveVals$sce, input$clusterCode),
                        useNA =
                          "ifany"))
  names(res) <- c("Cluster", "Size")
  
  total_cells <- sum(res$Size)
  res$Fraction <- round((res$Size / total_cells) * 100, 2)
  names(res) <- c("Cluster", "Size", "Proportion [%]")
  t(res)
},
caption = "Cluster Sizes",
caption.placement = "top", rownames = TRUE, colnames = FALSE)

output$clusterFrequencies <- renderTable({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun, input$clusterCode, cluster_ids(reactiveVals$sce, input$clusterCode))
  
  df <- calcClusterFreqBySample(reactiveVals$sce, k = input$clusterCode)
  res <- reshape2::dcast(df, sample_id ~ cluster_id, value.var = "Freq")
  names(res)[1] <- "Sample"
  
  return(res)
},
caption = "Cluster Frequencies per Sample [%]",
caption.placement = "top", rownames = FALSE, colnames = TRUE)

output$mergeClustersDT <- renderDT({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
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
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  runjs(
    "document.getElementById('clusteringVisualizationSelection').scrollIntoView();"
  )
  shinydashboard::box(
    div(HTML(
      'It is recommended to choose a meta-cluster where the plateau is reached, similarly to the `elbow method`.<br>
                          "The delta area represents the amount of extra cluster stability gained when clustering into k groups as compared to k-1 groups.<br>
                          It can be expected that high stability of clusters can be reached when clustering into the number of groups that best fits the data.<br>
                          The `natural` number of clusters present in the data should thus corresponds to the value of k where there is no longer a considerable increase in stability (plateau onset)." Crowell et al. (2020)',
      style = "text-align: center;vertical-align: middle;"
    )),
    renderPlot(
      CATALYST::delta_area(reactiveVals$sce)
    ),
    title = "1. Delta Area",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
  )
})

output$clusterMergingBox <- renderUI({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun, input$clusterCode)
  
  shinydashboard::box(
    HTML(
      "You can assign new cluster names in the <b>new_cluster</b> column by double-clicking.<br>
      <i>If you do not assign new names to all clusters the old cluster names will be kept.</i>"
    ),
    withSpinner(DTOutput("mergeClustersDT")),
    div(
      textInput("mergeClusteringNewName", "New Meta-Cluster Name", value = sprintf("merging_%s", isolate(input$clusterCode))),
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
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    add_meta <- list(cluster_ids(reactiveVals$sce, input$clusterCode))
    names(add_meta) <- input$clusterCode
    
    write.csv(cbind(colData(reactiveVals$sce), add_meta), file, row.names = FALSE)
    waiter_hide(id = "app")
  }
)

output$downloadFrequencies <- downloadHandler(
  filename = "Cluster_Frequencies_Per_Sample.csv",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    
    df <- calcClusterFreqBySample(reactiveVals$sce, k = input$clusterCode)
    res <- reshape2::dcast(df, sample_id ~ cluster_id, value.var = "Freq")
    names(res)[1] <- "Sample"
    
    write.csv(res, file, row.names = FALSE)
    waiter_hide(id = "app")
  }
)

output$clusteringOutput <- renderUI({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  abundanceByChoices <-
    c("Samples" = "sample_id", "Clusters" = "cluster_id")
  abundanceChoices <- names(colData(reactiveVals$sce))
  names(abundanceChoices) <- names(abundanceChoices)
  
  densityAssayChoices <-
    c("Transformed" = "exprs", "Raw" = "counts")
  densityAssayChoices <-
    densityAssayChoices[densityAssayChoices %in% assayNames(reactiveVals$sce)]
  
  sceEI <- ei(reactiveVals$sce)
  starMarkerFacets <- names(which(sapply(sceEI, function(feature) nlevels(as.factor(feature)) == 2)))
  names(starMarkerFacets) <- starMarkerFacets
  starMarkerFacets <- c("No Facetting" = "", starMarkerFacets)
  
  starMarkerSubselection <- colnames(sceEI)
  starMarkerSubselection <- c("No Subselection" = "", starMarkerSubselection[!starMarkerSubselection %in% c('n_cells', 'sample_id')])
  
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
                      droplevels(
                      SummarizedExperiment::rowData(reactiveVals$sce)$marker_class
                      )
                    )
                  ),
                selected = ifelse("type" %in% c(
                  "all",
                  levels(
                    droplevels(
                      SummarizedExperiment::rowData(reactiveVals$sce)$marker_class
                    )
                  )
                ), "type", "state")
              ),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div(HTML("Smoothed densities of marker intensities by cluster. <b style='color:#FF3358';>Attention:</b> This may take some time"), style = "text-align: center;vertical-align: middle;"),
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
              selectInput(
                "plotStarMarkerFacets",
                label = "Facet By",
                choices = starMarkerFacets
              ),
              selectInput(
                "plotStarMarkerSubselection",
                label = "Subselection By",
                choices = starMarkerSubselection
              ),
              uiOutput("plotStarMarkerSubselectionChoicesBox"),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div(
            HTML("Tree, where each node (cluster) is coloured depending on its median value for the given marker.
            <br>If you use subselection, be aware that the MST and the node sizes were computed on the whole dataset."),
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

output$plotStarMarkerSubselectionChoicesBox <- renderUI({
  req(input$plotStarMarkerSubselection, input$plotStarMarkerSubselection != "")
  
  choices <- levels(ei(reactiveVals$sce)[[input$plotStarMarkerSubselection]])
  
  selectInput(
    "plotStarMarkerSubselectionChoices",
    label = "Subselection",
    choices = choices
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
  req(SummarizedExperiment::rowData(reactiveVals$sce)$used_for_clustering)
  custom_colors <- reactiveVals$colorblind_palette
  if(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])) > length(reactiveVals$colorblind_palette)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$colorblind_palette)(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])))
  }
  custom_palette <- RColorBrewer::brewer.pal(length(metadata(reactiveVals$sce)$SOM$map$colsUsed), "Set3")
  if(length(metadata(reactiveVals$sce)$SOM$map$colsUsed) > 12){
    custom_palette <- grDevices::colorRampPalette(custom_palette)(length(metadata(reactiveVals$sce)$SOM$map$colsUsed))
  }
  reactiveVals$starCluster <-
    plotStarsCustom(metadata(reactiveVals$sce)$SOM, overall = TRUE, 
                    colorPalette = custom_palette,
                    backgroundValues = cluster_codes(reactiveVals$sce)[[input$clusterCode]], 
                    backgroundColors = custom_colors)
  print(reactiveVals$starCluster)
  reactiveVals$starCluster
})

output$clusterStarMarkerPlot <- renderPlot({
  custom_colors <- reactiveVals$colorblind_palette
  if(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])) > length(reactiveVals$colorblind_palette)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$colorblind_palette)(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])))
  }
  reactiveVals$starMarkerCluster <-
    plotMarkerCustom(
      reactiveVals$sce,
      input$plotStarMarkerFeatureIn,
      facet_by = input$plotStarMarkerFacets,
      subselection_col = isolate(input$plotStarMarkerSubselection),
      subselection = input$plotStarMarkerSubselectionChoices,
      assayType = names(metadata(reactiveVals$sce)$clusterRun$assayType),
      backgroundValues = cluster_codes(reactiveVals$sce)[[input$clusterCode]],
      backgroundColors = custom_colors
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
      shape_by = shape_by,
      k_pal = reactiveVals$colorblind_palette
    )
  reactiveVals$abundanceCluster
})

output$clusterHeatmapPlot <- renderPlot({
  req(nlevels(cluster_ids(reactiveVals$sce, input$clusterCode)) > 1)
  reactiveVals$heatmapCluster <-
    plotFreqHeatmapCustom(reactiveVals$sce,
                          input$clusterCode,
                          k_pal = reactiveVals$colorblind_palette)
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
    )+
    scale_color_manual(values = reactiveVals$colorblind_palette)+
    scale_fill_manual(values = reactiveVals$colorblind_palette)
  reactiveVals$exprsCluster
})

output$downloadPlotStar <- downloadHandler(
  filename = "Star_Charts_overall.pdf",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    pdf(file, width = 12, height = 8)
    print(reactiveVals$starCluster)
    dev.off()
    waiter_hide(id= "app")
  }
)

output$downloadPlotMarkerStar <- downloadHandler(
  filename = function() {sprintf("Star_Charts_%s.pdf", input$plotStarMarkerFeatureIn)},
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    if (input$plotStarMarkerFacets != "")
      my_width = 16
    else 
      my_width = 12
    pdf(file, width = my_width, height = 8)
    print(reactiveVals$starMarkerCluster)
    dev.off()
    waiter_hide(id= "app")
  }
)

output$downloadPlotAbundance <- downloadHandler(
  filename = function() {
    paste0("Population_Abundances", ".pdf")
  },
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(
      file,
      plot = reactiveVals$abundanceCluster,
      width = 12,
      height = 6
    )
    waiter_hide(id= "app")
  }
)

output$downloadPlotDensity <- downloadHandler(
  filename = function() {
    paste0("Cluster_Expression", ".pdf")
  },
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(
      file,
      plot = reactiveVals$exprsCluster,
      width = 16,
      height = 12
    )
    waiter_hide(id = "app")
  }
)

output$downloadPlotFrequency <- downloadHandler(
  filename = "Cluster_Heatmap.pdf",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    pdf(file, width = 12, height = 8)
    draw(reactiveVals$heatmapCluster)
    dev.off()
    waiter_hide(id = "app")
  }
)
