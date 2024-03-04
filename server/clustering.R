
resetClustering <- function(){
  reactiveVals$clusterFeatureNames <- NULL
  reactiveVals$mergingFrame <- NULL
  reactiveVals$starMarkerCluster <- NULL
  reactiveVals$starCluster <- NULL
  reactiveVals$abundanceCluster <- NULL
  reactiveVals$exprsCluster <- NULL
  reactiveVals$medianExprsCluster <- NULL
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

observeEvent(input$dropClusterButton, {
  clusters_to_drop <- isolate(input$selectedClustersToDrop)
  if(length(clusters_to_drop) == length(levels(cluster_codes(isolate(reactiveVals$sce))[, input$clusterCode]))){
    showNotification(HTML("You cannot drop all clusters!"), type = 'error')
    return(NULL)
  }
  resetVisualization()
  resetDE()
  sapply(reducedDimNames(reactiveVals$sce), function(dr){reducedDim(reactiveVals$sce, dr) <- NULL})
  reactiveVals$sce <- 
    evap(expression(reactiveVals$sce <- dropClusters(reactiveVals$sce, clusters_to_drop))[[1]],
                           params = list(clusters_to_drop = clusters_to_drop))
  msg <- paste0("Dropped cluster ",  paste0(isolate(input$selectedClustersToDrop), collapse = ', '), "! Note that the MST from the Star Charts is now no longer a spanning tree and that the plots from the metaclustering analysis are still the ones from the original cluster run.")
  showNotification(HTML(msg), type = "warning", duration = NULL)
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
             style = "float: left; margin-bottom: 10px;"
           ),
           width = 12,
           style = "overflow-x: scroll;"),
    fluidRow(withSpinner(uiOutput("metaClusteringAnalysis"))),
    fluidRow(withSpinner(uiOutput(
      "clusteringOutput"
    ))),
    fluidRow(withSpinner(uiOutput(
      "clusterModTabs"
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

output$metaClusteringAnalysis <- renderUI({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  runjs(
    "document.getElementById('clusteringVisualizationSelection').scrollIntoView();"
  )
  shinydashboard::box(
    fluidRow(
      shinydashboard::tabBox(
        tabPanel(
          "ECDF",
          div(
            HTML(
              'ConsensusClusterPlus re-samples from the original data and checks how often two cells end up in the same cluster. If the clustering is perfect, the resulting consensus index will only be zero (the two cells never end up in the same cluster) or one (the two cells always end up in the same cluster). This plot shows the cumulative distribution functions of the consensus indices for all metaclusters. <b> An optimal curve is hence horizontal between 0 and 1</b>. For more information, see <a href="https://doi.org/10.1023/A:1023949509487" target="_blank">Monti et al., 2003</a> or <a href="https://doi.org/10.1093/bioinformatics/btq170" target="_blank">Matthew D. Wilkerson, D. Neil Hayes, ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking, Bioinformatics, Volume 26, Issue 12, June 2010, Pages 1572–1573.</a>'
            ),
            style = "text-align: center; vertical-align: middle;"
          ),
          fluidRow(withSpinner(
            plotOutput('ecdf',
                       height = "800px")
          ))
        ),
        tabPanel(
          "PAC",
          div(
            HTML(
              'The PAC (Proportion of Ambiguous Clusters) is defined as the fraction of cell pairs with consensus indices between <it>x1, x2</it>. Here, x1=0.05 and x2=0.95. The ECDF plot has the consensus index values on the x-axis and CDF values on the y-axis. Since <it>CDF(c)</it> corresponds to the fraction of cell pairs with consensus index values &le; c, the PAC is given by <it>CDF(x2) - CDF(x1)</it>. <b>A low value of PAC indicates a flat middle segment. Hence, the optimal metacluster has the lowest PAC. The maximum curvature indicates the point with the highest change in slope, which could be the elbow point. </b>For more details, see <a href="https://doi.org/10.1038/srep06207" target="_blank">Șenbabaoğlu, Y., Michailidis, G. & Li, J. Critical limitations of consensus clustering in class discovery. Sci Rep 4, 6207 (2014).</a>'
              ),
            style = "text-align: center; vertical-align: middle;"
          ),
          fluidRow(withSpinner(
            plotlyOutput('pac',
                         height = "800px")
          ))
        ),
        tabPanel(
        "Delta Area",
        div(
          HTML(
            'The delta area represents the amount of extra cluster stability gained when clustering into k groups as compared to k-1 groups. It is the difference between the ECDF curves for metacluster k and k-1. It can be expected that high stability of clusters can be reached when clustering into the number of groups that best fit the data. <b>It is recommended to choose a meta-cluster where the plateau is reached, similarly to the <it>elbow method</it>. The maximum curvature indicates the point with the highest change in slope, which could be the elbow point</b>. The <it>natural</it> number of clusters present in the data should thus correspond to the value of k where there is no longer a considerable increase in stability (plateau onset)." Crowell et al. (2020)'
          ),
          style = "text-align: center; vertical-align: middle;"
        ),
        fluidRow(withSpinner(
          plotlyOutput('delta_area',
            height = "800px")
        ))
      ),
      title = "Metaclustering Visualization",
      id = "metaclusterVisTabBox",
      width = 12)
    ),
    title = "1. Metaclustering Analysis",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
  )
})

output$pac <- renderPlotly({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  plot_pac(reactiveVals$sce, interactive = TRUE)
})

output$delta_area <- renderPlotly({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  plot_delta_area(reactiveVals$sce, interactive = TRUE)
})

output$ecdf <- renderPlot({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun)
  
  plot_ecdf(reactiveVals$sce, interactive = FALSE, pal = reactiveVals$selected_palette)
})


output$clusterModTabs <- renderUI({
  shinydashboard::box(
    tabBox(
      tabPanel(
        uiOutput("clusterMergingBox"),
        value = "mergeClustersTab",
        title = "Merge Clusters"
      ),
      tabPanel(
        uiOutput("clusterDropBox"),
        value = "dropClustersTab",
        title = "Drop Clusters"
      ),
      id = "Modify Clusters",
      title = "Merge or delete clusters",
      width = 12
    ),
    title = "3. Modify Clusters",
    width = 12,
    collapsible = TRUE,
    collapsed = TRUE
  )
})

output$clusterMergingBox <- renderUI({
  req(reactiveVals$sce$cluster_id, cluster_codes(reactiveVals$sce), metadata(reactiveVals$sce)$clusterRun, input$clusterCode)
  
  div(
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
    )
  )
})

output$clusterDropBox <- renderUI({
  choices <- levels(CATALYST::cluster_codes(reactiveVals$sce)[, input$clusterCode])
  div(
    HTML(
      "You can drop whole clusters here, e.g., a cluster you identified to only contain dead cells."
    ),
    shinyWidgets::pickerInput(
      inputId = "selectedClustersToDrop",
      label = "Select which clusters you want to drop",
      choices = choices,
      selected = choices[1],
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 3",
        "dropup-auto" = T
      )
    ),
    div(
      bsButton(
        "dropClusterButton",
        "Drop Selected Clusters",
        icon("filter"),
        style = "warning"
      ),
      style = "margin-top: 10px; margin-bottom: 10px; float: right;"
    )
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
          "Marker Expression",
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
              selectizeInput(
                "medianClusterExpressionScale",
                "Scale:",
                c("never", "first", "last"),
                multiple = F
              ),
              circle = TRUE,
              status = "info",
              icon = icon("gear"),
              width = "400px",
              tooltip = tooltipOptions(title = "Click to see plot options")
            ),
            style = "position: relative; z-index: 99; float: left;"
          ),
          div(HTML("Heatmap of median marker expressions aggregated by cluster"), style = "text-align: center;vertical-align: middle;"),
          div(uiOutput("clusterMedianDownload"),
              style = "position: relative; z-index: 99; float: right;"),
          fluidRow(withSpinner(
            plotOutput("clusterMedianExprsPlot",
                       height = "800px")
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

output$clusterMedianDownload <- renderUI({
  req(reactiveVals$medianExprsCluster)
  
  downloadButton("downloadPlotMedian", "Download Plot")
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
  custom_colors <- reactiveVals$selected_palette
  if(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])) > length(reactiveVals$selected_palette)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])))
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
  custom_colors <- reactiveVals$selected_palette
  if(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])) > length(reactiveVals$selected_palette)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(length(levels(cluster_codes(reactiveVals$sce)[[input$clusterCode]])))
  }
  reactiveVals$starMarkerCluster <-
    plotMarkerCustom(
      reactiveVals$sce,
      input$plotStarMarkerFeatureIn,
      facet_by = input$plotStarMarkerFacets,
      subselection_col = input$plotStarMarkerSubselection,
      subselection = input$plotStarMarkerSubselectionChoices,
      assayType = metadata(reactiveVals$sce)$clusterRun$assayType,
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
      k_pal = reactiveVals$selected_palette
    )
  reactiveVals$abundanceCluster
})

output$clusterHeatmapPlot <- renderPlot({
  req(nlevels(cluster_ids(reactiveVals$sce, input$clusterCode)) > 1)
  reactiveVals$heatmapCluster <-
    plotFreqHeatmapCustom(reactiveVals$sce,
                          input$clusterCode,
                          k_pal = reactiveVals$selected_palette)
  reactiveVals$heatmapCluster
})

output$clusterMedianExprsPlot <- renderPlot({
  req(nlevels(cluster_ids(reactiveVals$sce, input$clusterCode)) > 1)
  features_tmp <- input$densityUseFeatureChoicesIn
  if (features_tmp == "all")
    features_tmp <- NULL
  reactiveVals$medianExprsCluster <-
    plotExprHeatmapCustom(
      reactiveVals$sce,
      features = features_tmp,
      by = 'cluster_id',
      k = input$clusterCode,
      assay = input$assayTypeVisIn,
      scale = input$medianClusterExpressionScale
    )
  reactiveVals$medianExprsCluster
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
    scale_color_manual(values = reactiveVals$selected_palette)+
    scale_fill_manual(values = reactiveVals$selected_palette)
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

output$downloadPlotMedian <- downloadHandler(
  filename = "Median_Cluster_Expression.pdf",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    pdf(file, width = 12, height = 8)
    draw(reactiveVals$medianExprsCluster)
    dev.off()
    waiter_hide(id = "app")
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
