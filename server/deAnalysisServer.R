library(diffcyt)

plotbox_height <- "45em"
methods_height <- "40em"

# checks which methods is selected and executes the diffcyt function accordingly
call_diffcyt <- function(){
  ei <- metadata(reactiveVals$sce)$experiment_info
  
  contrastVars <- isolate(input$contrastVars)
  
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DA-voom")){
    
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatix(contrastVars, design, designMatrix = T)
    
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      design = design,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-limma")){
    
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatix(contrastVars, design, designMatrix = T)
    
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      design = design,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-LMM")){
    
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    contrast <- createCustomContrastMatix(contrastVars, input$colsFixed, designMatrix = F)
    
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
  } else if (input$chosenDAMethod %in% c("diffcyt-DA-GLMM")){
    
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    contrast <- createCustomContrastMatix(contrastVars, input$colsFixed, designMatrix = F)
    
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
  }
  out
}

createCustomContrastMatix <- function(contrastVars, matrix, designMatrix = T){
  if(designMatrix){
    #the entries have to correspond to the columns of the design matrix
    cnames <- colnames(matrix)
    bool <- getBools(cnames, contrastVars)
    bool <- as.numeric(bool)
    contrast <- createContrast(bool)
    return(createContrast(contrast))
  }else{
    #the entries have to correspond to the levels of the fixed effect terms in the model formula
    lvlList <- lapply(matrix, function(x){levels(colData(reactiveVals$sce)[[x]])})
    names(lvlList) <- matrix
    bool <- getBools(matrix, contrastVars)
    bool <- as.numeric(bool)
    names(bool) <- matrix
    contrast <- unlist(lapply(names(lvlList), function(x){
      return( c(rep(bool[x], length(lvlList[[x]]))) ) 
      }))
    return(createContrast(unname(contrast)))
  }
}

getBools <- function(names, contrastVars){
  bool <- unlist(lapply(names, function(x){
    any(lapply(contrastVars, function(y){
      grepl(y,x, fixed = T )
    }))
  }))
  return(bool)
}

observe({
  if (reactiveVals$current_tab==6){
    shinyjs::hide("visDiffExp")
  }
})

observeEvent(input$deBoxFacet, {
  if(input$deBoxFacet == "cluster_id"){
    shinyjs::show("deBoxK")
  }else{
    shinyjs::hide("deBoxK")
  }
})

# displays available methods and selection of DA or DS
output$deMethodSelection <- renderUI({
   methodsDA <- c("edgeR" = "diffcyt-DA-edgeR", "Voom" = "diffcyt-DA-voom", "GLMM" = "diffcyt-DA-GLMM")
   methodsDS <- c("limma" = "diffcyt-DS-limma","LMM" = "diffcyt-DS-LMM")
   if(input$da_ds == "Differential Abundance"){
     choices <- methodsDA
     reactiveVals$methodType <- "DA"
   }else{
     choices <- methodsDS
     reactiveVals$methodType <- "DS"
   }
   div(
     selectizeInput(
      inputId = "chosenDAMethod",
      label = span("Available Methods", icon("question-circle"), id = "deMethodsQ"),
      choices = choices,
      multiple = F
    ),
    bsPopover(
      id = "deMethodsQ",
      title = "Available Methods",
      content = "Depending on what you want to analyse, there are different methods available. Please see their documentation for further explanation"
    )
   )
})

# box with cluster populations you want to compare
output$clusterSelection <- renderUI({
  selectizeInput(
    inputId = "deCluster",
    label = "Choose the cluster populations you want to compare",
    choices = names(cluster_codes(reactiveVals$sce)), 
    multiple = F
  )
})

# checks whether designMatrix or formula box must be visualized
output$modelSelection <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DS-limma","diffcyt-DA-voom")){
    uiOutput("designMatrixSelection")
  } else {
    uiOutput("formulaSelection")
  }
})

# if method using a design matrix is selected -> cols_design must be specified
output$designMatrixSelection <- renderUI({
  colsDesign <- colnames(metadata(reactiveVals$sce)$experiment_info)
  colsDesign <- colsDesign[!colsDesign %in% "n_cells"]
  div(
    pickerInput(
      "colsDesign",
      choices = colsDesign,
      selected = colsDesign[2],
      label = span(
        "Which columns to include in the design matrix?",
        icon("question-circle"),
        id = "deDesignMatrix"
      ),
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    ),
    bsPopover(
      id = "deDesignMatrix",
      title = "Design matrix for model fitting",
      content = "The selected columns will be included in the design matrix specifcying the models to be fitted. For example, this may include group IDs (e.g. groups for differential testing) or block IDs (e.g. patient IDs in a paired design)."
    )
  )
})

# if method using a formula is selected -> cols_fixed and cols_random must be specified
output$formulaSelection <- renderUI({
  cols <- colnames(metadata(reactiveVals$sce)$experiment_info)
  cols <- cols[!cols %in% "n_cells"]
  div(
    pickerInput(
      "colsFixed",
      choices = cols,
      selected = cols[2],
      label = span(
        "Which fixed effect terms to include in the model formula",
        icon("question-circle"),
        id = "deFormulaFix"
      ),
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    ),
    bsPopover(
      id = "deFormulaFix",
      title = "Fixed effect terms for the model formula",
      content = "Depending on the experimental design, this may include group IDs (e.g. groups for differential testing) or block IDs (e.g. patient IDs in a paired design)."
    ),
    pickerInput(
      "colsRandom",
      choices = cols,
      selected = cols[1],
      label = span(
        "Which random intercept terms to include in the model formula",
        icon("question-circle"),
        id = "deFormulaRandom"
      ),
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    ),
    bsPopover(
      id = "deFormulaRandom",
      title = "Random intercept terms for the model formula",
      content = "Depending on the experimental design, this may include group IDs (e.g. groups for differential testing) or block IDs (e.g. patient IDs in a paired design)."
    ),
  )
})

output$contrastSelection <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DS-limma","diffcyt-DA-voom")){
    choices <- input$colsDesign
  } else {
    choices <- input$colsFixed
  }
  div(
    pickerInput(
      "contrastVars",
      choices = choices,
      selected = choices[1],
      label = span(
        "What condition(s) do you want to analyse?",
        icon("question-circle"),
        id = "deContrastQ"
      ),
      multiple = TRUE
    ),
    bsPopover(
      id = "deContrastQ",
      title = "Contrast Matrix Design",
      content = "Here, you specify the comparison of interest, i.e. the combination of model parameters to test whether they are equal to zero."
    )
  )
})

# heatmap can be visualized for subset of clusters (DA)
output$visClusterSelection <- renderUI({
  out <- reactiveVals$DAruns[[input$deVisMethod]]
  div(
    pickerInput(
      "DEClusterSelection",
      choices = unique(rowData(out$res)$cluster_id),
      selected = unique(rowData(out$res)$cluster_id),
      label = "Visualize results for subset of clusters:",
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    )
  )
})

# heatmap can be visualized for selected features (DS)
output$visMarkerSelection <- renderUI({
  out <- reactiveVals$DAruns[[input$deVisMethod]]
  div(
    pickerInput(
      "DEMarkerSelection",
      choices = as.character(unique(rowData(out$res)$marker_id)),
      selected = as.character(unique(rowData(out$res)$marker_id)),
      label = "Visualize results for subset of markers:",
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    )
  )
})

# check whether cluster selection or marker selection should be displayed
output$visSelection <- renderUI({
  methodsDA <- c("diffcyt-DA-edgeR","diffcyt-DA-voom","diffcyt-DA-GLMM")
  if(input$deVisMethod %in% methodsDA){
    uiOutput("visClusterSelection")
  } else {
    uiOutput("visMarkerSelection")
  }
})

# choose method and parameter box
output$visDiffExp <- renderUI({
 runs <- names(reactiveVals$DAruns)
  
 shinydashboard::box(
  selectizeInput(
    inputId = "deVisMethod",
    label = "Successful Run",
    choices = runs,
    selected = runs[1],
    multiple = F
  ),
  uiOutput("visSelection"),
  numericInput("fdrThreshold",
               label = "FDR Threshold",
               value = 0.05,
               step = 0.05),
  numericInput("lfcThreshold",
               label = "Log2FC Threshold",
               value = 1,
               step = 0.5),
  selectizeInput(
    "heatmapSortBy",
    "Sort by:",
    choices = c("P-adjusted" = "padj", "None" = "none"),
    multiple = F
  ),
  selectizeInput(
    "heatmapNormalize",
    "Z-score normalization:",
    c("Yes" = "TRUE", "No" = "FALSE"),
    multiple = F
  ),
  div(
    bsButton(
      "visExpButton",
      "Visualize Differential Expression",
      icon = icon("palette"),
      style = "success"
    ),
    style = "float: right;"
  ),
  
  title = "Visualize Differential Expression Results",
  width = 6,
  height = methods_height,
 )
})

# if Start Analysis button is clicked -> diffcyt method should be performed
observeEvent(input$diffExpButton,{
  shinyjs:: disable("visExpButton")
  DAmethod <- isolate(input$chosenDAMethod)
 
  # update button and disable it
  updateButton(session,
               "diffExpButton",
               label = " Analysis...",
               disabled = TRUE)
  
  # check if a methods was already performed else create the list
  if (is.null(reactiveVals$DAruns)){
    reactiveVals$DAruns <- list()
  }
  

  # call diffcyt function
  out <- call_diffcyt()
  
  # add method to DAruns
  reactiveVals$DAruns[[DAmethod]] <- out
  
  # other method can be performed
  updateButton(session,
               "diffExpButton",
               label = " Start Analysis",
               disabled = FALSE)
  
  shinyjs::show("visDiffExp")
  shinyjs::enable("visExpButton")
})

# if Visualize Differential Analysis Button is clicked -> plotDiffHeatmap is called
observeEvent(input$visExpButton,{
  visMethod <- isolate(input$deVisMethod)
  fdrThreshold <- isolate(input$fdrThreshold)
  lfcThreshold <- isolate(input$lfcThreshold)
  heatmapSortBy <- isolate(input$heatmapSortBy)
  heatmapNormalize <- isolate(input$heatmapNormalize)
  deCluster <- isolate(input$deCluster)

  methodsDA <- c("diffcyt-DA-edgeR","diffcyt-DA-voom","diffcyt-DA-GLMM")
  
  if(visMethod %in% methodsDA){
    heatmapSelection <- isolate(input$DEClusterSelection)
  } else {
    heatmapSelection <- isolate(input$DEMarkerSelection)
  }
  
  ### HEATMAP FUNCTIONS
  
  # Box including Heatmap
  output$heatmapBox <- renderUI({
    shinydashboard::box(
      shinycssloaders::withSpinner(
        plotOutput("heatmapDEPlot", width = "100%", height = "550px")
      ),
      div(
        uiOutput("heatmapPlotDownload"),
        style = "position: absolute; bottom: 5px; right:5px"
      ),
      title = "Heatmap",
      width = 12,
      height = plotbox_height
    )
  })
  
  # Render Heatmap Plot
  output$heatmapDEPlot <- renderPlot({
    methodsDA <- c("diffcyt-DA-edgeR","diffcyt-DA-voom","diffcyt-DA-GLMM")
    methodsDS <- c("diffcyt-DS-limma","diffcyt-DS-LMM")
    
    if(visMethod %in% methodsDA){
      sub <- filterSCE(reactiveVals$sce, cluster_id %in% heatmapSelection, k=deCluster)
      x <- sub
    } else {
      x <- reactiveVals$sce[rownames(reactiveVals$sce) %in% heatmapSelection, ]
    }
    
    out <- reactiveVals$DAruns[[visMethod]]
    
    reactiveVals$diffHeatmapPlot <- plotDiffHeatmap(
      x=x,
      y=rowData(out$res), 
      fdr=as.numeric(fdrThreshold), 
      lfc=as.numeric(lfcThreshold), 
      sort_by = heatmapSortBy, 
      normalize=as.logical(heatmapNormalize ),
      all = TRUE
    )
    reactiveVals$diffHeatmapPlot
  })
  
  
  # ui for download button
  output$heatmapPlotDownload <- renderUI({
    req(reactiveVals$diffHeatmapPlot)
    downloadButton("downloadPlotDiffHeatmap", "Download Plot")
  })
  
  # function for downloading heatmap
  output$downloadPlotDiffHeatmap <- downloadHandler(
    filename = "DE_Heatmap.pdf", 
    content = function(file){
      pdf(file, width = 12, height = 8)
      draw(reactiveVals$diffHeatmapPlot)
      dev.off()
    }
  )
  
  ## TOP TABLE FUNCTIONS
  output$deTopTable <- renderUI({
    shinydashboard::box(
      dataTableOutput("topTable"),
      div(
        downloadButton("downloadTopTable", "Download Table Results"),
        style = "float: right;"
      ),
      title = "Table of results for top clusters or cluster-marker combinations",
      width = 12,
      height = plotbox_height
    )
    
  })
  
  output$topTable <- renderDataTable({
    out <- reactiveVals$DAruns[[visMethod]] 
    reactiveVals$topTable <- data.frame(diffcyt::topTable(out$res,all=TRUE,format_vals=TRUE))
    DT::datatable(reactiveVals$topTable,
                  rownames = FALSE,
                  options = list(pageLength=10, searching=FALSE))
  })
  
  output$downloadTopTable <- downloadHandler(
    filename = "Differential_Expression_Results.csv",
    content = function(file) {
      write.csv(reactiveVals$topTable, file, row.names = FALSE)
    }
  )
  
})

# check if clusters or markers are selected (to visualize heatmap)
observeEvent({
  input$DEClusterSelection
  input$DEMarkerSelection
}, {
  methodsDA <-
    c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM")
  if (input$deVisMethod %in% methodsDA) {
    if (length(input$DEClusterSelection) == 0) {
      print("disable")
      shinyjs::disable("visExpButton")
    } else {
      shinyjs::enable("visExpButton")
    }
  } else {
    if (length(input$DEMarkerSelection) == 0) {
      shinyjs::disable("visExpButton")
    } else {
      shinyjs::enable("visExpButton")
    }
  }
})


### BOXPLOT FUNCTIONS

output$deBoxPlots <- renderUI({
  uiOutput("deExprsCluster")
  #shinydashboard::tabBox(
  #  tabPanel(
  #    plotOutput("deExprsBoxPlot"), 
  #    title = "Overall marker expression",
  #    value = "deExprsTab", 
  #    with = 12
  #    ),
  #  tabPanel(
  #    uiOutput("deExprsCluster"),
  #    title = "Cluster marker expression",
  #    value = "deExprsClusterTab", 
  #    with = 12
  #  ),
  #  width = 12
  #)
  
})

#output$deExprsBoxPlot <- renderPlot(
#  plotPbExprs(reactiveVals$sce, features = "state", shape_by = "patient_id")
#)

output$deExprsCluster <- renderUI({
  factors <- names(colData(reactiveVals$sce))[!names(colData(reactiveVals$sce)) %in% c("patient_id", "sample_id")]
  fluidRow(column(1,
                  div(dropdownButton(
                    tags$h3("Plot Options"),
                    selectizeInput("deBoxFacet",
                                   "Facet by:",
                                   c("antigen", "cluster_id"), 
                                   multiple = F, 
                                   selected = "antigen"),
                    hidden(selectizeInput(
                      "deBoxK",
                      "Clusters",
                      names(cluster_codes(reactiveVals$sce)),
                      multiple = F,
                      selected = "meta9"
                    )),
                    selectizeInput("deBoxFeatures",
                                   "Markers:",
                                   unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class), 
                                   multiple = F),
                    selectizeInput(
                      "deBoxColor",
                      "Color by:",
                      c(names(colData(reactiveVals$sce)), names(cluster_codes(reactiveVals$sce))),
                      selected = factors[1],
                      multiple = F
                    ),
                    selectizeInput(
                      "deBoxShape",
                      "Shape By: ",
                      c(names(colData(reactiveVals$sce))),
                      multiple = F
                    ),
                    circle = TRUE,
                    status = "info",
                    icon = icon("gear"),
                    width = "400px",
                    tooltip = tooltipOptions(title = "Click to see plot options")
                  ),
                  style = "position: relative; height: 550px;"
                  )
                ),
           column(11, shinycssloaders::withSpinner(
             plotOutput("clusterDEPlot", width = "100%", height = "500px")
           )),
           div(
             uiOutput("pbExprsPlotDownload"),
             style = "position: absolute; bottom: 5px; right:5px"
           ),
           )
})

output$clusterDEPlot <- renderPlot({
  reactiveVals$pbExprsPlot <- plotPbExprs(reactiveVals$sce, 
              k = input$deBoxK, 
              features = input$deBoxFeatures, 
              color_by = input$deBoxColor, 
              facet_by = input$deBoxFacet, 
              shape_by = input$deBoxShape)
  reactiveVals$pbExprsPlot
})

# ui for download button
output$pbExprsPlotDownload <- renderUI({
  req(reactiveVals$pbExprsPlot)
  downloadButton("downloadPlotPbExprs", "Download Plot")
})

# function for downloading MDS plot
output$downloadPlotPbExprs <- downloadPlotFunction("Pb_Exprs_plot", reactiveVals$pbExprsPlot)








