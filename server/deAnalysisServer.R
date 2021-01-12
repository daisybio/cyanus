library(diffcyt)

plotbox_height <- "45em"
methods_height <- "35em"

# checks which methods is selected and executes the diffcyt function accordingly
call_diffcyt <- function(){
  ei <- metadata(reactiveVals$sce)$experiment_info
  #contrast <- createContrast(c(0, 1))  # TODO: contrast matrix must be reactive
  contrastVars <- isolate(input$contrastVars)
  nr_samples <- length(levels(colData(reactiveVals$sce)$sample_id))
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DA-voom")){
    
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatix(contrastVars, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      design = design,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
    }
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-limma")){
    
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatix(contrastVars, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      design = design,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
    }
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-LMM")){
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    contrast <- createCustomContrastMatix(contrastVars, input$colsFixed, designMatrix = F)
    if(nrow(contrast) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
    }
  } else if (input$chosenDAMethod %in% c("diffcyt-DA-GLMM")){
    
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    contrast <- createCustomContrastMatix(contrastVars, input$colsFixed, designMatrix = F)
    if(nrow(contrast) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
    )
    }
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
    shinyjs::hide("heatmapBox")
    shinyjs::hide("deTopTable")
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
  colsDesign <- colsDesign[!colsDesign %in% c("n_cells", "sample_id")]
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
  cols <- cols[!cols %in% c("n_cells", "sample_id")]
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

# choose method and parameter box
output$visDiffExp <- renderUI({
  if (reactiveVals$methodType == "DS") {
    sort_by <-  c("P-adjusted" = "padj", "None" = "none")
  } else {
    sort_by <- c("P-adjusted" = "padj",
                 "LogFC" = "lfc",
                 "None" = "none")
  }
  
 shinydashboard::box(
  selectizeInput(
    inputId = "deVisMethod",
    label = "Successful Run",
    choices = names(reactiveVals$DAruns),
    multiple = F
  ),
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
    choices = c(
      "P-adjusted" = "padj",
      "LogFC" = "lfc",
      "None" = "none"
    ),
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
  if(!is.null(out)){
    # add method to DAruns
    reactiveVals$DAruns[[DAmethod]] <- out
  
    # other method can be performed
    updateButton(session,
               "diffExpButton",
               label = " Start Analysis",
               disabled = FALSE)
  
    shinyjs::show("visDiffExp")
  }else{
    updateButton(session,
                 "diffExpButton",
                 label = " Start Analysis",
                 disabled = FALSE)
  }
})

# if Visualize Differential Analysis Button is clicked -> plotDiffHeatmap is called
observeEvent(input$visExpButton,{
  reactiveVals$visMethod <- isolate(input$deVisMethod)
  reactiveVals$fdrThreshold <- isolate(input$fdrThreshold)
  reactiveVals$lfcThreshold <- isolate(input$lfcThreshold)
  reactiveVals$heatmapSortBy <- isolate(input$heatmapSortBy)
  reactiveVals$heatmapNormalize <- isolate(input$heatmapNormalize)
  shinyjs::show("heatmapBox")
  shinyjs::show("deTopTable")
})


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
  out <- reactiveVals$DAruns[[reactiveVals$visMethod]] 
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


### HEATMAP FUNCTIONS

# Box including Heatmap
output$heatmapBox <- renderUI({
  shinydashboard::box(
    shinycssloaders::withSpinner(
      plotOutput("heatmapDEPlot", width = "100%", height = "550px")
    ),
    title = "Heatmap",
    width = 12,
    height = plotbox_height
  )
})

# Render Heatmap Plot
output$heatmapDEPlot <- renderPlot({
  out <- reactiveVals$DAruns[[reactiveVals$visMethod]] 
  reactiveVals$diffHeatmapPlot <- plotDiffHeatmap(
    x=reactiveVals$sce,
    y=rowData(out$res), 
    fdr=as.numeric(reactiveVals$fdrThreshold), 
    lfc=as.numeric(reactiveVals$lfcThreshold), 
    sort_by = reactiveVals$heatmapSortBy, 
    normalize=as.logical(reactiveVals$heatmapNormalize ),
  )
  reactiveVals$diffHeatmapPlot
})

# ui for download button
output$heatmapPlotDownload <- renderUI({
  req(reactiveVals$diffHeatmapPlot)
  downloadButton("downloadPlotDiffHeatmap", "Download Plot")
})

# function for downloading MDS plot
output$downloadPlotDiffHeatmap <- downloadPlotFunction("Heatmap_DiffExp_Plot", reactiveVals$diffHeatmapPLot)


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








