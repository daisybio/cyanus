library(diffcyt)

plotbox_height <- "48em"
methods_height <- "40em"

# checks which methods is selected and executes the diffcyt function accordingly
call_diffcyt <- function(){
  ei <- metadata(reactiveVals$sce)$experiment_info
  
  contrastVars <- isolate(input$contrastVars)
  nr_samples <- length(levels(colData(reactiveVals$sce)$sample_id))


  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR")){
    if(input$normalizeDE == "Yes"){
      normalize <- TRUE
    }else{
      normalize <- FALSE
    }

    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatrix(contrastVars, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
      reactiveVals$methodsInfo[["diffcyt-DA-edgeR"]] <- data.frame(
        analysis_type = "Differential Abundance", 
        method = "edgeR",
        designmatrix = toString(isolate(input$colsDesign)),
        conditions = toString(contrastVars),
        clustering = toString(isolate(input$deCluster)),
        trend_method = toString(isolate(input$edgeR_trendMethod)),
        normalize = toString(isolate(input$normalizeDE))
      )
      
      out <- diffcyt::diffcyt(
        d_input = reactiveVals$sce,
        design = design,
        contrast = contrast,
        analysis_type = reactiveVals$methodType,
        method_DA = input$chosenDAMethod,
        clustering_to_use = input$deCluster,
        normalize = normalize,
        trend_method = input$edgeR_trendMethod
      )
    }
  } else if(input$chosenDAMethod %in% c("diffcyt-DA-voom")){
    if(input$normalizeDE == "Yes"){
      normalize <- TRUE
    }else{
      normalize <- FALSE
    }
    if(input$blockID_voom != ""){
      blockID <- metadata(reactiveVals$sce)$experiment_info[[input$blockID_voom]]
    }else{
      blockID <- NULL
    }
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatrix(contrastVars, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else if(input$blockID_voom %in% input$colsDesign){
      showNotification("Please don't put your blocking variable in the design matrix. See our tooltip for more information", type = "error")
      out <- NULL
    }else{
      reactiveVals$methodsInfo[["diffcyt-DA-voom"]] <- data.frame(
        analysis_type = "Differential Abundance", 
        method = "Voom",
        designmatrix = toString(isolate(input$colsDesign)),
        conditions = toString(contrastVars),
        clustering = toString(isolate(input$deCluster)),
        block_id = toString(isolate(input$blockID_voom)),
        normalize = toString(isolate(input$normalizeDE))
      )
        out <- diffcyt::diffcyt(
          d_input = reactiveVals$sce,
          design = design,
          contrast = contrast,
          analysis_type = reactiveVals$methodType,
          method_DA = input$chosenDAMethod,
          clustering_to_use = input$deCluster,
          normalize = normalize,
          block_id = blockID
        )
    }
  }else if (input$chosenDAMethod %in% c("diffcyt-DS-limma")){
    if(input$blockID_limma != ""){
      blockID <- metadata(reactiveVals$sce)$experiment_info[[input$blockID_limma]]
    }else{
      blockID <- NULL
    }
    if(input$trend_limma == "Yes"){
      trend <- TRUE
    }else{
      trend <- FALSE
    }
    
    markersToTest <- input$DEFeaturesIn
    is_marker <- rowData(reactiveVals$sce)$marker_class %in% c("type", "state")
    if (input$DEMarkerToTest == "Marker by Class") {
      markersToTest <- (rowData(reactiveVals$sce)$marker_class %in% markersToTest)[is_marker] # type and state (but not none)
    }else{
      markersToTest <- rownames(reactiveVals$sce)[is_marker] %in% markersToTest
    }
    
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatrix(contrastVars, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else if(input$blockID_limma %in% input$colsDesign){
      showNotification("Please don't put your blocking variable in the design matrix. See our tooltip for more information", type = "error")
      out <- NULL
    }else{
      reactiveVals$methodsInfo[["diffcyt-DS-limma"]] <- data.frame(
        analysis_type = "Differential States", 
        method = "limma",
        designmatrix = toString(isolate(input$colsDesign)),
        conditions = toString(contrastVars),
        clustering = toString(isolate(input$deCluster)),
        features = toString(isolate(input$DEFeaturesIn)),
        block_id = toString(isolate(input$blockID_limma)),
        trend_method = toString(isolate(input$trend_limma))
      )
      out <- diffcyt::diffcyt(
        d_input = reactiveVals$sce,
        design = design,
        contrast = contrast,
        analysis_type = reactiveVals$methodType,
        method_DS = input$chosenDAMethod,
        clustering_to_use = input$deCluster,
        block_id = blockID,
        trend = trend,
        markers_to_test = markersToTest
      )
    }
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-LMM")){
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    contrast <- createCustomContrastMatrix(contrastVars, diffcyt::createDesignMatrix(ei, cols_design = input$colsFixed), designMatrix = T)
    
    markersToTest <- isolate(input$DEFeaturesIn)
    is_marker <- rowData(reactiveVals$sce)$marker_class %in% c("type", "state")
    if (input$DEMarkerToTest == "Marker by Class") {
      markersToTest <- (rowData(reactiveVals$sce)$marker_class %in% markersToTest)[is_marker] # type and state (but not none)
    }else{
      markersToTest <- rownames(reactiveVals$sce)[is_marker] %in% markersToTest
    }
    
    if(nrow(contrast) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
      reactiveVals$methodsInfo[["diffcyt-DS-LMM"]] <- data.frame(
        analysis_type = "Differential States", 
        method = "LMM",
        fixed_effects = toString(isolate(input$colsFixed)),
        random_effects = toString(isolate(input$colsRandom)),
        conditions = toString(contrastVars),
        clustering = toString(isolate(input$deCluster)),
        features = toString(isolate(input$DEFeaturesIn))
      )
      
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      markers_to_test = markersToTest,
    )
    }
  } else if (input$chosenDAMethod %in% c("diffcyt-DA-GLMM")){
    if(input$normalizeDE == "Yes"){
      normalize <- TRUE
    }else{
      normalize <- FALSE
    }
    
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    contrast <- createCustomContrastMatrix(contrastVars, diffcyt::createDesignMatrix(ei, cols_design = input$colsFixed), designMatrix = T)
    if(nrow(contrast) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      out <- NULL
    }else{
      reactiveVals$methodsInfo[["diffcyt-DA-GLMM"]] <- data.frame(
        analysis_type = "Differential Abundance", 
        method = "GLMM",
        fixed_effects = toString(isolate(input$colsFixed)),
        random_effects = toString(isolate(input$colsRandom)),
        conditions = toString(contrastVars),
        clustering = toString(isolate(input$deCluster)),
        normalize = toString(isolate(input$normalizeDE))
      )
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      normalize = normalize
    )
    }
  }
}

createCustomContrastMatrix <- function(contrastVars, matrix, designMatrix = T){
  if(designMatrix){
    #the entries have to correspond to the columns of the design matrix
    cnames <- colnames(matrix)
    bool <- getBools(cnames, contrastVars)
    contrast <- createContrast(bool)
    print(contrast)
    return(contrast)
  }else{
    #the entries have to correspond to the levels of the fixed effect terms in the model formula
    lvlList <- lapply(matrix, function(x){levels(colData(reactiveVals$sce)[[x]])})
    names(lvlList) <- matrix
    bool <- getBools(matrix, contrastVars)
    names(bool) <- matrix
    contrast <- unlist(lapply(names(lvlList), function(x){
      return( c( 0, rep(bool[x], length(lvlList[[x]]) -1 )) ) 
      }))
    print(contrast)
    return(createContrast(unname(contrast)))
  }
}

getBools <- function(names, contrastVars){
  bool <- unlist(lapply(names, function(x){
    any(lapply(contrastVars, function(y){
      grepl(y,x, fixed = T )
    }))
  }))
  bool <- as.numeric(bool)
  return(bool)
}

## RENDERER

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
  clusters <- rev(names(cluster_codes(reactiveVals$sce)))
  if(reactiveVals$methodType == "DA"){
    clusters <- clusters[!clusters %in% c("all")]
  }
  selectizeInput(
    inputId = "deCluster",
    label = "Choose the cluster populations you want to compare",
    choices = clusters, 
    multiple = F
  )
})

output$normalizeSelection <- renderUI({
  req(input$chosenDAMethod)
  if(input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DA-GLMM","diffcyt-DA-voom")){
    div(
      radioButtons(
        inputId = "normalizeDE",
        label = span("Normalize?", icon("question-circle"), id = "normalizeDEQ"),
        choices = c("Yes", "No"),
        inline = T
      ),
      bsPopover(
        id = "normalizeDEQ",
        title = "Composition Effects",
        content = "Whether to include optional normalization factors to adjust for composition effects. Only relevant for Differential Abundance methods."
      )
    )
  }
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
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR", "diffcyt-DS-limma", "diffcyt-DA-voom")) {
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

output$extraFeatures <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod == "diffcyt-DA-edgeR") {
    div(
      selectizeInput(
        inputId = "edgeR_trendMethod",
        label = span("Trend Method", icon("question-circle"), id = "trendMethodQ"),
        choices = c("locfit", "none", "movingave", "loess", "locfit.mixed"),
        multiple = F
      ),
      bsPopover(
        id = "trendMethodQ",
        title = "Method for estimating dispersion trend",
        content = "EdgeR specific parameter for estimating the dispersion trend. See more in their estimateDisp function documentation.",
        placement = "top"
      )
    )
  } else if (input$chosenDAMethod == "diffcyt-DA-voom") {
    cols <- colnames(metadata(reactiveVals$sce)$experiment_info)
    cols <- cols[!cols %in% c("n_cells", "sample_id")]
    div(
      selectizeInput(
        inputId = "blockID_voom",
        label = span("Block ID", icon("question-circle"), id = "blockIDvoomQ"),
        choices = cols,
        options = list(
          placeholder = "Select your block IDs or nothing",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = F
      ),
      bsPopover(
        id = "blockIDvoomQ",
        title = "Block IDs if you have replicates",
        content = "Block IDs (e.g. patient IDs if you have replicates in the same conditions) for paired experimental designs, to be included as random effects (for method testDA_voom or testDS_limma). If provided, the block IDs will be included as random effects using the limma duplicateCorrelation methodology. Alternatively, block IDs can be included as fixed effects in the design matrix.",
        placement = "top"
      )
    )
  } else if (input$chosenDAMethod == "diffcyt-DS-limma") {
    cols <- colnames(metadata(reactiveVals$sce)$experiment_info)
    cols <- cols[!cols %in% c("n_cells", "sample_id")]
    div(
      selectizeInput(
        inputId = "blockID_limma",
        label = span("Block ID", icon("question-circle"), id = "blockIDlimmaQ"),
        choices = cols,
        options = list(
          placeholder = "Select your block IDs or nothing",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = F
      ),
      bsPopover(
        id = "blockIDlimmaQ",
        title = "Block IDs if you have replicates",
        content = "Block IDs (e.g. patient IDs if you have replicates in the same conditions) for paired experimental designs, to be included as random effects (for method testDA_voom or testDS_limma). If provided, the block IDs will be included as random effects using the limma duplicateCorrelation methodology. Alternatively, block IDs can be included as fixed effects in the design matrix.",
        placement = "top"
      ),
      radioButtons(
        inputId = "trend_limma",
        label = span("Limma Trend Method", icon("question-circle"), id = "trendlimmaQ"),
        choices = c("Yes", "No"),
        inline = T
      ),
      bsPopover(
        id = "trendlimmaQ",
        title = "Limma Trend Method",
        content = "Whether to fit a mean-variance trend when calculating moderated tests with function eBayes from limma package (for method testDS_limma). When trend = TRUE, this is known as the limma-trend method (Law et al., 2014; Phipson et al., 2016).",
        placement = "top"
      )
    )
  }
})

# if diffcyt should be exectued on selected markers (markers_to_test)
output$markerToTestSelection <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod %in% c("diffcyt-DS-limma", "diffcyt-DS-LMM")) {
    div(
      selectInput(
        "DEMarkerToTest",
        label = "Features to choose from",
        choices = c("Marker by Class",
                    "Marker by Name")
      ),
      uiOutput("DEFeatureSelection")
    )
  }
})

# pickerinput with all markers (for markers_to_test)
output$DEFeatureSelection <- renderUI({
  req(input$DEMarkerToTest)
  if (input$DEMarkerToTest == "Marker by Class") {
    choices <-
      levels(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
    choices <- choices[choices %in% c("type", "state")]
    if("state" %in% choices){
      selected <- "state"
    }else{
      selected <- choices[1]
    }
  } else if (input$DEMarkerToTest == "Marker by Name") {
    is_marker <- rowData(reactiveVals$sce)$marker_class %in% c("type", "state")
    choices <- rownames(reactiveVals$sce)
    names(choices) <-
      sprintf("%s (%s)", choices, as.character(marker_classes(reactiveVals$sce)))
    if("state" %in% marker_classes(reactiveVals$sce)){
      selected <-
        rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) == "state"]
      choices <- sortMarkerNames(choices, as.character(marker_classes(reactiveVals$sce)), first = "state")
    }else{
      selected <- rownames(reactiveVals$sce)[marker_classes(reactiveVals$sce) == "type"]
      choices <- sortMarkerNames(choices, as.character(marker_classes(reactiveVals$sce)), first = "type")
    }
    choices <- choices[is_marker]
  } else
    stop("by name or by class?")
  shinyWidgets::pickerInput(
    inputId = "DEFeaturesIn",
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

# fluidRow with visualization parameter and heatmap
output$DEVisualization <- renderUI({
  shinydashboard::box(
    column(
      width = 3,
      div(uiOutput("visDiffExp"),
          style = "position: relative; height: 500px;"),
    ),
    column(
      width = 8,
      shinycssloaders::withSpinner(plotOutput(
        "heatmapDEPlot", width = "100%", height = "550px"
      )),
    ),
    column(
      width = 1,
      uiOutput("infoDE")
    ),
    title = "Visualize Differential Expression Results",
    height = plotbox_height,
    width = 12,
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

output$topNSelection <- renderUI({
  out <- reactiveVals$DAruns[[input$deVisMethod]]
  methodsDA <- c("diffcyt-DA-edgeR","diffcyt-DA-voom","diffcyt-DA-GLMM")
  if(input$deVisMethod %in% methodsDA){
    label <- "Number of top clusters to display:"
  } else {
    label <- "Number of top cluster-marker combinations to display:"
  }
  
  div(
    numericInput(
      "topN",
      label = label,
      value = min(20,nrow(rowData(out$res))),
      min = 1,
      max = nrow(rowData(out$res)),
      step = 1
    ), 
  )
})

# choose method and parameter box
output$visDiffExp <- renderUI({
  runs <- names(reactiveVals$DAruns)
  div(
    selectizeInput(
      inputId = "deVisMethod",
      label = "Successful Run",
      choices = runs,
      selected = runs[1],
      multiple = F
    ),
    uiOutput("visSelection"),
    uiOutput("topNSelection"),
    numericInput(
      "fdrThreshold",
      label = "FDR Threshold",
      value = 0.05,
      step = 0.05
    ),
    numericInput(
      "lfcThreshold",
      label = "Log2FC Threshold",
      value = 1,
      step = 0.5
    ),
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
        "Visualize Results",
        icon = icon("palette"),
        style = "success"
      ),
      style = "float: left;"
    ),
    div(uiOutput("heatmapPlotDownload"),
        style = "float:right;"
    ),
    
    style = "position: relative; height: 500px;"
  )
})


## OBSERVER

observe({
  if (reactiveVals$current_tab==6){
    shinyjs::hide("DEVisualization")
  }
})

observeEvent(input$deBoxFacet, {
  if(input$deBoxFacet == "cluster_id"){
    shinyjs::show("deBoxK")
  }else{
    shinyjs::hide("deBoxK")
  }
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
    output$heatmapDEPlot <- renderPlot({
      ggplotObject <- ggplot() + theme_void()
      return(ggplotObject)
    })
  }
  
  if(is.null(reactiveVals$methodsInfo)){
    reactiveVals$methodsInfo <- list()
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
  
    shinyjs::show("DEVisualization")
    shinyjs::enable("visExpButton")
  }else{
    updateButton(session,
                 "diffExpButton",
                 label = " Start Analysis",
                 disabled = FALSE)
  }
})

# if Visualize Differential Analysis Button is clicked -> plotDiffHeatmap is called
observeEvent(input$visExpButton,{
  visMethod <- isolate(input$deVisMethod)
  fdrThreshold <- isolate(input$fdrThreshold)
  lfcThreshold <- isolate(input$lfcThreshold)
  heatmapSortBy <- isolate(input$heatmapSortBy)
  heatmapNormalize <- isolate(input$heatmapNormalize)
  deCluster <- isolate(input$deCluster)
  runs <- isolate(reactiveVals$DAruns)
  topN <- isolate(input$topN)

  methodsDA <- c("diffcyt-DA-edgeR","diffcyt-DA-voom","diffcyt-DA-GLMM")
  
  if(visMethod %in% methodsDA){
    heatmapSelection <- isolate(input$DEClusterSelection)
  } else {
    heatmapSelection <- isolate(input$DEMarkerSelection)
  }
  
  shinyjs::show("heatmapDEPlot")
  shinyjs::show("heatmapPlotDownload")
  
  ### HEATMAP FUNCTIONS
  
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
    
    out <- runs[[visMethod]]
    
    reactiveVals$diffHeatmapPlot <- plotDiffHeatmap(
      x=x,
      y=rowData(out$res), 
      fdr=as.numeric(fdrThreshold), 
      lfc=as.numeric(lfcThreshold), 
      sort_by = heatmapSortBy, 
      normalize=as.logical(heatmapNormalize ),
      all = TRUE,
      top_n = as.numeric(topN),
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
  
  ## Info button
  output$infoDE <- renderUI({
    table <- isolate(reactiveVals$methodsInfo[[visMethod]])
    value <- renderTable(
      checkNullTable(table),
      caption.placement = "top"
    )
    
    dropdownButton(
      shinydashboard::box(value, title = "Info", width = 12),
      icon = icon("info-circle"),
      status = "info",
      right = TRUE
    )
  })
  
  ## TOP TABLE FUNCTIONS
  output$deTopTable <- renderUI({
    shinydashboard::box(
      dataTableOutput("topTable"),
      div(
        downloadButton("downloadTopTable", "Download Table Results"),
        style = "float: right; bottom:5px"
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
}, ignoreNULL = FALSE, ignoreInit=TRUE)

# if no markers to test are selected -> analysis cant be performed
observeEvent(input$DEFeaturesIn,{
  if (length(input$DEFeaturesIn)==0){
    shinyjs:: disable("diffExpButton")
  } else {
    shinyjs::enable("diffExpButton")
  }
  }, ignoreNULL=FALSE, ignoreInit=TRUE)



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
  markers <- unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
  if("state" %in% markers){
    selected_marker <- "state"
  }else{
    selected_marker <- markers[1]
  }
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
                                   markers, 
                                   selected = selected_marker,
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
  if(input$deBoxK == "all"){
    k <- NULL
    facet_by <- NULL
  }else{
    k <- input$deBoxK
    facet_by <- input$deBoxFacet
  }
  reactiveVals$pbExprsPlot <- plotPbExprs(reactiveVals$sce, 
              k = k, 
              features = input$deBoxFeatures, 
              color_by = input$deBoxColor, 
              facet_by = facet_by, 
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








