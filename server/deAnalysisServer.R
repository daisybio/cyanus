source("server/sceEMD.R")
library(diffcyt)

plotbox_height <- "48em"
methods_height <- "40em"

# checks which methods is selected and executes the diffcyt function accordingly

call_DE <- function(){
  contrastVars <- isolate(input$contrastVars)
  if(is.null(input$deSubselection)){
    subselection <- "No"
  }else{
    subselection <- isolate(input$deSubselection)
  }
  sce <- isolate(reactiveVals$sce)
  if(subselection != "No"){
    catCount <- sapply(colnames(metadata(sce)$experiment_info)[!colnames(metadata(sce)$experiment_info) %in% c("n_cells", "sample_id")], function(x){
      return(0)
    })
    names(catCount) <- colnames(metadata(sce)$experiment_info)[!colnames(metadata(sce)$experiment_info) %in% c("n_cells", "sample_id")]
    exclude <- list()
    for(s in subselection){
      category <- reactiveVals$subselectionMap[[s]]
      catCount[[category]] <- catCount[[category]] + 1 
      if(!is.null(exclude[[category]])){
        exclude[[category]] <- c(exclude[[category]], s)
      }else{
        exclude[[category]] <- s
      }
    }
    for(x in names(catCount)){
      if(x %in% contrastVars & catCount[[x]] == 1){
        showNotification("You want to analyse a condition you subsetted. That is not meaningful. Try again.", type = "error")
        return(NULL)
      }
    }
    reactiveVals$exclusionList <- exclude
    for(cat in names(exclude)){
      sub <- exclude[[cat]]
      print(sprintf("only using %s from the condition %s", sub, cat))
      sce <- filterSCE(sce, get(cat) %in% sub)
    }
  }

  # toggle_inputs()
  ei <- ei(sce)
  nr_samples <- nlevels(sample_ids(sce))


  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR")){
    if(input$normalizeDE == "Yes"){
      normalize <- TRUE
    }else{
      normalize <- FALSE
    }
    for(designCols in input$colsDesign){
      if(length(levels(ei[[designCols]])) < 2){
        msg <- paste("Your design column", designCols, "has less than 2 levels left. Please do not include this in the design matrix.")
        showNotification(msg, type = "error")
        return(NULL)
      }
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
        normalize = toString(isolate(input$normalizeDE)), 
        filter = toString(subselection)
      )
      
      out <- diffcyt::diffcyt(
        d_input = sce,
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
      blockID <- metadata(sce)$experiment_info[[input$blockID_voom]]
    }else{
      blockID <- NULL
    }
    for(designCols in input$colsDesign){
      if(length(levels(ei[[designCols]])) < 2){
        msg <- paste("Your design column", designCols, "has less than 2 levels left. Please do not include this in the design matrix.")
        showNotification(msg, type = "error")
        return(NULL)
      }
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
        normalize = toString(isolate(input$normalizeDE)),
        filter = toString(subselection)
      )
        out <- diffcyt::diffcyt(
          d_input = sce,
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
      blockID <- metadata(sce)$experiment_info[[input$blockID_limma]]
    }else{
      blockID <- NULL
    }
    if(input$trend_limma == "Yes"){
      trend <- TRUE
    }else{
      trend <- FALSE
    }
    
    markersToTest <- input$DEFeaturesIn
    is_marker <- rowData(sce)$marker_class %in% c("type", "state")
    if (input$DEMarkerToTest == "Marker by Class") {
      markersToTest <- (rowData(sce)$marker_class %in% markersToTest)[is_marker] # type and state (but not none)
    }else{
      markersToTest <- rownames(sce)[is_marker] %in% markersToTest
    }
    for(designCols in input$colsDesign){
      if(length(levels(ei[[designCols]])) < 2){
        msg <- paste("Your design column", designCols, "has less than 2 levels left. Please do not include this in the design matrix.")
        showNotification(msg, type = "error")
        return(NULL)
      }
    }
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    contrast <- createCustomContrastMatrix(contrastVars, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples left which is not meaningful. Try again.", type = "error")
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
        trend_method = toString(isolate(input$trend_limma)),
        filter = toString(subselection)
      )
      out <- diffcyt::diffcyt(
        d_input = sce,
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
    is_marker <- rowData(sce)$marker_class %in% c("type", "state")
    if (input$DEMarkerToTest == "Marker by Class") {
      markersToTest <- (rowData(sce)$marker_class %in% markersToTest)[is_marker] # type and state (but not none)
    }else{
      markersToTest <- rownames(sce)[is_marker] %in% markersToTest
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
        features = toString(isolate(input$DEFeaturesIn)),
        filter = toString(subselection)
      )
      
    out <- diffcyt::diffcyt(
      d_input = sce,
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
        normalize = toString(isolate(input$normalizeDE)),
        filter = toString(subselection)
      )
    out <- diffcyt::diffcyt(
      d_input = sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      normalize = normalize
    )
    }
  } else if (input$chosenDAMethod == "sceEMD") {
    markersToTest <- isolate(input$DEFeaturesIn)
    if (isolate(input$DEMarkerToTest) == "Marker by Class") {
      sce <- filterSCE(sce, marker_classes(sce) %in% markersToTest)
    }else{
      sce <- filterSCE(sce, rownames(sce) %in% markersToTest)
    }
    
    binSize <- isolate(input$emdBinwidth)
    reactiveVals$methodsInfo[["sceEMD"]] <- data.frame(
      analysis_type = "Differential States", 
      method = "EMD",
      condition = toString(input$emdCond),
      binSize = binSize,
      nperm = isolate(input$emdNperm),
      clustering = toString(isolate(input$deCluster)),
      features = toString(markersToTest),
      filter = subselection
    )
    
    if (binSize == 0) binSize <- NULL
    
    showNotification(
      ui =
        HTML(
          "<div id='emdProgress'><b>EMD Progress:</b><div>"
        ),
      duration = NULL,
      id = "emdProgressNote"
    )
    
    withCallingHandlers({
      out <-
        sceEMD(
          sce = sce,
          k = isolate(input$deCluster),
          condition = isolate(input$emdCond),
          binSize = binSize,
          nperm = isolate(input$emdNperm)
        )
    },
    message = function(m) {
      shinyjs::html(id = "emdProgress",
                    html = sprintf("<br>%s", HTML(m$message)),
                    add = TRUE)
    })
    
  }
  # toggle_inputs(enable_inputs = TRUE)
  removeNotification("emdProgressNote")
  return(out)
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
  reactiveVals$continue <- TRUE
   methodsDA <- c("edgeR" = "diffcyt-DA-edgeR", "Voom" = "diffcyt-DA-voom", "GLMM" = "diffcyt-DA-GLMM")
   methodsDS <- c("limma" = "diffcyt-DS-limma","LMM" = "diffcyt-DS-LMM", "EMD" = "sceEMD")
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
  req(input$chosenDAMethod, input$chosenDAMethod != "sceEMD")
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
      selected = colsDesign[1],
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
      selected = cols[1],
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
      content = "Fixed effects are variables that we expect will have an effect on the dependent/response variable: they’re what you call explanatory variables in a standard linear regression."),
    pickerInput(
      "colsRandom",
      choices = cols,
      selected = cols[2],
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
      content = "Random effects are usually grouping factors for which we are trying to control. Note that the golden rule is that you generally want your random effect to have at least five levels. For example if you are analysing a condition (A vs. B) on samples (patient1_A, patient1_B) belonging to the same patient (patient1), you can include patient_id as random effect."),
  )
})

output$contrastSelection <- renderUI({
  req(input$chosenDAMethod, input$chosenDAMethod != "sceEMD")
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR", "diffcyt-DS-limma", "diffcyt-DA-voom")) {
    choices <- input$colsDesign
  } else {
    choices <- input$colsFixed
  }
  div(
    selectInput(
      "contrastVars",
      choices = choices,
      selected = choices[1],
      label = span(
        "What condition do you want to analyse?",
        icon("question-circle"),
        id = "deContrastQ"
      ),
      multiple = F
    ),
    bsPopover(
      id = "deContrastQ",
      title = "Contrast Matrix Design",
      content = "Here, you specify the comparison of interest, i.e. the combination of model parameters to test whether they are equal to zero."
    )
  )
})

output$emdInput <- renderUI({
  req(input$chosenDAMethod == "sceEMD")
  sceEI <- ei(reactiveVals$sce)
  condChoices <- which(sapply(sceEI, function(feature) nlevels(as.factor(feature)) == 2))
  if (length(condChoices) == 0){
    showNotification("No condition with exactly two levels found. EMD is not applicable, please select another method.", duration = NULL, type = "error")
    return(NULL)
    }
  condChoices <- names(condChoices)
  list(div(
    selectizeInput(
      "emdCond",
      choices = condChoices,
      label = span(
        "What condition do you want to analyse?",
        icon("question-circle"),
        id = "emdCondQ"
      )
    ),
  bsPopover(
    id = "emdCondQ",
    title = "Condition for EMD analysis",
    content = HTML("Here, you specify the comparison of interest, i.e. the group on which to split the expression distributions.<br><b>Currently only conditions with two levels are supported.</b>")
  )
  ),
  uiOutput("emdNpermInput"),
  div(
    numericInput(
      "emdBinwidth",
      label = span(
        "Bin width for comparing histograms",
        icon("question-circle"),
        id = "emdBinwidthQ"
      ),
      value = 0,
      min = 0,
      max = 1,
      step = .1
    ),
    bsPopover(
      id = "emdBinwidthQ",
      title = "Bin width for comparing histograms",
      content = HTML("You can set a custom binwidth but we recommend to leave this at zero.<br><b>Set this to 0 to compute the binwidth for each marker based on the Freedman-Diaconis rule.</b>")
    )
  ))
})

output$emdNpermInput <- renderUI({
  req(input$emdCond)
  maxPerm <- as.numeric(RcppAlgos::permuteCount(ei(reactiveVals$sce)[[input$emdCond]]))
  div(
    numericInput(
      "emdNperm",
      label = span(
        "Number of permutations for p-value estimation",
        icon("question-circle"),
        id = "emdNpermQ"
      ),
      value = min(100, maxPerm),
      min = 0,
      max = maxPerm,
      step = 10
    ),
    bsPopover(
      id = "emdNpermQ",
      title = "Number of permutations for p-value estimation",
      content = HTML("Note that meaningful results require many permutations. For an unadjusted pvalue smaller than 0.01 at least 100 permutations are necessary.<br><b>This value must not exceed the factorial of the number of samples.</b>")
    )
  )
})

output$deSubselection <- renderUI({
  choices <- isolate(colnames(metadata(reactiveVals$sce)$experiment_info))
  choices <- choices[!choices %in% c("n_cells", "sample_id", "patient_id")]

  map <- unlist(sapply(choices, function(x){
    lvls <- isolate(levels(metadata(reactiveVals$sce)$experiment_info[[x]]))
    return(rep(x, length(lvls)))
  }))

  choices <- unlist(sapply(choices, function(x){
    lvls <- isolate(levels(metadata(reactiveVals$sce)$experiment_info[[x]]))
    return(lvls)
  }))
  names(map) <- choices
  names(choices) <- paste("only", choices)
  reactiveVals$subselectionMap <- map
  div(
    checkboxGroupInput(
      inputId = "deSubselection",
      label = span("Do you want to analyse this condition just on a subset?", icon("question-circle"), id = "subSelectQ"),
      choices = c(choices), 
      inline = T
    ),
    bsPopover(
      id = "subSelectQ",
      title = "Run differential expression on a subset of your data",
      content = "Sometimes it might make sense to compare differential expression just in a subset of your data, e.g. you have two different treatment groups and want to investigate the effect of an activation agent separately. You can do the subselection right at the beginning (Preprocessing) or here.",
      placement = "top"
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
  if (input$chosenDAMethod %in% c("diffcyt-DS-limma", "diffcyt-DS-LMM", "sceEMD")) {
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
    label = "Features to use for differential expression analysis",
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
  out <- reactiveVals$DEruns[[input$deVisMethod]]

  if (input$deVisMethod != "sceEMD")
    out <- rowData(out$res)
  div(
    pickerInput(
      "DEClusterSelection",
      choices = unique(out$cluster_id),
      selected = unique(out$cluster_id),
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
  out <- reactiveVals$DEruns[[input$deVisMethod]]
  
  if (input$deVisMethod != "sceEMD")
    out <- rowData(out$res)
  div(
    pickerInput(
      "DEMarkerSelection",
      choices = as.character(unique(out$marker_id)),
      selected = as.character(unique(out$marker_id)),
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
  out <- reactiveVals$DEruns[[input$deVisMethod]]
  if (input$deVisMethod != "sceEMD")
    out <- rowData(out$res)
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
      value = min(20,nrow(out)),
      min = 1,
      max = nrow(out),
      step = 1
    ), 
  )
})

# choose method and parameter box
output$visDiffExp <- renderUI({
  runs <- names(reactiveVals$DEruns)
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
  if (is.null(reactiveVals$DEruns)){
    reactiveVals$DEruns <- list()
    output$heatmapDEPlot <- renderPlot({
      ggplotObject <- ggplot() + theme_void()
      return(ggplotObject)
    })
  }
  
  if(is.null(reactiveVals$methodsInfo)){
    reactiveVals$methodsInfo <- list()
  }
  

  # call diffcyt function
  out <- call_DE()
  if(!is.null(out)){
    # add method to DAruns
    reactiveVals$DEruns[[DAmethod]] <- out
  
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
  runs <- isolate(reactiveVals$DEruns)
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
    methodsDS <- c("diffcyt-DS-limma","diffcyt-DS-LMM", "sceEMD")
    subselection <- isolate(reactiveVals$methodsInfo[[visMethod]]$filter)
    sce <- isolate(reactiveVals$sce)
    
    if(subselection != "No"){
      exclude <- isolate(reactiveVals$exclusionList)
      for(cat in names(exclude)){
        sub <- exclude[[cat]]
        print(sprintf("only using %s from the condition %s", sub, cat))
        sce <- filterSCE(sce, get(cat) %in% sub)
      }
      #category <- isolate(reactiveVals$subselectionMap[[subselection]])
      #sce <- filterSCE(sce, get(category) == subselection)
    }
    
    if(visMethod %in% methodsDA){
      sub <- filterSCE(sce, cluster_id %in% heatmapSelection, k=deCluster)
      x <- sub
    } else {
      x <- sce[rownames(sce) %in% heatmapSelection, ]
    }
    
    out <- runs[[visMethod]]
    if (visMethod != "sceEMD")
      out <- rowData(out$res)
    out$p_val[is.na(out$p_val)] <- 1
    out$p_adj[is.na(out$p_adj)] <- 1
    
    reactiveVals$diffHeatmapPlot <- plotDiffHeatmapCustom(
      x=x,
      y=out, 
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
    table <- reactiveVals$methodsInfo[[visMethod]]
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
    
    out <- reactiveVals$DEruns[[visMethod]] 
    if (visMethod != "sceEMD")
      out <- diffcyt::topTable(out$res,all=TRUE,format_vals=TRUE)
    out$p_val[is.na(out$p_val)] <- 1
    out$p_adj[is.na(out$p_adj)] <- 1
    reactiveVals$topTable <- data.frame(out)
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
  #if we only have counts -> copy counts to exprs
  if(length(assays(reactiveVals$sce)) == 1){
    showNotification("You have not normalized your data. We will assume that you have given us expression data as input.", type = "warning", duration = 10)
    assays(reactiveVals$sce)$exprs <- assays(reactiveVals$sce)$counts
  }
  
  factors <- names(colData(reactiveVals$sce))[!names(colData(reactiveVals$sce)) %in% c("patient_id", "sample_id")]
  markers <- unique(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
  if("state" %in% markers){
    selected_marker <- "state"
  }else{
    selected_marker <- markers[1]
  }
  choices <- isolate(colnames(metadata(reactiveVals$sce)$experiment_info))
  choices <- choices[!choices %in% c("n_cells", "sample_id", "patient_id")]
  
  choices <- unlist(sapply(choices, function(x){
    lvls <- isolate(levels(metadata(reactiveVals$sce)$experiment_info[[x]]))
    return(lvls)
  }))
  names(choices) <- paste("only", choices)

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
                    selectizeInput(
                      "deBoxSubselect",
                      "Subselection",
                      choices = c("No", choices), 
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
  sce <- isolate(reactiveVals$sce)
  if(input$deBoxSubselect != "No"){
    category <- isolate(reactiveVals$subselectionMap[[input$deBoxSubselect]])
    sce <- filterSCE(sce, get(category) == input$deBoxSubselect)
  }
  reactiveVals$pbExprsPlot <- plotPbExprs(sce, 
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
output$downloadPlotPbExprs <- downloadHandler(
  filename = function(){
    paste0("Pb_Exprs_plot", ".pdf")
  },
  content = function(file){
    ggsave(file, plot =reactiveVals$pbExprsPlot, width=10, height=12)
  }
)

plotDiffHeatmapCustom <- function (x, y, k = NULL, top_n = 20, fdr = 0.05, lfc = 1, all = FALSE, 
                                   sort_by = c("padj", "lfc", "none"), y_cols = list(padj = "p_adj", 
                                                                                     lfc = "logFC", target = "marker_id"), assay = "exprs", 
                                   fun = c("median", "mean", "sum"), normalize = TRUE, col_anno = TRUE, 
                                   row_anno = TRUE, hm_pal = NULL, fdr_pal = c("lightgrey", 
                                                                               "lightgreen"), lfc_pal = c("blue3", "white", "red3")) 
{
  fun <- match.arg(fun)
  sort_by <- match.arg(sort_by)
  args <- as.list(environment())
  defs <- as.list(formals("plotDiffHeatmap")$y_cols[-1])
  miss <- !names(defs) %in% names(args$y_cols)
  if (any(miss)) 
    y_cols <- args$y_cols <- c(args$y_cols, defs[miss])[names(defs)]
  CATALYST:::.check_args_plotDiffHeatmap(args)
  stopifnot(y_cols[[sort_by]] %in% names(y))
  y_cols <- y_cols[y_cols %in% names(y)]
  if (is.null(k)) {
    kids <- levels(y$cluster_id)
    same <- vapply(cluster_codes(x), function(u) identical(levels(u), 
                                                           kids), logical(1))
    if (!any(same)) 
      stop("Couldn't match any clustering", " in input data 'x' with results in 'y'.")
    k <- names(cluster_codes(x))[same][1]
  }
  else {
    k <- CATALYST:::.check_k(x, k)
  }
  x$cluster_id <- cluster_ids(x, k)
  y <- data.frame(y, check.names = FALSE)
  y <- mutate_if(y, is.factor, as.character)
  if (any(rownames(x) %in% unlist(y))) {
    features <- intersect(rownames(x), y[[y_cols$target]])
    if (length(features) == 0) 
      stop("Couldn't match features between", " results 'y' and input data 'x'.")
    i <- y[[y_cols$target]] %in% rownames(x)
    type <- "ds"
  }
  else {
    i <- TRUE
    type <- "da"
  }
  y <- dplyr::rename(y, target = y_cols$target, padj = y_cols$padj, 
                     lfc = y_cols$lfc)
  i <- i & !is.na(y$padj) & y$cluster_id %in% levels(x$cluster_id)
  if (!all) {
    i <- i & y$padj < fdr
    if (!is.null(y$lfc)) 
      i <- i & abs(y$lfc) > lfc
  }
  y <- y[i, , drop = FALSE]
  if (nrow(y) == 0) 
    stop("No results remaining;", " perhaps 'x' or 'y' has been filtered,", 
         " or features couldn't be matched.")
  if (sort_by != "none") {
    o <- order(abs(y[[sort_by]]), decreasing = (sort_by == 
                                                  "lfc"))
    y <- y[o, , drop = FALSE]
  }
  if (top_n > nrow(y)) 
    top_n <- nrow(y)
  top <- y[seq_len(top_n), ]
  if (!isFALSE(col_anno)) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, 
                              "column")
  }
  else top_anno <- NULL
  if (is.null(hm_pal)) 
    hm_pal <- rev(RColorBrewer::brewer.pal(11, ifelse(type == "ds", "RdYlBu", 
                                                      "RdBu")))
  if (row_anno) {
    s <- factor(ifelse(top$padj < fdr, "yes", "no"), levels = c("no", 
                                                                "yes"))
    if (!is.null(top$lfc)) {
      lfc_lims <- range(top$lfc, na.rm = TRUE)
      if (all(lfc_lims > 0)) {
        lfc_brks <- c(0, lfc_lims[2])
        lfc_pal <- lfc_pal[-1]
      }
      else if (all(lfc_lims < 0)) {
        lfc_brks <- c(lfc_lims[1], 0)
        lfc_pal <- lfc_pal[-3]
      }
      else lfc_brks <- c(lfc_lims[1], 0, lfc_lims[2])
      lfc_anno <- top$lfc
      anno_cols <- list(logFC = circlize::colorRamp2(lfc_brks, lfc_pal))
    }
    else {
      lfc_anno <- NULL
      anno_cols <- list()
    }
    names(fdr_pal) <- levels(s)
    anno_cols$significant <- fdr_pal
    right_anno <- ComplexHeatmap::rowAnnotation(logFC = lfc_anno, significant = s, 
                                foo = ComplexHeatmap::row_anno_text(scales::scientific(top$padj, 2), gp = gpar(fontsize = 8)), 
                                col = anno_cols, gp = gpar(col = "white"), show_annotation_name = FALSE, 
                                simple_anno_size = unit(4, "mm"))
  }
  else right_anno <- NULL
  switch(type, da = {
    ns <- table(x$cluster_id, x$sample_id)
    fq <- prop.table(ns, 2)
    fq <- fq[top$cluster_id, ]
    y <- as.matrix(unclass(fq))
    if (normalize) y <- CATALYST:::.z_normalize(asin(sqrt(y)))
    ComplexHeatmap::Heatmap(matrix = y, name = paste0("normalized\n"[normalize], 
                                      "frequency"), col = hm_pal, na_col = "lightgrey", 
            rect_gp = gpar(col = "white"), cluster_rows = FALSE, 
            cluster_columns = FALSE, row_names_side = "left", 
            top_annotation = top_anno, right_annotation = right_anno)
  }, ds = {
    y <- assay(x, assay)
    cs <- CATALYST:::.split_cells(x, c("cluster_id", "sample_id"))
    z <- t(mapply(function(k, g) vapply(cs[[k]], function(cs) {
      if (length(cs) == 0) return(NA)
      get(fun)(y[g, cs, drop = FALSE])
    }, numeric(1)), k = top$cluster_id, g = top$target))
    rownames(z) <- sprintf("%s(%s)", top$target, top$cluster_id)
    if (normalize) z <- CATALYST:::.z_normalize(z)
    ComplexHeatmap::Heatmap(matrix = z, name = paste0("z-normalized\n"[normalize], 
                                                      "expression"), col = hm_pal, cluster_rows = FALSE, 
                            cluster_columns = FALSE, top_annotation = top_anno, 
                            row_names_side = "left", rect_gp = gpar(col = "white"), 
                            right_annotation = right_anno, heatmap_legend_param = list(title_gp = gpar(fontsize = 10, 
                                                                                                       fontface = "bold", lineheight = 0.8)))
  })
}
