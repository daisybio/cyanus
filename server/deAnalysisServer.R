library(diffcyt)

resetDE <- function(){
  reactiveVals$methodType <- NULL
  reactiveVals$exclusionList <- NULL
  reactiveVals$subselectionMap <- NULL
  reactiveVals$methodsInfo <- NULL
  reactiveVals$eff_r <- NULL
  reactiveVals$DEruns <- NULL
  reactiveVals$pbExprsPlot <- NULL
  reactiveVals$diffHeatmapPlot <- NULL
  reactiveVals$topTable <- NULL
}

plotbox_height <- "48em"
methods_height <- "40em"

methodsDA <- c("edgeR" = "diffcyt-DA-edgeR", "Voom" = "diffcyt-DA-voom", "GLMM" = "diffcyt-DA-GLMM")
methodsDS <- c("limma" = "diffcyt-DS-limma","LMM" = "diffcyt-DS-LMM", "EMD" = "sceEMD", "CytoGLMM" = "CytoGLMM", "CytoGLM" = "CytoGLM", "Wilcoxon rank-sum test" = "wilcoxon_median", "Student's t-test" = "t_test")


# Main function: 
# checks which methods is selected and executes the runDS function accordingly
 
call_DE <- function() {
  
  # read/initialize input   ------------------------------------------

  method <- isolate(input$chosenDAMethod)
  methodType <- isolate(reactiveVals$methodType)
  sce <- isolate(reactiveVals$sce)
  ei <- ei(sce)
  condition <- isolate(input$conditionIn)
  groupCol <- isolate(input$groupCol)
  if (groupCol == "") groupCol <- NULL
  
  clustering_to_use <- isolate(input$deCluster)
  nr_samples <- nlevels(sample_ids(sce))
  
  addTerms <- isolate(input$addTerms)
  
  normalize <- NULL
  featuresIn <- NULL
  includeWeights <- NULL
  sceEMD_nperm <- NULL
  sceEMD_binsize <- NULL
  CytoGLM_num_boot <- NULL
  trend_edgeR <- NULL
  trend_limma <- NULL
  blockID <- NULL
  
  #check if subselection is valid, set subselection:
  if (is.null(input$deSubselection)) {
    subselection <- "No"
  } else{
    subselection <- isolate(input$deSubselection)
  }
  if (subselection != "No") {
    if (is.null(reactiveVals$exclusionList)) {
      reactiveVals$exclusionList <- list()
    }
    showNotification(
      ui =
        HTML("<div id='subselection'><b>Subselecting...</b><div>"),
      duration = NULL,
      id = "subselectionNote",
      type = "error"
    )
    withCallingHandlers({
      returned <-
        doConditionSubselection(
          sce,
          subselection,
          isolate(reactiveVals$subselectionMap),
          condition,
          reactiveVals$exclusionList,
          method
        )
    },
    message = function(m) {
      shinyjs::html(id = "subselection",
                    html = sprintf("<br>%s", HTML(m$message)),
                    add = TRUE)
    })
    removeNotification("subselection")
    if(is.null(returned)) return(NULL)
    sce <- returned[["sce"]]
    reactiveVals$exclusionList <- returned[["exclusionList"]]
    remove(returned)
  }
  
  # error handling, make parameters ------------------------------------------
  
  parameters <- list()
  
  if (method %in% c("diffcyt-DA-edgeR", "diffcyt-DA-voom")) {
    design_matrix <- c(groupCol, condition, addTerms)
    parameters[[method]] <- prepDiffExp(
      sce = sce,
      contrastVars = condition,
      colsDesign = design_matrix,
      method = method
    )
    for (designCols in design_matrix) {
      if (length(levels(ei[[designCols]])) < 2) {
        msg <-
          paste(
            "Your design column",
            designCols,
            "has less than 2 levels left. Please do not include this in the design matrix."
          )
        showNotification(msg, type = "error")
        return(NULL)
      }
    }
    if (ncol(parameters[[method]][["design"]]) >= nr_samples) {
      showNotification(
        "You selected more conditions than there are samples which is not meaningful. Try again.",
        type = "error"
      )
      return(NULL)
    }
    if (method == "diffcyt-DA-voom" &&
        input$blockID_voom %in% design_matrix) {
      showNotification(
        "Please don't put your blocking variable in the design matrix. See our tooltip for more information",
        type = "error"
      )
      return(NULL)
    }
    # if (method == "diffcyt-DS-limma" &&
    #     input$blockID_limma %in% design_matrix) {
    #   showNotification(
    #     "Please don't put your blocking variable in the design matrix. See our tooltip for more information",
    #     type = "error"
    #   )
    #   return(NULL)
    # }
    
  } else if (method %in% c("diffcyt-DA-GLMM")) {
    fixed_effects = c(condition, addTerms)
    random_effects = groupCol
    parameters[[method]] <- prepDiffExp(
      sce = sce,
      contrastVars = condition,
      colsFixed = fixed_effects,
      colsRandom = random_effects,
      method = method
    )
    if (nrow(parameters[[method]][["contrast"]]) >= nr_samples) {
      showNotification(
        "You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.",
        type = "error"
      )
      return(NULL)
    }
  }
  
  # set method dependent parameters ------------------------------------------
  
  
  if (method %in% methodsDA) {
    normalize <- ifelse(input$normalizeDE == "Yes", TRUE, FALSE)
  } 
  if (method == "diffcyt-DA-edgeR") {
    trend_edgeR <- isolate(input$edgeR_trendMethod)
  }
  if (method == "diffcyt-DA-voom" && input$blockID_voom != "") {
    blockID <- metadata(sce)$experiment_info[[input$blockID_voom]]
  } else if (method == "diffcyt-DS-limma" &&
             input$blockID_limma != "") {
    blockID <- metadata(sce)$experiment_info[[input$blockID_limma]]
  }
  if (method %in% methodsDS) {
    featuresIn <- isolate(input$DEFeaturesIn)
    if (method == "diffcyt-DS-limma") {
      trend_limma <- ifelse(input$trend_limma == "Yes", TRUE, FALSE)
    }
    if (method %in% c("diffcyt-DS-limma", "diffcyt-DS-LMM")) {
      includeWeights <- isolate(input$weightsSelection)
      includeWeights <- ifelse(includeWeights == "Yes", TRUE, FALSE)
    } else if (input$chosenDAMethod == "sceEMD") {
      sceEMD_binsize <- isolate(input$emdBinwidth)
      if (sceEMD_binsize == 0)
        sceEMD_binsize <- NULL
      emdNperm <- isolate(input$emdNperm)
    } else if (input$chosenDAMethod == "CytoGLM") {
      CytoGLM_num_boot <- isolate(input$CytoGLM_num_boot)
    }else if(method == "CytoGLMM" & is.null(groupCol)){
      shiny::showNotification("Results of CytoGLMM are not meaningful when no grouping variable like patient ID is selected.", 
                              type = "warning", duration = 10)
    }
  }
  
  # make info data table ------------------------------------------
  
  reactiveVals$methodsInfo[[method]] <- data.frame(
    analysis_type = isolate(input$da_ds),
    method = method,
    condition = condition,
    grouping_columns = toString(groupCol),
    additional_fixed_effects = toString(addTerms),
    clustering = toString(clustering_to_use),
    normalize = toString(normalize),
    filter = toString(subselection),
    features = toString(featuresIn),
    trend_edgeR = toString(trend_edgeR),
    trend_limma = toString(trend_limma),
    block_id = toString(blockID),
    include_weights = toString(includeWeights),
    sceEMD_binsize = toString(sceEMD_binsize),
    sceEMD_nperm = toString(sceEMD_nperm),
    CytoGLM_num_boot = toString(CytoGLM_num_boot)
  )
  
  # MAIN: run method ------------------------------------------

    showNotification(
      ui =
        HTML("<div id='deProgress'><b>DE Progress:</b><div>"),
      duration = NULL,
      id = "progressNote"
    )
    
    withCallingHandlers({
      #extra args: sceEMD_condition, binSize, nperm
      if (methodType == "DA") {
        results <- runDA(
          sce = sce, 
          da_methods = method,
          parameters = parameters,
          clustering_to_use = clustering_to_use, 
          normalize = normalize, 
          trend_edgeR = trend_edgeR, 
          blockID = blockID
        )
        out <- results[[method]]
        reactiveVals$eff_r[[method]] <- NULL
      } else {
        #extra args: parameters, blockID, trend_limma, markersToTest, includeWeights
        results <- runDS(
          sce = sce,
          ds_methods = method,
          clustering_to_use = clustering_to_use,
          contrast_vars = condition,
          markers_to_test = featuresIn,
          parameters = parameters,
          blockID = blockID,
          trend_limma = trend_limma,
          include_weights = includeWeights,
          design_matrix_vars = c(condition, addTerms, groupCol), 
          fixed_effects = c(condition, addTerms), 
          random_effects = groupCol,
          sceEMD_nperm = emdNperm, 
          sceEMD_binsize = sceEMD_binsize,
          cytoGLMM_num_boot = CytoGLM_num_boot,
          time_methods = FALSE,
          parallel = FALSE
        )
        out <- results[[method]]
        #the effect sizes do not have to be computed multiple times
        reactiveVals$eff_r[[method]] <- findEffectSize(sce, condition, groupCol, clustering_to_use)
      }
    },
    message = function(m) {
      shinyjs::html(id = "deProgress",
                    html = sprintf("<br>%s", HTML(m$message)),
                    add = TRUE)
    })
    
  
  removeNotification("progressNote")
  return(out)
}

findEffectSize <- function(sce, condition, groupCol, clustering_to_use){
  #the effect sizes do not have to be computed multiple times
  found_effect_size <- FALSE
  if(!is.null(reactiveVals$eff_r)){
    for(method in names(reactiveVals$eff_r)){
      eff_size <- reactiveVals$eff_r[[method]]
      if(is.null(groupCol) & !("grouped" %in% eff_size$overall_group)){
        found_effect_size <- TRUE
        return(eff_size)
      }else if(!is.null(groupCol) & "grouped" %in% eff_size$overall_group){
        if(reactiveVals$methodsInfo[[method]]$grouping_columns == groupCol){
          found_effect_size <- TRUE
          return(eff_size)
        }
      }
    }
  }
  if(is.null(reactiveVals$eff_r) | !found_effect_size){
    return(effectSize(sce = sce,
                      condition = condition,
                      group = groupCol, k=clustering_to_use, 
                      use_assay="exprs", use_mean=FALSE)
           )
  }
  
}

# Renderer ----

output$selectionBoxDE <- renderUI({
  shinydashboard::box(
    column(
    radioButtons(
      inputId = "da_ds",
      label = span(
        "What type of testing do you want to perform?",
        icon("question-circle"),
        id = "da_dsQ"
      ),
      choices = c("Differential Cluster Abundance", "Differential Marker Expression"),
      inline = T
    ),
    bsPopover(
      id = "da_dsQ",
      title = "Analysis type",
      content = HTML(
        "Before doing this, you should have done clustering, preferrably by type. <br> <b>Differential Cluster Abundance:</b> Differential analysis of cell population abundance regarding the clusters. Compares the proportions of cell types across experimental condition and aims to highlight populations that are present at different ratios. <br> <b>Differential Marker Expression:</b> Differential analysis of the marker expression in each cell population (i.e. cluster or overall)."
      )
    ),
    uiOutput("deMethodSelection"),
    uiOutput("conditionSelection"),
    uiOutput("groupSelection"),
    uiOutput("additionalTermsSelection"),
    # uiOutput("contrastSelection"),
    # uiOutput("modelSelection"),
    uiOutput("emdInput"),
    uiOutput("CytoGLM_num_boot"),
    uiOutput("deSubselection"),
    width = 6
  ),
  column(
    uiOutput("clusterSelection"),
    uiOutput("markerToTestSelection"),
    uiOutput("extraFeatures"),
    uiOutput("normalizeSelection"),
    uiOutput("weightSelection"),
    width = 6),
  div(
    bsButton(
      "diffExpButton",
      "Start Analysis",
      icon = icon("tools"),
      style = "success"
    ),
    style = "float: right; bottom:5px"
  ),
  title = "Choose Method and Parameters",
  width = 12,
  height = methods_height
  )
})

# displays available methods and selection of DA or DS
output$deMethodSelection <- renderUI({
   if(input$da_ds == "Differential Cluster Abundance"){
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

output$conditionSelection <- renderUI({
  sceEI <- CATALYST::ei(reactiveVals$sce)
  condChoices <- which(sapply(sceEI, function(feature) nlevels(as.factor(feature)) == 2))
  if (length(condChoices) == 0) {
    showNotification("No condition with exactly two levels found. Unfortunately we currently only support comparisons between two conditions. You might want to subset your data.", duration = NULL, type = "error")
    return(NULL)
  }
  condChoices <- names(condChoices)
  list(div(
    selectizeInput(
      "conditionIn",
      choices = condChoices,
      label = span(
        "What condition do you want to analyse?",
        icon("question-circle"),
        id = "conditionInQ"
      )
    ),
    bsPopover(
      id = "conditionInQ",
      title = "Condition for DE analysis",
      content = HTML("Here, you specify the comparison of interest.<br><b>Currently only conditions with two levels are supported.</b>")
    )))
})

output$groupSelection <- renderUI({
  all_methods <- c(methodsDA, methodsDS)
  req(input$conditionIn, input$chosenDAMethod %in% all_methods[all_methods != 'sceEMD'])
  sceEI <- data.table::as.data.table(CATALYST::ei(reactiveVals$sce))
  groupCol <- names(sceEI)[!names(sceEI) %in% c("n_cells", "sample_id")]
  groupCol <- groupCol[sapply(groupCol, function(x) sceEI[, .(e2 = data.table::uniqueN(get(input$conditionIn)) == 2),, by=get(x)][, all(e2)])]
  names(groupCol) <- groupCol
  groupCol <- c('unpaired samples' = '', groupCol)
  div(
    pickerInput(
      "groupCol",
      choices = groupCol,
      label = span(
        "Do you have paired samples? Which column identifies the group e.g. patient_id?",
        icon("question-circle"),
        id = "groupColQ"
      ),
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = FALSE
    ),
    bsPopover(
      id = "groupColQ",
      title = "TODO",
      content = "TODO: explain paired data"
    )
  )
})
output$additionalTermsSelection <- renderUI({
  req(input$chosenDAMethod, startsWith(input$chosenDAMethod, 'diffcyt')) # this means this is a linear model and additional terms are allowed
  addTerms <- names(ei(reactiveVals$sce))
  addTerms <- addTerms[!addTerms %in% c("n_cells", "sample_id", input$conditionIn, input$groupCol)]
  div(
    pickerInput(
      "addTerms",
      choices = addTerms,
      label = span(
        "Additional fixed terms to include in the Model",
        icon("question-circle"),
        id = "addTermsQ"
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
      id = "addTermsQ",
      title = "Additional terms to include in the Model",
      content = "Since you have specified a method using a linear model, you can include additional terms that might have an effect on expression other than the condition, such as batches."
    )
  )
})

# selection of weights for analysis
output$weightSelection <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod %in% c("diffcyt-DS-LMM", "diffcyt-DS-limma")){
    div(
      radioButtons(
        inputId = "weightsSelection",
        label = span("Do you want to include precision weights (cell counts) within the model?", icon("question-circle"), id = "weightSelectQ"),
        choices = c("Yes", "No"), 
        inline = T
      ),
      bsPopover(
        id = "weightSelectQ",
        title = "Whether to include precision weights within each model (across samples).",
        content = "These represent the relative uncertainty in calculating each median value. The cell counts of each sample are incorporated in the analysis.",
        placement = "top"
      )
    )
  }
})

# box with cluster populations you want to compare
output$clusterSelection <- renderUI({
  req(reactiveVals$methodType)
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
        content = "Whether to include optional normalization factors to adjust for composition effects. Only relevant for Differential Cluster Abundance methods.",
        placement = "top"
      )
    )
  }
})

output$emdInput <- renderUI({
  req(input$chosenDAMethod == "sceEMD")
  sceEI <- ei(reactiveVals$sce)
  # condChoices <- which(sapply(sceEI, function(feature) nlevels(as.factor(feature)) == 2))
  # if (length(condChoices) == 0){
  #   showNotification("No condition with exactly two levels found. EMD is not applicable, please select another method.", duration = NULL, type = "error")
  #   return(NULL)
  # }
  # condChoices <- names(condChoices)
  list(
  #   div(
  #   selectizeInput(
  #     "emdCond",
  #     choices = condChoices,
  #     label = span(
  #       "What condition do you want to analyse?",
  #       icon("question-circle"),
  #       id = "emdCondQ"
  #     )
  #   ),
  # bsPopover(
  #   id = "emdCondQ",
  #   title = "Condition for EMD analysis",
  #   content = HTML("Here, you specify the comparison of interest, i.e. the group on which to split the expression distributions.<br><b>Currently only conditions with two levels are supported.</b>")
  # )
  # ),
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
  req(input$conditionIn)
  maxPerm <- as.numeric(RcppAlgos::permuteCount(ei(reactiveVals$sce)[[input$conditionIn]]))
  div(
    numericInput(
      "emdNperm",
      label = span(
        "Number of permutations for p-value estimation",
        icon("question-circle"),
        id = "emdNpermQ"
      ),
      value = min(500, maxPerm),
      min = 0,
      max = maxPerm,
      step = 100
    ),
    bsPopover(
      id = "emdNpermQ",
      title = "Number of permutations for p-value estimation",
      content = HTML("Note that meaningful results require many permutations. E.g. for an unadjusted pvalue smaller than 0.01 at least 100 permutations are necessary.<br><b>This value must not exceed the factorial of the number of samples.</b>")
    )
  )
})

output$CytoGLM_num_boot <- renderUI({
  req(input$chosenDAMethod == "CytoGLM")
  div(
    numericInput(
      "cytoNBoot",
      label = span(
        "Number of bootstrap samples",
        icon("question-circle"),
        id = "cytoNBootQ"
      ),
      value = 1000,
      min = 0,
      max = 10000,
      step = 100
    ),
    bsPopover(
      id = "cytoNBootQ",
      title = "Number of bootstrap samples",
      content = HTML('CytoGLM uses bootstrapping with replacement to preserve the cluster structure in donors. For more information refer to <a href="https://doi.org/10.1186/s12859-021-04067-x" target="_blank">Seiler et al.</a>')
    )
  )
})

output$deSubselection <- renderUI({
  choices <- colnames(metadata(reactiveVals$sce)$experiment_info)
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
      title = "Run differential analysis on a subset of your data",
      content = "Sometimes it might make sense to compare differential expression/abundance just in a subset of your data, e.g. you have two different treatment groups and want to investigate the effect of an activation agent separately. You can do the subselection right at the beginning (Preprocessing) or here.",
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
        choices = c("none", "locfit", "movingave", "loess", "locfit.mixed"),
        multiple = F
      ),
      bsPopover(
        id = "trendMethodQ",
        title = "Method for estimating dispersion trend",
        content = "EdgeR specific parameter for estimating the dispersion trend. See more in their estimateDisp function documentation. By default, we use the option none to calculate dispersions, since the dispersion-mean relationship typically does not resemble RNA-sequencing data.",
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
      content = "Here, you specify the comparison of interest. The p-values will be calculated on the basis of this variable, i.e. it will be tested whether the coefficient of this parameter in the model is equal to zero. "
    )
  )
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
      content = "The selected columns will be included in the design matrix specifying the models to be fitted. A design matrix includes all variables that are relevant/interesting for your analysis because the linear model will be built using these variables. That means that you MUST include the condition you want to analyze here."
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
      content = "Fixed effects are variables that we expect will have an effect on the dependent/response variable: they’re what you call explanatory variables in a standard linear regression. You MUST include the condition variable here. "),
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

# if diffcyt should be exectued on selected markers (markers_to_test)
output$markerToTestSelection <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod %in% c("diffcyt-DS-limma", "diffcyt-DS-LMM", "sceEMD", "CytoGLMM", "CytoGLM", "wilcoxon_median", "t_test")) {
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
    is_marker <- SummarizedExperiment::rowData(reactiveVals$sce)$marker_class %in% c("type", "state")
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
    label = "Features to use for differential analysis",
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
    title = "Visualize Differential Analysis Results",
    height = plotbox_height,
    width = 12,
  )
})

# heatmap can be visualized for subset of clusters (DA)
output$visClusterSelection <- renderUI({
  out <- reactiveVals$DEruns[[input$deVisMethod]]

  #if (input$deVisMethod != "sceEMD")
  #  out <- rowData(out$res)
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
  req(!is.null(reactiveVals$DEruns[[input$deVisMethod]]))
  out <- reactiveVals$DEruns[[input$deVisMethod]]
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
    # TODO: add the possibilty to alter Fold Change
    # numericInput(
    #   "lfcThreshold",
    #   label = "Log2FC Threshold",
    #   value = 1,
    #   step = 0.5
    # ),
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

### BOXPLOT FUNCTIONS

output$deBoxPlots <- renderUI({
  uiOutput("deExprsCluster")
})

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
                      multiple = F
                      # selected = "meta9" #TODO: take first 
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
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
  #shinyjs::show("selectionBoxDE")
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
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot =reactiveVals$pbExprsPlot, width=10, height=12)
    waiter_hide(id="app")
  }
)


# Observer-----

observe({
  if (reactiveVals$current_tab==6){
    shinyjs::hide("DEVisualization")
    shinyjs::hide("selectionBoxDE")
    req(reactiveVals$pbExprsPlot)
    shinyjs::show("selectionBoxDE")
    req(reactiveVals$DEruns)
    shinyjs::show("DEVisualization")
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
  waiter_show(
    id = "app",
    html = tagList(
      spinner$logo,
      HTML("<br>DE Analysis in Progress...<br>Please be patient")
    ),
    color = spinner$color
  )
  out <- call_DE()
  waiter_hide(id = "app")
  if(!is.null(out)){
    # add method to DAruns
    reactiveVals$DEruns[[DAmethod]] <- out
  
    # other method can be performed
    updateButton(session,
               "diffExpButton",
               label = " Start Analysis",
               disabled = FALSE)
  
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
  # lfcThreshold <- isolate(input$lfcThreshold)
  heatmapSortBy <- isolate(input$heatmapSortBy)
  heatmapNormalize <- isolate(input$heatmapNormalize)
  deCluster <- isolate(input$deCluster)
  runs <- isolate(reactiveVals$DEruns)
  topN <- isolate(input$topN)
  shinyjs::disable("previousTab")
  shinyjs::disable("nextTab")
  shinyjs::disable("diffExpButton")

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
    subselection <- isolate(reactiveVals$methodsInfo[[visMethod]]$filter)
    sce <- isolate(reactiveVals$sce)
    if(subselection != "No"){
      exclude <- isolate(reactiveVals$exclusionList[[visMethod]])
      for(cat in names(exclude)){
        sub <- exclude[[cat]]
        print(sprintf("only using %s from the condition %s", sub, cat))
        sce <- filterSCE(sce, get(cat) %in% sub)
      }
    }
    
    if(visMethod %in% methodsDA){
      sub <- filterSCE(sce, cluster_id %in% heatmapSelection, k=deCluster)
      x <- sub
    } else {
      x <- sce[rownames(sce) %in% heatmapSelection, ]
    }
    
    out <- runs[[visMethod]]
    eff_r <- isolate(reactiveVals$eff_r[[visMethod]])
    #if (visMethod != "sceEMD")
    #  out <- rowData(out$res)
    #out$p_val[is.na(out$p_val)] <- 1
    #out$p_adj[is.na(out$p_adj)] <- 1
    reactiveVals$diffHeatmapPlot <- plotDiffHeatmapCustom(
      x=x,
      y=out, 
      fdr=as.numeric(fdrThreshold), 
      # lfc=as.numeric(lfcThreshold), 
      sort_by = heatmapSortBy, 
      normalize=as.logical(heatmapNormalize ),
      all = TRUE,
      top_n = as.numeric(topN),
      eff_r = eff_r
    )
    reactiveVals$diffHeatmapPlot
  })
  
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
  shinyjs::enable("diffExpButton")
  
  # ui for download button
  output$heatmapPlotDownload <- renderUI({
    req(reactiveVals$diffHeatmapPlot)
    downloadButton("downloadPlotDiffHeatmap", "Download Plot")
  })
  
  # function for downloading heatmap
  output$downloadPlotDiffHeatmap <- downloadHandler(
    filename = "DE_Heatmap.pdf", 
    content = function(file){
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      pdf(file, width = 12, height = 8)
      draw(reactiveVals$diffHeatmapPlot)
      dev.off()
      waiter_hide(id="app")
    }
  )
  
  ## Info button
  output$infoDE <- renderUI({
    table <- reactiveVals$methodsInfo[[visMethod]]
    columnsToKeep <- sapply(table, function(x){
      ifelse(x[1] == "", FALSE, TRUE)
    })
    table <- table[, columnsToKeep]
    value <- renderTable(
      checkNullTable(table),
      caption.placement = "top"
    )
    
    dropdownButton(
      value,
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
    reactiveVals$topTable <- data.frame(out)
    DT::datatable(reactiveVals$topTable, rownames = FALSE, 
                  options = list(pageLength = 10, searching = FALSE, 
                                 columnDefs = list(list( targets = c(1,2), 
                                                         render = JS("function(data, type, row, meta) {","return data === null ? 'NA' : data;","}")))))
  })
  
  output$downloadTopTable <- downloadHandler(
    filename = "Differential_Expression_Results.csv",
    content = function(file) {
      waiter_show(id = "app",html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      write.csv(reactiveVals$topTable, file, row.names = FALSE)
      waiter_hide(id="app")
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
    shinyjs::disable("diffExpButton")
  } else {
    shinyjs::enable("diffExpButton")
  }
  }, ignoreNULL=FALSE, ignoreInit=TRUE)
