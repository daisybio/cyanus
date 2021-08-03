library(diffcyt)

resetDE <- function(){
  shinyjs::reset("selectionBoxDE")
  shinyjs::reset("da_ds")
  shinyjs::reset("conditionIn")
  #reactiveVals$methodType <- NULL
  reactiveVals$exclusionList <- NULL
  reactiveVals$subselectionMap <- NULL
  reactiveVals$methodsInfo <- NULL
  reactiveVals$eff_r <- NULL
  reactiveVals$DEruns <- NULL
  reactiveVals$pbExprsPlot <- NULL
  reactiveVals$diffHeatmapPlot <- NULL
  reactiveVals$topTable <- NULL
  reactiveVals$visExpClicked <- NULL
}

plotbox_height <- "48em"
methods_height <- "48em"

methodsDA <- c("edgeR" = "diffcyt-DA-edgeR", "Voom" = "diffcyt-DA-voom", "GLMM" = "diffcyt-DA-GLMM")
methodsDS <- c("limma" = "diffcyt-DS-limma","LMM" = "diffcyt-DS-LMM", "CyEMD" = "CyEMD", "CytoGLMM" = "CytoGLMM", "CytoGLM" = "CytoGLM", "Wilcoxon rank-sum test" = "wilcoxon_median", "Student's t-test" = "t_test")


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
  cyEMD_nperm <- NULL
  cyEMD_binsize <- NULL
  CytoGLM_num_boot <- NULL
  trend_edgeR <- NULL
  trend_limma <- NULL
  blockID <- NULL
  
  #check if downsampling should be performed
  downsampling <- isolate(input$downsampling_Yes_No_DE)
  if(downsampling == "Yes"){
    downsamplingNumber <- isolate(input$downsamplingNumberDE)
    downsamplingPerSample <- isolate(input$downsampling_per_sampleDE)
    if(downsamplingPerSample == "Yes"){
      downsamplingPerSample <- TRUE
    }else{
      downsamplingPerSample <- FALSE
    }
    downsamplingSeed <- isolate(input$downsamplingSeedDE)
    showNotification("Subsampling the data...", type = "warning")
    sce <- performDownsampling(sce, downsamplingPerSample, downsamplingNumber, downsamplingSeed)
    if(is.null(sce)){
      return(NULL)
    }
  }
  
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
    } else if (input$chosenDAMethod == "CyEMD") {
      cyEMD_binsize <- isolate(input$emdBinwidth)
      if (cyEMD_binsize == 0)
        cyEMD_binsize <- NULL
      emdNperm <- isolate(input$emdNperm)
    } else if (input$chosenDAMethod == "CytoGLM") {
      CytoGLM_num_boot <- isolate(input$CytoGLM_num_boot)
    }else if(method == "CytoGLMM" & is.null(groupCol)){
      shiny::showNotification("Results of CytoGLMM are not meaningful when no grouping variable like patient ID is selected.", 
                              type = "warning", duration = NULL)
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
    cyEMD_binsize = toString(cyEMD_binsize),
    cyEMD_nperm = toString(cyEMD_nperm),
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
      #extra args: cyEMD_condition, binSize, nperm
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
          cyEMD_nperm = emdNperm, 
          cyEMD_binsize = cyEMD_binsize,
          cytoGLMM_num_boot = CytoGLM_num_boot,
          time_methods = FALSE,
          parallel = FALSE
        )
        out <- results[[method]]
        reactiveVals$eff_r[[method]] <- effectSize(sce = sce,
                                                   condition = condition,
                                                   group = groupCol, k=clustering_to_use, 
                                                   use_assay="exprs", use_mean=FALSE)
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
      choices = c("Differential Marker Expression", "Differential Cluster Abundance"),
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
    uiOutput("emdInput"),
    uiOutput("CytoGLM_num_boot"),
    uiOutput("deSubselection"),
    uiOutput("downsamplingDE"),
    width = 6
  ),
  column(
    uiOutput("clusterSelection"),
    uiOutput("markerToTestSelection"),
    uiOutput("extraFeatures"),
    uiOutput("normalizeSelection"),
    uiOutput("weightSelection"),
    div(
      bsButton(
        "diffExpButton",
        "Start Analysis",
        icon = icon("tools"),
        style = "success"
      ),
      style = "float: right; bottom:5px"
    ),
    width = 6),
  title = "Choose Method and Parameters",
  width = 12,
  height = methods_height,
  id = "selectionBoxDE_box"
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
  req(input$conditionIn, input$chosenDAMethod %in% all_methods[all_methods != 'CyEMD'])
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
      title = "Paired study design",
      content = "In a paired study design, every sample has another corresponding sample. These samples are in different conditions. Here, you can select the column which identifies what samples are paired. Often this is something like the PatientID."
    )
  )
})
output$additionalTermsSelection <- renderUI({
  req(input$chosenDAMethod, (startsWith(input$chosenDAMethod, 'diffcyt')||startsWith(input$chosenDAMethod, 'CytoGLM')) ) # this means this is a linear model and additional terms are allowed
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
  req(input$chosenDAMethod == "CyEMD")
  sceEI <- ei(reactiveVals$sce)
  list(
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
      content = HTML('CytoGLM uses bootstrapping with replacement to preserve the cluster structure in donors. For more information refer to <a href="https://doi.org/10.1186/s12859-021-04067-x" target="_blank">Seiler et al.</a><br><b>Setting this number very high has a great influence on runtime.</b>')
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

output$downsamplingDE <- renderUI({
  req(reactiveVals$sce)
  smallest_n <- min(CATALYST::ei(reactiveVals$sce)$n_cells)
  sum_n <- sum(CATALYST::ei(reactiveVals$sce)$n_cells)
  div(
    column(
      radioButtons(
        "downsampling_Yes_No_DE",
        label = span("Do you want to perform downsampling?", icon("question-circle"), id="dsDEPopover"),
        choices = c("Yes", "No"),
        selected = "No",
        inline = TRUE
      ),
      bsPopover(
        id="dsDEPopover",
        title = "Downsample your data",
        content = "If you have a big dataset and do not want to wait too long for your analyses, you can perform a downsampling on your dataset. If you choose to downsample per sample, the number of cells you specify will be randomly picked from each sample. Otherwise, the number you specify will be divided by the number of samples and this number will be randomly picked from each sample. If the number is bigger than the sample size, all cells from this sample will be taken.",
        placement = "top"
      ),
      width = 3,
      style='padding:0px;'
    ),
    column(
      numericInput(
        "downsamplingNumberDE",
        label=sprintf("How many cells? (Lowest #cells/sample: %s)", smallest_n),
        value=ifelse(smallest_n > 20000, 20000, smallest_n),
        min=1000,
        max=sum_n,
        step=1000
        ),
      width = 3
      ),
    column(
      radioButtons(
        "downsampling_per_sampleDE",
        label = "Per sample?",
        choices = c("Yes", "No"),
        inline = TRUE
      ),
      width = 3
      ),
    column(
      numericInput(
        "downsamplingSeedDE",
        label = "Set Seed",
        value = 1234,
        min=1,
        max=100000,
        step=1
      ),
      width = 3,
      style='padding:0px;'
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

# if diffcyt should be exectued on selected markers (markers_to_test)
output$markerToTestSelection <- renderUI({
  req(input$chosenDAMethod)
  if (input$chosenDAMethod %in% c("diffcyt-DS-limma", "diffcyt-DS-LMM", "CyEMD", "CytoGLMM", "CytoGLM", "wilcoxon_median", "t_test")) {
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
    title = "Visualize Differential Analysis Results",
    height = plotbox_height,
    width = 12,
  )
})

# heatmap can be visualized for subset of clusters (DA)
output$visClusterSelection <- renderUI({
  out <- reactiveVals$DEruns[[input$deVisMethod]]

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
    showNotification("You have not transformed your data. We will assume that you have given us already transformed data as input.", type = "warning", duration = 10)
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
                                   c("marker", "cluster_id"), 
                                   multiple = F),
                    uiOutput("deBoxK"),
                    selectizeInput("deBoxFeatures",
                                   "Markers:",
                                   markers, 
                                   selected = selected_marker,
                                   multiple = F),
                    selectizeInput(
                      "deBoxColor",
                      "Color and Group by:",
                      c(names(colData(reactiveVals$sce))),# names(cluster_codes(reactiveVals$sce))),
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

output$deBoxK <- renderUI({
  req(input$deBoxFacet == "cluster_id")
  selectizeInput(
    "deBoxK",
    "Clusters",
    rev(names(cluster_codes(reactiveVals$sce))),
    multiple = F
  )
})

output$clusterDEPlot <- renderPlot({
  req(input$deBoxFacet)
  facet_by <- input$deBoxFacet
  if (facet_by == "marker"){
    facet_by <- "antigen"
    k <- 'all'
  } else {
    k <- input$deBoxK
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
    my_height <- 3
    my_width <- 2.5
    facet_info <- reactiveVals$pbExprsPlot %>% ggplot2::ggplot_build() %>% magrittr::extract2('layout') %>% magrittr::extract2('layout')
    nr_row <- max(facet_info$ROW)
    nr_col <- max(facet_info$COL)
    ggsave(file, plot = reactiveVals$pbExprsPlot, width=my_width*nr_col, height=my_height*nr_row)
    waiter_hide(id="app")
  }
)


# Observer-----

observe({
  if (reactiveVals$current_tab==6){
    shinyjs::hide("DEVisualization")
    shinyjs::hide("selectionBoxDE")
    shinyjs::hide("deTopTable")
    req(reactiveVals$pbExprsPlot)
    shinyjs::show("selectionBoxDE")
    req(reactiveVals$DEruns)
    shinyjs::show("DEVisualization")
    req(reactiveVals$visExpClicked)
    shinyjs::show("deTopTable")
  }
})

observeEvent(input$downsampling_Yes_No_DE, {
  req(input$downsampling_Yes_No_DE)
  if(input$downsampling_Yes_No_DE == "No"){
    shinyjs::disable("downsamplingNumberDE")
    shinyjs::disable("downsampling_per_sampleDE")
    shinyjs::disable("downsamplingSeedDE")
  }else{
    shinyjs::enable("downsamplingNumberDE")
    shinyjs::enable("downsampling_per_sampleDE")
    shinyjs::enable("downsamplingSeedDE")
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
  reactiveVals$visExpClicked <- TRUE

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
    topTableOut <- data.frame(out)
    eff_r <- isolate(reactiveVals$eff_r)[[visMethod]]
    eff_r[, marker_id := sapply(strsplit(eff_r$group2,'::'), "[", 1)]
    topTableOut <- merge(topTableOut, eff_r[, c("cluster_id", "marker_id", "overall_group","effsize", "magnitude")], by = c("cluster_id", "marker_id"), all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE)
    colnames(topTableOut) <- c(colnames(data.frame(out)), "overall_group","cohens_d", "magnitude")
    
    reactiveVals$topTable <- topTableOut
    topTableOut$p_val <- formatC(topTableOut$p_val)
    topTableOut$p_adj <- formatC(topTableOut$p_adj)
    topTableOut$cohens_d <- formatC(topTableOut$cohens_d)
    if(visMethod == "diffcyt-DS-limma"){
      topTableOut$logFC <- formatC(topTableOut$logFC)
      topTableOut$AveExpr <- formatC(topTableOut$AveExpr)
      topTableOut$t <- formatC(topTableOut$t)
      topTableOut$B <- formatC(topTableOut$B)
      
    }
    DT::datatable(topTableOut, rownames = FALSE, 
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
