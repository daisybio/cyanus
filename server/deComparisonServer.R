library(ggvenn)

resetDEComparison <- function(){
  shinyjs::reset("column1_comparison")
  shinyjs::reset("column2_comparison")
  reactiveVals$ds_bool <- NULL
  reactiveVals$lastVenn <- NULL
  reactiveVals$lastAllResults <- NULL
}

# --------------------------------------------------------------
# Main function: -------------------
# parses inputs and runs all methods by calling runDA / runDS
# ---------------------------------------------------------------------------------

runMethods <- function(){
  resultVenn <- list()
  parameters <- list()
  reactiveVals$ds_bool <- T
  sce <- isolate(reactiveVals$sce)
  ei <- metadata(sce)$experiment_info
  nr_samples <- length(levels(colData(sce)$sample_id))
  
  condition <- isolate(input$conditionInComp)
  group <- isolate(input$groupColComp)
  if(group == "") group <- NULL
  addTerms <- isolate(input$addTermsComp)
  clusters <- isolate(input$deClusterVenn)
  
  #check if downsampling should be performed
  if(input$downsampling_Yes_No_Comp == "Yes"){
    downsamplingNumberComp <- isolate(input$downsamplingNumberComp)
    downsampling_per_sampleComp <- isolate(input$downsampling_per_sampleComp)
    if(downsampling_per_sampleComp == "Yes"){
      downsampling_per_sampleComp <- TRUE
    }else{
      downsampling_per_sampleComp <- FALSE
    }
    downsamplingSeedComp <- isolate(input$downsamplingSeedComp)
    showNotification("Subsampling the data...", type = "warning")
    sce <- performDownsampling(sce, downsampling_per_sampleComp, downsamplingNumberComp, downsamplingSeedComp)
    if(is.null(sce)){
      return(NULL)
    }
    
  }
  
  if(is.null(input$deSubselectionComp)){
    subselection <- "No"
  }else{
    subselection <- isolate(input$deSubselectionComp)
  }
  
  if (subselection != "No") {
    excluded <- list()
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
          excluded,
          method # TODO there is no method here. Test this again
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
    remove(returned)
  }
  
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    reactiveVals$ds_bool <- F
    # ------------
    # DA: get parameters -------------------
    # -------------------------------
    colsDesignDA <- c(condition, group, addTerms)
    edgeRTrend <- isolate(input$edgeR_trendMethodVenn)
    colsFixedDA <- c(condition, addTerms)
    colsRandomDA <- group
    blockIDVoom <- isolate(input$blockID_voomVenn)
    if(blockIDVoom != ""){
      blockIDVoom <- metadata(sce)$experiment_info[[blockIDVoom]]
    }else{
      blockIDVoom <- NULL
    }
    normalize <- isolate(input$normalizeDEVenn)
    if(normalize == "Yes"){
      normalize <- TRUE
    }else{
      normalize <- FALSE
    }

    parameters[["diffcyt-DA-edgeR"]] <- prepDiffExp(
      sce = sce, 
      contrastVars = condition,
      colsDesign = colsDesignDA,
      method = "diffcyt-DA-edgeR"
    )
    parameters[["diffcyt-DA-voom"]] <-  prepDiffExp(
      sce = sce, 
      contrastVars = condition,
      colsDesign = colsDesignDA,
      method = "diffcyt-DA-voom"
    )
    
    parameters[["diffcyt-DA-GLMM"]] <- prepDiffExp(
      sce = sce, 
      contrastVars = condition,
      colsFixed = colsFixedDA, 
      colsRandom = colsRandomDA,
      method = "diffcyt-DA-GLMM"
    )
    
    # -------------
    # DA: control inputs ------------------
    # -------------------------------
    
    if(ncol(parameters[["diffcyt-DA-edgeR"]][["design"]]) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      return(NULL)
    }
    
    if(nrow(parameters[["diffcyt-DA-GLMM"]][["contrast"]]) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      return(NULL)
    }
    
    if(isolate(input$blockID_voomVenn) %in% colsDesignDA){
      showNotification("Please don't put your blocking variable in the design matrix. See our tooltip for more information", type = "error")
      return(NULL)
    }
    
    # -------------------------------
    # DA: run all methods
    # -------------------------------
  
    resultVenn <- runDA(sce = sce, 
          parameters = parameters,
          da_methods = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"),
          clustering_to_use = clusters, 
          normalize = normalize, 
          trend_edgeR = edgeRTrend, 
          blockID = blockIDVoom
    )

    return(resultVenn)
    
  }else{
    # -----------
    # DS: get parameters --------------------
    # -------------------------------
    markersToTest <- isolate(input$DEFeaturesInVenn)
    methods <- isolate(input$chosenDAMethodComp)
    
    blockIDLimma <- NULL
    limmaTrend <- NULL
    cyEMD_binsize <- NULL
    if ('limma' %in% input$chosenDAMethodComp){
      blockIDLimma <- isolate(input$blockID_limmaVenn)
      if(blockIDLimma != ""){
        blockIDLimma <- metadata(sce)$experiment_info[[blockIDLimma]]
      }
      limmaTrend <- ifelse(isolate(input$trend_limmaVenn) == "Yes", TRUE, FALSE)
    }
    includeWeights <- isolate(input$weightsSelectionVenn)
    includeWeights <- ifelse(includeWeights == "Yes", TRUE, FALSE)
    if ("CyEMD" %in% methods) {
      cyEMD_binsize <- isolate(input$emdBinwidthComp)
      if (cyEMD_binsize == 0)
        cyEMD_binsize <- NULL
      emdNperm <- isolate(input$emdNpermComp)
    }
    if ("CytoGLM" %in% methods) {
      CytoGLM_num_boot <- isolate(input$CytoGLM_num_bootComp)
    }
    if("CytoGLMM" %in% methods & is.null(group)){
      shiny::showNotification("Results of CytoGLMM are not meaningful when no grouping variable like patient ID is selected.", 
                              type = "warning", duration = NULL)
    }
    
    
    # -------------------------------
    # DS: run all methods
    # -------------------------------
    
    showNotification(
      ui =
        HTML(
          "<div id='emdProgress'><b>Progress:</b><div>"
        ),
      duration = NULL,
      id = "emdProgressNote"
    )

    withCallingHandlers({
      resultVenn <- runDS(
        sce = sce,
        ds_methods = methods,
        clustering_to_use = clusters,
        contrast_vars = condition,
        markers_to_test = markersToTest,
        parameters = parameters,
        blockID = blockIDLimma,
        trend_limma = limmaTrend,
        include_weights = includeWeights,
        design_matrix_vars = c(condition, addTerms, group), 
        fixed_effects = c(condition, addTerms), 
        random_effects = group,
        cyEMD_nperm = emdNperm, 
        cyEMD_binsize = cyEMD_binsize,
        cytoGLMM_num_boot = CytoGLM_num_boot,
        time_methods = FALSE,
        parallel = FALSE
      )
      #the effect sizes do not have to be computed multiple times
      reactiveVals$eff_r[["comparison"]] <- findEffectSize(sce, condition, group, clusters)
      resultVenn[["effect_size"]] <- reactiveVals$eff_r[["comparison"]]
    },
    message = function(m) {
      shinyjs::html(id = "emdProgress",
                    html = sprintf("<br>%s", HTML(m$message)),
                    add = TRUE)
    })
    
    removeNotification("emdProgressNote")
    return(resultVenn)
  }
}

# ---------------------------------------------------------------------------------
# Renderer
# ---------------------------------------------------------------------------------

output$vennDiagrams <- renderPlot({
  ggplotObject <- ggplot() + theme_void()
  return(ggplotObject)
})

output$vennTitle <- renderUI({
  div(
    ""
  )
})
shinyjs::hide("vennDiagramsBox")

# column 1 -----

output$DSVenn <- renderUI({
  req(input$da_dsVenn == "Differential Marker Expression")
  div(
  pickerInput(
    "chosenDAMethodComp",
    choices = methodsDS,
    label = span("Available Methods", icon("question-circle"), id = "deMethodsCompQ"),
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  ),
  bsPopover(
    id = "deMethodsCompQ",
    title = "Available Methods",
    content = "Depending on what you want to analyse, there are different methods available. Please see their documentation in the DE tab for further explanation."
  ),
  id = "DS_Venn_div"
  )
})


output$conditionSelectionComp <- renderUI({
  sceEI <- CATALYST::ei(reactiveVals$sce)
  condChoices <- which(sapply(sceEI, function(feature) nlevels(as.factor(feature)) == 2))
  if (length(condChoices) == 0) {
    showNotification("No condition with exactly two levels found. Unfortunately we currently only support comparisons between two conditions. You might want to subset your data.", duration = NULL, type = "error")
    return(NULL)
  }
  condChoices <- names(condChoices)
  list(div(
    selectizeInput(
      "conditionInComp",
      choices = condChoices,
      label = span(
        "What condition do you want to analyse?",
        icon("question-circle"),
        id = "conditionInQComp"
      )
    ),
    bsPopover(
      id = "conditionInQComp",
      title = "Condition for DE analysis",
      content = HTML("Here, you specify the comparison of interest.<br><b>Currently only conditions with two levels are supported.</b>")
    )))
})


output$groupSelectionComp <- renderUI({
  all_methods <- c(methodsDA, methodsDS)
  req(input$conditionInComp)
  if (input$da_dsVenn == "Differential Marker Expression")
    any(input$chosenDAMethodComp %in% all_methods[all_methods != 'CyEMD'])
  sceEI <- data.table::as.data.table(CATALYST::ei(reactiveVals$sce))
  groupCol <- names(sceEI)[!names(sceEI) %in% c("n_cells", "sample_id")]
  groupCol <- groupCol[sapply(groupCol, function(x) sceEI[, .(e2 = data.table::uniqueN(get(input$conditionInComp)) == 2),, by=get(x)][, all(e2)])]
  names(groupCol) <- groupCol
  groupCol <- c('unpaired samples' = '', groupCol)
  div(
    pickerInput(
      "groupColComp",
      choices = groupCol,
      label = span(
        "Do you have paired samples? Which column identifies the group e.g. patient_id?",
        icon("question-circle"),
        id = "groupColQComp"
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
      id = "groupColQComp",
      title = "Grouping Variables",
      content = "All methods except for CyEMD are able to take a grouping variable like patient ID into account. CytoGLMM only yields meaningful results when a grouping variable is specified. For limma, the grouping variable is included as an additional fixed effect. LMM, CytoGLMM and CytoGLM include this variable as a random effect. When a grouping variable is specified for the statistical tests, a paired test and a Wilcoxon signed rank test are computed, respectively."
    )
  )
})


output$additionalTermsSelectionComp <- renderUI({
  if (input$da_dsVenn == "Differential Marker Expression"){
    req(input$chosenDAMethodComp, any(startsWith(input$chosenDAMethodComp, 'diffcyt') | startsWith(input$chosenDAMethodComp, 'CytoGLM'))) # this means this is a linear model and additional terms are allowed
  }
  addTerms <- names(ei(reactiveVals$sce))
  addTerms <- addTerms[!addTerms %in% c("n_cells", "sample_id", input$conditionInComp, input$groupColComp)]
  div(
    pickerInput(
      "addTermsComp",
      choices = addTerms,
      label = span(
        "Additional fixed terms to include in the Model",
        icon("question-circle"),
        id = "addTermsQComp"
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
      id = "addTermsQComp",
      title = "Additional terms to include in the Model",
      content = "Since you have specified a method using a linear model (limma, LMM, CytoGLMM, CytoGLM), you can include additional terms that might have an effect on expression other than the condition, such as batches."
    )
  )
})


output$emdInputComp <- renderUI({
  req("CyEMD" %in% input$chosenDAMethodComp)
  sceEI <- ei(reactiveVals$sce)
  list(
    uiOutput("emdNpermInputComp"),
    div(
      numericInput(
        "emdBinwidthComp",
        label = span(
          "CyEMD: Bin width for comparing histograms",
          icon("question-circle"),
          id = "emdBinwidthQComp"
        ),
        value = 0,
        min = 0,
        max = 1,
        step = .1
      ),
      bsPopover(
        id = "emdBinwidthQComp",
        title = "Bin width for comparing histograms",
        content = HTML("You can set a custom binwidth but we recommend to leave this at zero.<br><b>Set this to 0 to compute the binwidth for each marker based on the Freedman-Diaconis rule.</b>")
      )
    ))
})


output$emdNpermInputComp <- renderUI({
  req(input$conditionInComp)
  maxPerm <- as.numeric(RcppAlgos::permuteCount(ei(reactiveVals$sce)[[input$conditionInComp]]))
  div(
    numericInput(
      "emdNpermComp",
      label = span(
        "CyEMD: Number of permutations for p-value estimation",
        icon("question-circle"),
        id = "emdNpermQComp"
      ),
      value = min(500, maxPerm),
      min = 0,
      max = maxPerm,
      step = 100
    ),
    bsPopover(
      id = "emdNpermQComp",
      title = "Number of permutations for p-value estimation",
      content = HTML("Note that meaningful results require many permutations. E.g. for an unadjusted pvalue smaller than 0.01 at least 100 permutations are necessary.<br><b>This value must not exceed the factorial of the number of samples.</b>")
    )
  )
})


output$CytoGLM_num_bootComp <- renderUI({
  req("CytoGLM" %in% input$chosenDAMethodComp)
  div(
    numericInput(
      "cytoNBootComp",
      label = span(
        "CytoGLM: Number of bootstrap samples",
        icon("question-circle"),
        id = "cytoNBootQComp"
      ),
      value = 1000,
      min = 0,
      max = 10000,
      step = 100
    ),
    bsPopover(
      id = "cytoNBootQComp",
      title = "Number of bootstrap samples",
      content = HTML('CytoGLM uses bootstrapping with replacement to preserve the cluster structure in donors. For more information refer to <a href="https://doi.org/10.1186/s12859-021-04067-x" target="_blank">Seiler et al.</a><br><b>Setting this number very high has a great influence on runtime.</b>')
    )
  )
})


output$deSubselectionComp <- renderUI({
  choices <- isolate(colnames(metadata(reactiveVals$sce)$experiment_info))
  choices <- choices[!choices %in% c("n_cells", "sample_id", "patient_id")]
  
  choices <- unlist(sapply(choices, function(x){
    lvls <- isolate(levels(metadata(reactiveVals$sce)$experiment_info[[x]]))
    return(lvls)
  }))
  names(choices) <- paste("only", choices)
  div(
    checkboxGroupInput(
      inputId = "deSubselectionComp",
      label = span("Do you want to analyse this condition just on a subset?", icon("question-circle"), id = "subSelectVennQ"),
      choices = choices, 
      inline = T
    ),
    bsPopover(
      id = "subSelectVennQ",
      title = "Run differential expression on a subset of your data",
      content = "Sometimes it might make sense to compare differential expression just in a subset of your data, e.g. you have two different treatment groups and want to investigate the effect of an activation agent separately. You can do the subselection right at the beginning (Preprocessing) or here.",
      placement = "top"
    )
  )
})

# until here ----
# begin column 2 ----

output$downsamplingComp <- renderUI({
  req(reactiveVals$sce)
  smallest_n <- min(CATALYST::ei(reactiveVals$sce)$n_cells)
  sum_n <- sum(CATALYST::ei(reactiveVals$sce)$n_cells)
  fluidRow(
    column(
      radioButtons(
        "downsampling_Yes_No_Comp",
        label = span("Do you want to perform downsampling?", icon("question-circle"), id="dsCompPopover"),
        choices = c("Yes", "No"),
        selected = "No",
        inline = TRUE
      ),
      bsPopover(
        id="dsCompPopover",
        title = "Downsample your data",
        content = "If you have a big dataset and do not want to wait too long for your analyses, you can perform a downsampling on your dataset. If you choose to downsample per sample, the number of cells you specify will be randomly picked from each sample. Otherwise, the number you specify will be divided by the number of samples and this number will be randomly picked from each sample. If the number is bigger than the sample size, all cells from this sample will be taken."
      ),
      width = 3
    ),
    column(
      numericInput(
        "downsamplingNumberComp",
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
        "downsampling_per_sampleComp",
        label = "Per sample?",
        choices = c("Yes", "No"),
        inline = TRUE
      ),
      width = 3
    ),
    column(
      numericInput(
        "downsamplingSeedComp",
        label = "Set Seed",
        value = 1234,
        min=1,
        max=100000,
        step=1
      ),
      width = 3
    )
  )
})


output$clusterSelectionComp <- renderUI({
  clusters <- rev(names(cluster_codes(reactiveVals$sce)))
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    clusters <- clusters[!clusters %in% c("all")]
  }
  selectizeInput(
    inputId = "deClusterVenn",
    label = "Choose the cluster populations you want to compare",
    choices = clusters, 
    multiple = F
  )
})


output$markerToTestSelectionComp <- renderUI({
  req(input$da_dsVenn)
  if(input$da_dsVenn == "Differential Marker Expression"){
    div(
      selectInput(
        "DEMarkerToTestVenn",
        label = "Features to choose from",
        choices = c("Marker by Class",
                    "Marker by Name")
      ),
      uiOutput("DEFeatureSelectionVenn")
    )
  }
})


output$DEFeatureSelectionVenn <- renderUI({
  req(input$DEMarkerToTestVenn)
  if (input$DEMarkerToTestVenn == "Marker by Class") {
    choices <-
      levels(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
    choices <- choices[choices %in% c("type", "state")]
    if("state" %in% choices){
      selected <- "state"
    }else{
      selected <- choices[1]
    }
  } else if (input$DEMarkerToTestVenn == "Marker by Name") {
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
    inputId = "DEFeaturesInVenn",
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


output$extraFeaturesComp <- renderUI({
  req(input$da_dsVenn)
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    cols <- colnames(metadata(reactiveVals$sce)$experiment_info)
    cols <- cols[!cols %in% c("n_cells", "sample_id")]
    div(
      selectizeInput(
        inputId = "edgeR_trendMethodVenn",
        label = span("EdgeR Trend Method", icon("question-circle"), id = "trendMethodVennQ"),
        choices = c( "none", "locfit", "movingave", "loess", "locfit.mixed"),
        multiple = F
      ),
      bsPopover(
        id = "trendMethodVennQ",
        title = "Method for estimating dispersion trend",
        content = "EdgeR specific parameter for estimating the dispersion trend. See more in their estimateDisp function documentation. By default, we use the option none to calculate dispersions, since the dispersion-mean relationship typically does not resemble RNA-sequencing data.",
        placement = "top"
      ),
      selectizeInput(
        inputId = "blockID_voomVenn",
        label = span("Voom: Block ID", icon("question-circle"), id = "blockIDvoomVennQ"),
        choices = cols,
        options = list(
          placeholder = "Select your block IDs or nothing",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = F
      ),
      bsPopover(
        id = "blockIDvoomVennQ",
        title = "Block IDs if you have replicates",
        content = "Block IDs (e.g. patient IDs if you have replicates in the same conditions) for paired experimental designs, to be included as random effects (for method testDA_voom or testDS_limma). If provided, the block IDs will be included as random effects using the limma duplicateCorrelation methodology. Alternatively, block IDs can be included as fixed effects in the design matrix.",
        placement = "top"
      ),
      
    )
  }else{
    req('diffcyt-DS-limma' %in% input$chosenDAMethodComp)
    cols <- colnames(metadata(reactiveVals$sce)$experiment_info)
    cols <- cols[!cols %in% c("n_cells", "sample_id")]
    div(
      selectizeInput(
        inputId = "blockID_limmaVenn",
        label = span("Limma: Block ID", icon("question-circle"), id = "blockIDlimmaQVenn"),
        choices = cols,
        options = list(
          placeholder = "Select your block IDs or nothing",
          onInitialize = I("function() { this.setValue(''); }")
        ),
        multiple = F
      ),
      bsPopover(
        id = "blockIDlimmaQVenn",
        title = "Block IDs if you have replicates",
        content = "Block IDs (e.g. patient IDs if you have replicates in the same conditions) for paired experimental designs, to be included as random effects (for method testDA_voom or testDS_limma). If provided, the block IDs will be included as random effects using the limma duplicateCorrelation methodology. Alternatively, block IDs can be included as fixed effects in the design matrix.",
        placement = "top"
      ),
      radioButtons(
        inputId = "trend_limmaVenn",
        label = span("Limma Trend Method", icon("question-circle"), id = "trendlimmaQVenn"),
        choices = c("Yes", "No"),
        inline = T
      ),
      bsPopover(
        id = "trendlimmaQVenn",
        title = "Limma Trend Method",
        content = "Whether to fit a mean-variance trend when calculating moderated tests with function eBayes from limma package (for method testDS_limma). When trend = TRUE, this is known as the limma-trend method (Law et al., 2014; Phipson et al., 2016).",
        placement = "top"
      )
    )
  }
})


output$normalizeSelectionComp <- renderUI({
  req(input$da_dsVenn)
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    div(
      radioButtons(
        inputId = "normalizeDEVenn",
        label = span("Normalize?", icon("question-circle"), id = "normalizeDEQVenn"),
        choices = c("Yes", "No"),
        inline = T
      ),
      bsPopover(
        id = "normalizeDEQVenn",
        title = "Composition Effects",
        content = "Whether to include optional normalization factors to adjust for composition effects. Only relevant for Differential Cluster Abundance methods."
      )
    )
  }
})


# selection of weights for analysis
output$weightSelectionComp <- renderUI({
  req(input$da_dsVenn)
  if (input$da_dsVenn == "Differential Marker Expression"){
    div(
      radioButtons(
        inputId = "weightsSelectionVenn",
        label = span("Do you want to include precision weights (cell counts) within the model?", icon("question-circle"), id = "weightSelectVennQ"),
        choices = c("Yes", "No"), 
        inline = T
      ),
      bsPopover(
        id = "weightSelectVennQ",
        title = "Whether to include precision weights within each model (across samples).",
        content = "These represent the relative uncertainty in calculating each median value. The cell counts of each sample are incorporated in the analysis.",
        placement = "top"
      )
    )
  }
})


output$fdrComp <- renderUI({
  div(
    numericInput(
      inputId = "fdrThresholdVenn",
      label = "FDR threshold",
      value = 0.05,
      min = 0.0,
      max = 1.0,
      step = 0.01
    )
  )
})

### Download buttons

output$downloadVenn <- renderUI({
  req(reactiveVals$lastVenn)
  div(
    downloadButton("downloadVennButton", "Download Plot"),
    style = "position: absolute; bottom: 5px; right:5px"
  )
})

output$downloadVennButton <- downloadHandler(
  filename = function(){
    paste0("VennDiagram", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$lastVenn, width=12, height=12)
    waiter_hide(id="app")
  }
)


output$downloadTableVenn <- renderUI({
  req(reactiveVals$lastAllResults)
  fluidRow(
    div(
      downloadButton("downloadTableVennAll", "Download All Results"),
      style = "position: absolute; top: 5px; right:5px"
    )
  )
})


output$downloadTableVennAll <- downloadHandler(
  filename = "AllResults.csv",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    write.csv(reactiveVals$lastAllResults, file, row.names = FALSE)
    waiter_hide(id="app")
  }
)

# ---------------------------------------------------------------------------------
# Observer
# ---------------------------------------------------------------------------------

observeEvent(input$downsampling_Yes_No_Comp, {
  req(input$downsampling_Yes_No_Comp)
  if(input$downsampling_Yes_No_Comp == "No"){
    shinyjs::disable("downsamplingNumberComp")
    shinyjs::disable("downsampling_per_sampleComp")
    shinyjs::disable("downsamplingSeedComp")
  }else{
    shinyjs::enable("downsamplingNumberComp")
    shinyjs::enable("downsampling_per_sampleComp")
    shinyjs::enable("downsamplingSeedComp")
  }
})

observeEvent(input$diffExpButtonVenn, {
  req(input$da_dsVenn)
  if (input$da_dsVenn == "Differential Marker Expression" & is.null(isolate(input$chosenDAMethodComp))) {
    showNotification("No method selected. Try again.", type = "error")
    return(NULL)
  }
  waiter_show(id = "app",html = tagList(spinner$logo, 
                             HTML("<br>DE Analysis in Progress...<br>Please be patient")), 
              color=spinner$color)
  resultVenn <- runMethods()
  if(!is.null(resultVenn)){
      ds_bool <- isolate(reactiveVals$ds_bool)
    output$vennTitle <- renderUI({
      div(
        paste("Significant results for FDR threshold", isolate(input$fdrThresholdVenn))
      )
    })
    output$vennDiagrams <- renderPlot({
      venn <- createVennHeatmap(resultVenn, ds_bool, isolate(input$fdrThresholdVenn))
      reactiveVals$lastVenn <- venn
      venn
    })
    output$vennTable <- renderUI({
      req(resultVenn)
      if(ds_bool){
        firstCol <- c("marker_id", "cluster_id")
      }else{
        firstCol <- c("cluster_id")
      }
      library(data.table)
      listDT <- lapply(resultVenn[names(resultVenn) != "effect_size"], function(x) {
        vec <- c(firstCol, "p_val", "p_adj")
        as.data.table(x)[, ..vec]
      })
      allResultsDT <- rbindlist(listDT, idcol = "method")
      if(ds_bool){
        eff_r <- resultVenn[["effect_size"]]
        eff_r[, marker_id := sapply(strsplit(eff_r$group2,'::'), "[", 1)]
        allResultsDT <- merge(allResultsDT, eff_r[, c("cluster_id", "marker_id", "overall_group","effsize", "magnitude")], by = c("cluster_id", "marker_id"), all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE)
        colnames(allResultsDT) <- c("cluster_id", "marker_id", "method", "p_val", "p_adj", "overall_group","cohens_d", "magnitude")
      }
      reactiveVals$lastAllResults <- allResultsDT
      
      shinydashboard::box(
        div(
          renderDataTable(
            DT::datatable(
              allResultsDT,
              rownames = F,
              options = list(pageLength = 10, searching = FALSE, 
                             columnDefs = list(list( targets = "_all", 
                                                     render = JS("function(data, type, row, meta) {","return data === null ? 'NA' : data;","}"))))
            )
          )
        ),
        div(
          uiOutput("downloadTableVenn")
        ),
        id = "vennResultsTable",
        title = "Results",
        width = 12
      )
    })
  }
  waiter_hide(id = "app")
  shinyjs::show("vennDiagramsBox")
  shinyjs::enable("diffExpButtonVenn")
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
})




