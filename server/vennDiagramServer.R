library(ggvenn)

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
#shinyjs::hidden("vennTable")

output$modelSelectionVenn <- renderUI({
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    uiOutput("DAVenn")
  }else{
    uiOutput("DSVenn")
  }
})

output$DAVenn <- renderUI({
  colsDesign <- colnames(metadata(reactiveVals$sce)$experiment_info)
  colsDesign <- colsDesign[!colsDesign %in% c("n_cells", "sample_id")]
  div(
    pickerInput(
      "colsDesignDA",
      choices = colsDesign,
      selected = colsDesign[1],
      label = span(
        "EdgeR, Voom: Which columns to include in the design matrix?",
        icon("question-circle"),
        id = "deDesignMatrixVenn"
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
      id = "deDesignMatrixVenn",
      title = "Design matrix for model fitting",
      content = "edgeR and voom work with a design matrix: The selected columns will be included in the design matrix specifying the models to be fitted. A design matrix includes all variables that are relevant/interesting for your analysis because the linear model will be built using these variables. That means that you MUST include the condition you want to analyze here."
    ),
    pickerInput(
      "colsFixedDA",
      choices = colsDesign,
      selected = colsDesign[1],
      label = span(
        "GLMM: Which fixed effect terms to include in the model formula?",
        icon("question-circle"),
        id = "deFormulaFixVenn"
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
      id = "deFormulaFixVenn",
      title = "Fixed effect terms for the model formula",
      content = "GLMM works with fixed and random effects: Fixed effects are variables that we expect will have an effect on the dependent/response variable: they’re what you call explanatory variables in a standard linear regression. You MUST include the condition variable here."
    ),
    pickerInput(
      "colsRandomDA",
      choices = colsDesign,
      selected = colsDesign[2],
      label = span(
        "GLMM: Which random intercept terms to include in the model formula",
        icon("question-circle"),
        id = "deFormulaRandomVenn"
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
      id = "deFormulaRandomVenn",
      title = "Random intercept terms for the model formula",
      content = "GLMM works with fixed and random effects: Depending on the experimental design, this may include group IDs (e.g. groups for differential testing) or block IDs (e.g. patient IDs in a paired design)."
    )
  )
})

output$DSVenn <- renderUI({
  colsDesign <- colnames(metadata(reactiveVals$sce)$experiment_info)
  colsDesign <- colsDesign[!colsDesign %in% c("n_cells", "sample_id")]
  div(
    pickerInput(
      "colsDesignDS",
      choices = colsDesign,
      selected = colsDesign[1],
      label = span(
        "Limma: Which columns to include in the design matrix?",
        icon("question-circle"),
        id = "deDesignMatrixVennDS"
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
      id = "deDesignMatrixVennDS",
      title = "Design matrix for model fitting",
      content = "limma works with a design matrix: The selected columns will be included in the design matrix specifying the models to be fitted. A design matrix includes all variables that are relevant/interesting for your analysis because the linear model will be built using these variables. That means that you MUST include the condition you want to analyze here."
    ),
    pickerInput(
      "colsFixedDS",
      choices = colsDesign,
      selected = colsDesign[1],
      label = span(
        "LMM: Which fixed effect terms to include in the model formula?",
        icon("question-circle"),
        id = "deFormulaFixVennDS"
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
      id = "deFormulaFixVennDS",
      title = "Fixed effect terms for the model formula",
      content = "LMM works with fixed and random effects: Fixed effects are variables that we expect will have an effect on the dependent/response variable: they’re what you call explanatory variables in a standard linear regression. You MUST include the condition variable here."
    ),
    pickerInput(
      "colsRandomDS",
      choices = colsDesign,
      selected = colsDesign[2],
      label = span(
        "LMM: Which random intercept terms to include in the model formula",
        icon("question-circle"),
        id = "deFormulaRandomVennDS"
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
      id = "deFormulaRandomVennDS",
      title = "Random intercept terms for the model formula",
      content = "LMM works with fixed and random effects: Depending on the experimental design, this may include group IDs (e.g. groups for differential testing) or block IDs (e.g. patient IDs in a paired design)."
    )
  )
})

output$contrastSelectionVenn <- renderUI({
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    uiOutput("DAContrastVenn")
  }else{
    div(
      uiOutput("emdInputVenn"),
      uiOutput("DSContrastVenn")
    )
  }
})

output$DAContrastVenn <- renderUI({
  req(input$colsDesignDA)
  req(input$colsFixedDA)
  choices <- intersect(input$colsDesignDA, input$colsFixedDA)
  div(
    selectInput(
      "contrastVarsDA",
      choices = choices,
      selected = choices[1],
      label = span(
        "What condition do you want to analyse?",
        icon("question-circle"),
        id = "deContrastVennQ"
      ),
      multiple = F
    ),
    bsPopover(
      id = "deContrastVennQ",
      title = "Contrast Matrix Design",
      content = "Here, you specify the comparison of interest. The p-values will be calculated on the basis of this variable, i.e. it will be tested whether the coefficient of this parameter in the model is equal to zero."
    )
  )
})

output$DSContrastVenn <- renderUI({
  req(input$colsDesignDS)
  req(input$colsFixedDS)
  choices <- intersect(input$colsDesignDS, input$colsFixedDS)
  div(
    selectInput(
      "contrastVarsDS",
      choices = choices,
      selected = choices[1],
      label = span(
        "What condition do you want to analyse?",
        icon("question-circle"),
        id = "deContrastVennDSQ"
      ),
      multiple = F
    ),
    bsPopover(
      id = "deContrastVennDSQ",
      title = "Contrast Matrix Design",
      content = "Here, you specify the comparison of interest. The p-values will be calculated on the basis of this variable, i.e. it will be tested whether the coefficient of this parameter in the model is equal to zero."
    )
  )
})

output$emdInputVenn <- renderUI({
  req(input$contrastVarsDS)
  maxPerm <- as.numeric(RcppAlgos::permuteCount(ei(reactiveVals$sce)[[input$contrastVarsDS]]))
  div(
    numericInput(
      "emdNpermVenn",
      label = span(
        "EMD: Number of permutations for p-value estimation",
        icon("question-circle"),
        id = "emdNpermQVenn"
      ),
      value = min(100, maxPerm),
      min = 0,
      max = maxPerm,
      step = 10
    ),
    bsPopover(
      id = "emdNpermQVenn",
      title = "EMD: Number of permutations for p-value estimation",
      content = HTML("Note that meaningful results require many permutations. For an unadjusted pvalue smaller than 0.01 at least 100 permutations are necessary.<br><b>This value must not exceed the factorial of the number of samples.</b>")
    ),
    numericInput(
      "emdBinwidthVenn",
      label = span(
        "Bin width for comparing histograms",
        icon("question-circle"),
        id = "emdBinwidthQVenn"
      ),
      value = 0,
      min = 0,
      max = 1,
      step = .1
    ),
    bsPopover(
      id = "emdBinwidthQVenn",
      title = "Bin width for comparing histograms",
      content = HTML("You can set a custom binwidth but we recommend to leave this at zero.<br><b>Set this to 0 to compute the binwidth for each marker based on the Freedman-Diaconis rule.</b>")
    )
  )
})

output$deSubselectionVenn <- renderUI({
  choices <- isolate(colnames(metadata(reactiveVals$sce)$experiment_info))
  choices <- choices[!choices %in% c("n_cells", "sample_id", "patient_id")]
  
  choices <- unlist(sapply(choices, function(x){
    lvls <- isolate(levels(metadata(reactiveVals$sce)$experiment_info[[x]]))
    return(lvls)
  }))
  names(choices) <- paste("only", choices)
  div(
    checkboxGroupInput(
      inputId = "deSubselectionVenn",
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

output$clusterSelectionVenn <- renderUI({
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

output$markerToTestSelectionVenn <- renderUI({
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

output$extraFeaturesVenn <- renderUI({
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

output$normalizeSelectionVenn <- renderUI({
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

output$fdrVenn <- renderUI({
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

observeEvent(input$diffExpButtonVenn, {
  waiter_show(html = tagList(spinner$logo, 
                             HTML("<br>DE Analysis in Progress...<br>Please be patient")), 
              color=spinner$color)
  shinyjs::disable("diffExpButtonVenn")
  shinyjs::disable("previousTab")
  shinyjs::disable("nextTab")
  req(input$da_dsVenn)
  resultVenn <- runMethods()
  if(!is.null(resultVenn)){
    output$vennTitle <- renderUI({
      div(
        paste("Significant results for FDR threshold", isolate(input$fdrThresholdVenn))
      )
    })
    output$vennDiagrams <- renderPlot({
      ds_bool <- isolate(reactiveVals$ds_bool)
      venn <- createVennDiagram(resultVenn, ds_bool)
      reactiveVals$lastVenn <- venn
      venn
    })
    output$vennTable <- renderUI({
      req(resultVenn)
      if("diffcyt-DA-edgeR" %in% names(resultVenn)){
        firstCol <- c("cluster_id")
      }else{
        firstCol <- c("marker_id", "cluster_id")
      }
      library(data.table)
      
      listDT <- lapply(resultVenn, function(x) {
        vec <- c(firstCol, "p_val", "p_adj")
        as.data.table(x)[, ..vec]
      })
      allResultsDT <- rbindlist(listDT, idcol = "method")
      reactiveVals$lastAllResults <- allResultsDT
      
      fdrThreshold <- isolate(input$fdrThresholdVenn)
      allResultsSign <- allResultsDT[p_adj < fdrThreshold]
      if("diffcyt-DA-edgeR" %in% names(resultVenn)){
        allResultsSign <- dcast(allResultsSign, cluster_id ~ method, value.var = "p_adj")
      }else{
        allResultsSign <- dcast(allResultsSign, marker_id + cluster_id ~ method, value.var = "p_adj")
      }
      allResultsSign <- allResultsSign[rev(order(allResultsSign[, Reduce(`+`, lapply(.SD,function(x) !is.na(x)))]))]
      reactiveVals$lastAllResultsSign <- allResultsSign
      
      shinydashboard::tabBox(
        shiny::tabPanel(
          renderDataTable(
            DT::datatable(
              allResultsSign,
              rownames = F,
              options = list(pageLength = 10, searching = FALSE, 
                             columnDefs = list(list( targets = "_all", 
                                                     render = JS("function(data, type, row, meta) {","return data === null ? 'NA' : data;","}"))))
            )
          ),
          value = "resultsIntersections",
          title = "Venn Diagram Results"
        ),
        shiny::tabPanel(
          renderDataTable(
            DT::datatable(
              allResultsDT,
              rownames = F,
              options = list(pageLength = 10, searching = FALSE, 
                             columnDefs = list(list( targets = "_all", 
                                                     render = JS("function(data, type, row, meta) {","return data === null ? 'NA' : data;","}"))))
            )
          ),
          value = "resultsAllVenn",
          title = "All results"
        ),
        id = "vennResultsTable",
        title = "Results",
        width = 12
      )
    })
  }
  waiter_hide()
  shinyjs::show("vennDiagramsBox")
  shinyjs::enable("diffExpButtonVenn")
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
})

output$downloadVenn <- renderUI({
  req(reactiveVals$lastVenn)
  div(
    downloadButton("downloadVennButton", "Download Plot"),
    style = "float:right;"
  )
})

output$downloadVennButton <- downloadHandler(
  filename = function(){
    paste0("VennDiagram", ".pdf")
  },
  content = function(file){
    ggsave(file, plot = reactiveVals$lastVenn, width=12, height=12)
  }
)

output$downloadTableVenn <- renderUI({
  req(reactiveVals$lastAllResults)
  req(reactiveVals$lastAllResultsSign)
  fluidRow(
    div(
      downloadButton("downloadTableVennAll", "Download All Results"),
      style = "float:right;"
    ),
    div(
      downloadButton("downloadTableSign", "Download Venn Diagram Results"),
      style = "float:right;"
    )
  )
})

output$downloadTableSign <- downloadHandler(
  filename = "VennDiagramResults.csv",
  content = function(file) {
    write.csv(reactiveVals$lastAllResultsSign, file, row.names = FALSE)
  }
)

output$downloadTableVennAll <- downloadHandler(
  filename = "AllResults.csv",
  content = function(file) {
    write.csv(reactiveVals$lastAllResults, file, row.names = FALSE)
  }
)

runMethods <- function(){
  resultVenn <- list()
  reactiveVals$ds_bool <- T
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    contrast <- isolate(input$contrastVarsDA)
  }else{
    contrast <- isolate(input$contrastVarsDS)
  }
  clusters <- isolate(input$deClusterVenn)
  
  if(is.null(input$deSubselectionVenn)){
    subselection <- "No"
  }else{
    subselection <- isolate(input$deSubselectionVenn)
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
      if(x %in% contrast & catCount[[x]] == 1){
        showNotification("You want to analyse a condition you subsetted. That is not meaningful. Try again.", type = "error")
        return(NULL)
      }
    }
    for(cat in names(exclude)){
      sub <- exclude[[cat]]
      print(sprintf("only using %s from the condition %s", sub, cat))
      sce <- filterSCE(sce, get(cat) %in% sub)
    }
  }
  
  ei <- metadata(sce)$experiment_info
  nr_samples <- length(levels(colData(sce)$sample_id))
  
  if(input$da_dsVenn == "Differential Cluster Abundance"){
    reactiveVals$ds_bool <- F
    colsDesignDA <- isolate(input$colsDesignDA)
    
    design <- createDesignMatrix(ei, cols_design = colsDesignDA)
    contrastMatrix <- createCustomContrastMatrix(contrast, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      return(NULL)
    }
    
    colsFixedDA <- isolate(input$colsFixedDA)
    colsRandomDA <- isolate(input$colsRandomDA)
    
    formulaGLMM <- createFormula(ei, cols_fixed = colsFixedDA, cols_random = colsRandomDA)
    contrastGLMM <- createCustomContrastMatrix(contrast, diffcyt::createDesignMatrix(ei, cols_design = colsFixedDA), designMatrix = T)
    if(nrow(contrastGLMM) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      out <- NULL
    }
    
    edgeRTrend <- isolate(input$edgeR_trendMethodVenn)
    
    blockIDVoom <- isolate(input$blockID_voomVenn)
    if(blockIDVoom %in% colsDesignDA){
      showNotification("Please don't put your blocking variable in the design matrix. See our tooltip for more information", type = "error")
      return(NULL)
    }
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
    
    for(method in c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM")){
      print(method)
      showNotification(sprintf("Calculating %s", method))
      if(method == "diffcyt-DA-GLMM"){
        out <- diffcyt::diffcyt(
          d_input = sce,
          formula = formulaGLMM,
          contrast = contrastGLMM,
          analysis_type = "DA",
          method_DA = method,
          clustering_to_use = clusters,
          normalize = normalize
        )
      }else{
        out <- diffcyt::diffcyt(
          d_input = sce,
          design = design,
          contrast = contrastMatrix,
          analysis_type = "DA",
          method_DA = method,
          clustering_to_use = clusters,
          normalize = normalize,
          trend_method = edgeRTrend,
          block_id = blockIDVoom
        )
      }
      resultVenn[[method]] <- rowData(out$res)
    }
    return(resultVenn)
    
  }else{
    colsDesignDS <- isolate(input$colsDesignDS)
    
    design <- createDesignMatrix(ei, cols_design = colsDesignDS)
    contrastMatrix <- createCustomContrastMatrix(contrast, design, designMatrix = T)
    if(ncol(design) >= nr_samples){
      showNotification("You selected more conditions than there are samples which is not meaningful. Try again.", type = "error")
      return(NULL)
    }
    
    colsFixedDS <- isolate(input$colsFixedDS)
    colsRandomDS <- isolate(input$colsRandomDS)
    
    formulaLMM <- createFormula(ei, cols_fixed = colsFixedDS, cols_random = colsRandomDS)
    contrastLMM <- createCustomContrastMatrix(contrast, diffcyt::createDesignMatrix(ei, cols_design = colsFixedDS), designMatrix = T)
    if(nrow(contrastLMM) >= nr_samples){
      showNotification("You selected more conditions than there are samples as fixed effects which is not meaningful. Try again.", type = "error")
      return(NULL)
    }
    
    markersToTest <- isolate(input$DEFeaturesInVenn)
    is_marker <- rowData(sce)$marker_class %in% c("type", "state")
    if (input$DEMarkerToTestVenn == "Marker by Class") {
      markersToTest <- (rowData(sce)$marker_class %in% markersToTest)[is_marker] # type and state (but not none)
    }else{
      markersToTest <- rownames(sce)[is_marker] %in% markersToTest
    }
    
    blockIDLimma <- isolate(input$blockID_limmaVenn)
    if(blockIDLimma %in% input$colsDesignDS){
      showNotification("Please don't put your blocking variable in the design matrix. See our tooltip for more information", type = "error")
      return(NULL)
    }
    if(blockIDLimma != ""){
      blockIDLimma <- metadata(sce)$experiment_info[[blockIDLimma]]
    }else{
      blockIDLimma <- NULL
    }
    
    limmaTrend <- isolate(input$trend_limmaVenn)
    if(limmaTrend == "Yes"){
      limmaTrend <- TRUE
    }else{
      limmaTrend <- FALSE
    }
    for(method in c("diffcyt-DS-limma", "diffcyt-DS-LMM")){
      print(method)
      showNotification(sprintf("Calculating %s", method))
      if(method == "diffcyt-DS-limma"){
        out <- diffcyt::diffcyt(
          d_input = sce,
          design = design,
          contrast = contrastMatrix,
          analysis_type = "DS",
          method_DS = method,
          clustering_to_use = clusters,
          block_id = blockIDLimma,
          trend = limmaTrend,
          markers_to_test = markersToTest
        )
      }else{
        out <- diffcyt::diffcyt(
          d_input = sce,
          formula = formulaLMM,
          contrast = contrastLMM,
          analysis_type = "DS",
          method_DS = method,
          clustering_to_use = clusters,
          markers_to_test = markersToTest,
        )
      }
      resultVenn[[method]] <- rowData(out$res)
    }
    #run EMD
    markersToTest <- isolate(input$DEFeaturesInVenn)
    if (isolate(input$DEMarkerToTestVenn) == "Marker by Class") {
      sce <- filterSCE(sce, marker_classes(sce) %in% markersToTest)
    }else{
      sce <- filterSCE(sce, rownames(sce) %in% markersToTest)
    }
    
    binSize <- isolate(input$emdBinwidthVenn)
    nperm <- isolate(input$emdNpermVenn)
    if (binSize == 0) binSize <- NULL
    
    showNotification(
      ui =
        HTML(
          "<div id='emdProgress'><b>EMD Progress:</b><div>"
        ),
      duration = NULL,
      id = "emdProgressNote"
    )
    #print(paste("k=", clusters, ", condition=", contrast, ", binsize =", binSize, ", nperm=", nperm))
    withCallingHandlers({
      out <-
        sceEMD(
          sce = sce,
          k = clusters,
          condition = contrast,
          binSize = binSize,
          nperm = nperm
        )
    },
    message = function(m) {
      shinyjs::html(id = "emdProgress",
                    html = sprintf("<br>%s", HTML(m$message)),
                    add = TRUE)
    })
    removeNotification("emdProgressNote")
    resultVenn[["sceEMD"]] <- out
    return(resultVenn)
  }
}

createVennDiagram <- function(res, DS = T) {
  input_venn <- list()
  if(DS){
    feature <- "marker_id"
  }else{
    feature <- "cluster_id"
  }
  # take of each method the data table containing the pvalues
  for (ds_method in names(res)) {
    if(feature == "cluster_id"){
      result <-
          data.frame(res[[ds_method]])[c(feature, "p_val", "p_adj")]
      featureNew <- "cluster_id"
    }else{
      result <-
        data.frame(res[[ds_method]])[c("cluster_id", feature, "p_val", "p_adj")]
      if(result$cluster_id[1] != "all"){
        result$marker_id_joined <- paste0(result$marker_id, "(", result$cluster_id, ")")
        featureNew <- "marker_id_joined"
      }else{
        featureNew <- "marker_id"
      }
      
    }
    result$significant <- result$p_adj < isolate(input$fdrThresholdVenn)
    significants <-
      unlist(subset(
        result,
        significant == TRUE,
        select = c(get(featureNew)),
        use.names = FALSE
      ))
    input_venn[[ds_method]] <- significants
  }
  
  venn <- ggvenn::ggvenn(input_venn, show_elements = TRUE, label_sep ="\n", fill_alpha = 0.3, set_name_size = 6, text_size = 4)
  return(venn)
}



