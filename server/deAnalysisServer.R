library(diffcyt)

observeEvent(input$deBoxFacet, {
  if(input$deBoxFacet == "cluster_id"){
    shinyjs::show("deBoxK")
  }else{
    shinyjs::hide("deBoxK")
  }
})

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

output$clusterSelection <- renderUI({
  selectizeInput(
    inputId = "deCluster",
    label = "Choose the cluster populations you want to compare",
    choices = names(cluster_codes(reactiveVals$sce)), 
    multiple = F
  )
})

output$markerSelection <- renderUI({
  selectizeInput(
    inputId = "deMarker",
    label = "Choose the markers you want to compare",
    # make reactiveVals of features where a cluster id is provided
    choices = input$featuresIn, 
    multiple = F
  )
})

output$modelSelection <- renderUI({
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DS-limma","diffcyt-DA-voom")){
    uiOutput("designMatrixSelection")
  } else {
    uiOutput("formulaSelection")
  }
})

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

observeEvent(input$diffExpButton,{
  shinyjs::disable("diffExpButton")
  ei <- metadata(reactiveVals$sce)$experiment_info
  contrast <- createContrast(c(0, 1))
  if (input$chosenDAMethod %in% c("diffcyt-DA-edgeR","diffcyt-DA-voom")){
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      design = design,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      markers_to_test = input$deMarker
    )
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-limma")){
    design <- createDesignMatrix(ei, cols_design = input$colsDesign)
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      design = design,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      markers_to_test = input$deMarker
    )
  } else if (input$chosenDAMethod %in% c("diffcyt-DS-LMM")){
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DS = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      markers_to_test = input$deMarker
    )
  } else if (input$chosenDAMethod %in% c("diffcyt-DA-GLMM")){
    formula <- createFormula(ei, cols_fixed = input$colsFixed, cols_random = input$colsRandom)
    out <- diffcyt::diffcyt(
      d_input = reactiveVals$sce,
      formula = formula,
      contrast = contrast,
      analysis_type = reactiveVals$methodType,
      method_DA = input$chosenDAMethod,
      clustering_to_use = input$deCluster,
      markers_to_test = input$deMarker
    )
  }
  
  reactiveVals$out <- out
  
  output$heatmapDEPlot <- renderPlot({
    reactiveVals$diffHeatmapPlot <- plotDiffHeatmap(x=reactiveVals$sce,y=rowData(reactiveVals$out$res), k=input$deCluster,fdr=as.numeric(input$fdrThreshold), lfc=as.numeric(input$lfcThreshold), sort_by = input$heatmapSortBy, all=TRUE, normalize=as.logical(input$heatmapNormalize))
    reactiveVals$diffHeatmapPlot
  })
  
  output$deHeatmapBox <- renderUI({
    fluidRow(column(4,
                    div(dropdownButton(
                      tags$h3("Plot Options"),
                      textInput("fdrThreshold", "Fdr Threshold:", value ="0.05"),
                      textInput("lfcThreshold", "Log2FC Threshold:", value ="1"),
                      selectizeInput(
                        "heatmapSortBy",
                        "Sort by:",
                        choices = c("P-adjusted"="padj", "LogFC"="lfc", "None"="none"),
                        multiple = F
                      ),
                      selectizeInput(
                        "heatmapNormalize",
                        "Z-score normalization:",
                        c("Yes"="TRUE","No"="FALSE"),
                        multiple=F
                      ),
                      circle = TRUE,
                      status = "info",
                      icon = icon("gear"),
                      width = "100%",
                      tooltip = tooltipOptions(title = "Click to see plot options")
                    ),
                    div(
                      uiOutput("heatmapPlotDownload"),
                      style = "position: absolute; bottom: 5px;"
                    ),
                    style = "position: relative; height: 550px;"
                    ),
              ),
             column(8, shinycssloaders::withSpinner(
               plotOutput("heatmapDEPlot", width = "100%", height = "550px")
             )), 
             )
    
  })
  
  # ui for download button
  output$heatmapPlotDownload <- renderUI({
    req(reactiveVals$diffHeatmapPlot)
    downloadButton("downloadPlotDiffHeatmap", "Download Plot")
  })
  
  # function for downloading MDS plot
  output$downloadPlotDiffHeatmap <- downloadPlotFunction("Heatmap_DiffExp_Plot", reactiveVals$diffHeatmapPLot)
  
  shinyjs::enable("diffExpButton")
  
})

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
  fluidRow(column(4,
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
                    width = "100%",
                    tooltip = tooltipOptions(title = "Click to see plot options")
                  ),
                  div(
                    uiOutput("pbExprsPlotDownload"),
                    style = "position: absolute; bottom: 5px;"
                  ),
                  style = "position: relative; height: 550px;"
                  )
                ),
           column(8, shinycssloaders::withSpinner(
             plotOutput("clusterDEPlot", width = "100%", height = "500px")
           )))
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










