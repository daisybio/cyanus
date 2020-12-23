# Preprocessing Server

# transform SingleCellExperiment
transformData <-
  function (sce,
            transform,
            cf = 5,
            ain = "counts",
            aout = "exprs") {
    if (transform == "no") {
      y <- assay(sce, ain)
      if (transform == "arcsinh") {
        chs <- channels(sce)
        stopifnot(is.numeric(cf), cf > 0)
        if (length(cf) == 1) {
          int_metadata(sce)$cofactor <- cf
          cf <- rep(cf, nrow(sce))
        }
        else {
          stopifnot(!is.null(names(cf)), chs %in% names(cf))
          cf <- cf[match(chs, names(cf))]
          int_metadata(sce)$cofactor <- cf
        }
        fun <- asinh
        op <- "/"
        y <- fun(sweep(y, 1, cf, op))
      } else if (transform == "log") {
        y <- log(y + 1)
      }
      assay(sce, aout, FALSE) <- y
    }
    return(sce)
  }

# render counts plot
output$countsPlot <- renderPlot({
  CATALYST::plotCounts(
    reactiveVals$sce,
    group_by = input$countsGroupBy,
    color_by = input$countsColorBy
  )
})

# render mds plot
output$mdsPlot <- renderPlot({
  CATALYST::pbMDS(reactiveVals$sce, color_by = input$mdsColorBy, label_by = input$mdsLabels, features = input$nrsFeatures)
  
})

# render nrs plot
output$nrsPlot <- renderPlot({
  CATALYST::plotNRS(reactiveVals$sce,
                    features = input$nrsFeatures,
                    color_by = input$nrsColorBy)
})

# render expression heatmap
output$exprsHeatmapPlot <- renderPlot({
  CATALYST::plotExprHeatmap(reactiveVals$sce, features = input$hpexprFeatures, scale = input$hpexprScale)
})

# render expression plot
output$exprsPlot <- renderPlot({
  CATALYST::plotExprs(reactiveVals$sce, color_by = input$exprsColorBy, features = input$exprsColorBy)
})


# marker selection -> if markers are selected -> visualization can start
observeEvent({
  input$markerSelection
  input$sampleSelection
}, {
  if ((length(input$markerSelection) != 0) &&
      (length(input$sampleSelection) != 0)) {
    updateActionButton(session, "continue", label = "Visualization")
    shinyjs::show("continue")
  }
})

# if start transformation button is clicked
observeEvent(input$prepButton, {
  # data transformation
  fcs <- dirname(input$fcsFiles$datapath)[1]   #later this will be done in upload data
  sce <- prepData(fcs, transform = F) #later this will be done in upload data
  sce <-transformData(
      sce = sce,
      transform = input$transformation,
      cf = as.numeric(input$cofactor)
    )
  reactiveVals$sce <- sce
  
  # show markers and samples selection
  updatePickerInput(
    session,
    'markerSelection',
    choices = names(channels(reactiveVals$sce)),
    selected = names(channels(reactiveVals$sce)) ,
  )
  updatePickerInput(
    session,
    'sampleSelection',
    choices = unique(sample_ids(reactiveVals$sce)),
    selected = unique(sample_ids(reactiveVals$sce)),
  )
})

output$designCounts <- renderUI({
  shinydashboard::box(
    selectizeInput("countsGroupBy",
                   "Group by:",
                   names(colData(
                     reactiveVals$sce
                   )), multiple = F),
    selectizeInput("countsColorBy",
                   "Color by:",
                   names(colData(
                     reactiveVals$sce
                   )), multiple = F),
    width = 4,
  )
})

output$designNrs <- renderUI({
  shinydashboard::box(
    selectizeInput(
      "nrsFeatures",
      "Features:",
      unique(rowData(reactiveVals$sce)$marker_class),
      multiple = F
    ),
    selectizeInput("nrsColorBy",
                   "Color by:",
                   names(colData(
                     reactiveVals$sce
                   )), multiple = F),
    width = 4
  )
})

output$designMds <- renderUI({
  shinydashboard::box(
    selectizeInput(
      "mdsFeatures",
      "Features::",
      unique(rowData(reactiveVals$sceT)$marker_class),
      multiple = F
    ),
    selectizeInput("mdsColorBy",
                   "Color by:",
                   names(colData(
                     reactiveVals$sce
                   )), multiple = F),
    selectizeInput("mdsLabels",
                   "Label by:",
                   unique(colData(
                     reactiveVals$sce
                   )), multiple = F),
    width = 4
  )
})

output$designExprs <- renderUI({
  shinydashboard::box(
    selectizeInput(
      "exprFeatures",
      "Features:",
      unique(rowData(reactiveVals$sce)$marker_class),
      multiple = F
    ),
    selectizeInput("exprsColorBy",
                   "Color by:",
                   names(colData(
                     reactiveVals$sce
                   )), multiple = F),
    width = 4
  )
})

output$designHeatmapExprs <- renderUI({
  shinydashboard::box(
    selectizeInput(
      "hpexprFeatures",
      "Features:",
      unique(rowData(reactiveVals$sce)$marker_class),
      multiple = F
    ),
    selectizeInput(
      "hpexprScale",
      "Scale:",
      c("first", "last", "never"),
      multiple = F
    ),
    width = 4
  )
  
})
