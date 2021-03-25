# Preprocessing Server

# transform SingleCellExperiment
transformData <-
  function (sce,
            cf = 5,
            ain = "counts",
            aout = "exprs") {
    y <- assay(sce, ain)
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
    assay(sce, aout, FALSE) <- y
    sce
  }

# marker, sample, patient selection -> if markers, patients, samples are selected -> prepValButton can be clicked
observeEvent({
  input$patientSelection
  input$sampleSelection
}, {
  if (length(input$sampleSelection) == 0 || length(input$patientSelection)==0){
    reactiveVals$continue <- TRUE
    shinyjs::disable("prepSelectionButton")
    shinyjs::disable("filterSelectionButton")
  } else {
    shinyjs::enable("prepSelectionButton")
    shinyjs::enable("filterSelectionButton")
  }
}, ignoreNULL = FALSE)

# check current tab
observe({
  if (reactiveVals$current_tab == 3) {
    if (!reactiveVals$preprocessingShowed){
      plotPreprocessing(reactiveVals$sce)
      reactiveVals$preprocessingShowed <- TRUE
    }
    if (!("patient_id" %in% colnames(colData(reactiveVals$sce)))) {
      shinyjs::hide("patientsBox")
    }
  }
})

# render markers box
output$markersBox <- renderUI({
  pickerInput(
    inputId = "markerSelection",
    label = "Markers",
    choices = names(channels(reactiveVals$sce)),
    selected = names(channels(reactiveVals$sce)),
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  )
})

# render samples box
output$samplesBox <- renderUI({
  pickerInput(
    "sampleSelection",
    choices = as.character(unique(colData(reactiveVals$sce)$sample_id)),
    selected = as.character(unique(colData(reactiveVals$sce)$sample_id)),
    label = "Samples",
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  )
})

# render samples box
output$patientsBox <- renderUI({
  pickerInput(
    "patientSelection",
    choices = as.character(unique(colData(reactiveVals$sce)$patient_id)),
    selected = as.character(unique(colData(reactiveVals$sce)$patient_id)),
    label = "Patients",
    options = list(
      `actions-box` = TRUE,
      size = 4,
      `selected-text-format` = "count > 3",
      "dropup-auto" = FALSE
    ),
    multiple = TRUE
  )
})

# render sortable lists for all conditions
output$reorderingTabs <- renderUI({
  library(sortable)
  conditions <- names(metadata(reactiveVals$sce)$experiment_info)
  conditions <- conditions[!conditions %in% c("sample_id", "patient_id", "n_cells")]
  lapply(conditions, function(condition){
    rank_list(
      text = paste0("Reorder the condition: ",condition),
      labels = levels(metadata(reactiveVals$sce)$experiment_info[[condition]]),
      input_id = condition
    )
  })
  
})

# if conditions are ordered
observeEvent(input$reorderButton, {
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("filterSelectionButton")
  shinyjs::disable("previousTab")
  shinyjs::disable("nextTab")
  
  # reorder levels 
  conditions <- names(metadata(reactiveVals$sce)$experiment_info)
  conditions <- conditions[!conditions %in% c("sample_id", "patient_id", "n_cells")]
  lapply(conditions, function(condition){
    ordered <- input[[condition]]
    levels(reactiveVals$sce[[condition]]) <- ordered
  })
  plotPreprocessing(reactiveVals$sce)
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("filterSelectionButton")
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")

})

# if start transformation button is clicked
observeEvent(input$prepButton, {
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("filterSelectionButton")
  shinyjs::disable("previousTab")
  shinyjs::disable("nextTab")
  # data transformation
  reactiveVals$sce <-
    transformData(sce = reactiveVals$sce,
                  cf = as.numeric(input$cofactor))
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("filterSelectionButton")
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
})

# if visualize selection button is clicked
observeEvent(input$prepSelectionButton, {
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("previousTab")
  shinyjs::disable("nextTab")
  shinyjs::disable("filterSelectionButton")
  allpatients <- length(as.character(unique(colData(reactiveVals$sce)$patient_id)))
  allsamples <- length(as.character(unique(colData(reactiveVals$sce)$sample_id)))
  if ((length(input$patientSelection) != allpatients) || (length(input$sampleSelection) != allsamples)){
    showNotification(HTML(
      "<b>Attention!</b><br>
      The unselected samples and patients are <b>deleted</b> from the data when pressing the <b>Confirm Selection</b> button. Further analysis is being performed only on the selected patients and samples!"
    ),
    duration = 10,
    type = "warning")
  }
  markers <- isolate(input$markerSelection)
  samples <- isolate(input$sampleSelection)
  patients <- isolate(input$patientSelection)
  sce <- filterSCE(reactiveVals$sce,sample_id %in% samples)
  if (("patient_id" %in% colnames(colData(reactiveVals$sce)))){
    sce <- filterSCE(sce,patient_id %in% patients)
  }
  
  sce <- sce[rownames(sce) %in% markers, ]
  plotPreprocessing(sce)
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("continue")
  shinyjs::enable("filterSelectionButton")
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
})

# if filtering button is clicked -> selection is applied to sce
observeEvent(input$filterSelectionButton,{
  shinyjs::disable("prepButton")
  shinyjs::disable("prepSelectionButton")
  shinyjs::disable("filterSelectionButton")
  shinyjs::disable("previousTab")
  shinyjs::disable("nextTab")
  allpatients <- length(as.character(unique(colData(reactiveVals$sce)$patient_id)))
  allsamples <- length(as.character(unique(colData(reactiveVals$sce)$sample_id)))
  
  if (length(input$sampleSelection) != allsamples){
    reactiveVals$sce <- filterSCE(reactiveVals$sce,sample_id %in% input$sampleSelection)
  }
  if (("patient_id" %in% colnames(colData(reactiveVals$sce)))){
    if (length(input$patientSelection) != allpatients){
      reactiveVals$sce <- filterSCE(reactiveVals$sce,patient_id %in% input$patientSelection)
    }
  }
  
  markers <- isolate(input$markerSelection)
  sce <- reactiveVals$sce[rownames(reactiveVals$sce) %in% markers, ]
  plotPreprocessing(sce)
  
  shinyjs::enable("prepButton")
  shinyjs::enable("prepSelectionButton")
  shinyjs::enable("filterSelectionButton")
  shinyjs::enable("previousTab")
  shinyjs::enable("nextTab")
})

observeEvent(reactiveVals$sce, {
  if ("exprs" %in% names(assays(reactiveVals$sce)))
    shinyjs::hide("noTransformationWarning")
  else 
    shinyjs::show("noTransformationWarning")
})

# method for plotting all kinds of preprocessing plots
plotPreprocessing <- function(sce) {
  groupColorLabelBy <- names(colData(sce))
  possAssays <- assayNames(sce)
  if (all(possAssays == c("counts", "exprs"))) {
    possAssays <- c("Normalized" = "exprs", "Raw" = "counts")
  }
  features <-
    c("all", as.character(unique(rowData(sce)$marker_class)))
  
  ## COUNTS
  
  # ui for counts
  output$designCounts <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("countsGroupBy",
                       "Group by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput("countsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "countsProp",
          "Stacked or dodged:",
          c(
            "dodged (total cell counts)" = FALSE,
            "stacked (relative abundance)" = TRUE
          ),
          multiple = F
        ),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),

      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("countsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("countsPlotDownload"),
      style = "position: absolute; bottom: 10px; right:10px;"
    ))
  })
  
  # render counts plot
  output$countsPlot <- renderPlot({
    reactiveVals$countsPlot <-  CATALYST::plotCounts(
      sce,
      group_by = input$countsGroupBy,
      color_by = input$countsColorBy,
      prop = as.logical(input$countsProp)
    )
    reactiveVals$countsPlot
  })
  
  # ui for download button
  output$countsPlotDownload <- renderUI({
    req(reactiveVals$countsPlot)
    downloadButton("downloadPlotCounts", "Download Plot")
  })
  
  # function for downloading count plot
  output$downloadPlotCounts <- downloadHandler(
    filename = function(){
      paste0("Counts_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$countsPlot, width=12, height=6)
    }
  )
  
  ## MDS 
  
  # ui for MDS
  output$designMDS <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("mdsLabelBy",
                       "Label by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput("mdsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "mdsAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput("mdsFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      ),
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("mdsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("mdsPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px"
    ),)
  })
  
  # render mds plot
  output$mdsPlot <- renderPlot({
    feature <- input$mdsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$mdsPlot <- CATALYST::pbMDS(
      sce,
      label_by = input$mdsLabelBy,
      color_by = input$mdsColorBy,
      features = feature,
      assay = input$mdsAssay,
    )
    reactiveVals$mdsPlot
    
  })
  
  # ui for download button
  output$mdsPlotDownload <- renderUI({
    req(reactiveVals$mdsPlot)
    downloadButton("downloadPlotMDS", "Download Plot")
  })
  
  # function for downloading MDS plot
  output$downloadPlotMDS <- downloadHandler(
    filename = function(){
      paste0("MDS_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$mdsPlot, width=16, height=11)
    }
  )
  
  ## NRS
  
  # ui for NRS
  output$designNRS <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("nrsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "nrsAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput("nrsFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("nrsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("nrsPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px;"
    ))
  })
  
  # render nrs plot
  output$nrsPlot <- renderPlot({
    feature <- input$nrsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$nrsPlot <- CATALYST::plotNRS(
      sce,
      color_by = input$nrsColorBy,
      features = feature,
      assay = input$nrsAssay
    )
    reactiveVals$nrsPlot
  })
  
  # ui for download button
  output$nrsPlotDownload <- renderUI({
    req(reactiveVals$nrsPlot)
    downloadButton("downloadPlotNRS", "Download Plot")
  })
  
  # function for downloading NRS plot
  output$downloadPlotNRS <- downloadHandler(
    filename = function(){
      paste0("NRS_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$nrsPlot, width=12, height=6)
    }
  )
  
  ## Exprs

  # ui for expr
  output$designExprs <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput("exprsColorBy",
                       "Color by:",
                       groupColorLabelBy, multiple = F),
        selectizeInput(
          "exprsAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput("exprsFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("exprsPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("exprsPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px;"
    ),)
  })
  
  # render exprs plot
  output$exprsPlot <- renderPlot({
    feature <- input$exprsFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$exprsPlot <- CATALYST::plotExprs(
      sce,
      color_by = input$exprsColorBy,
      features = feature,
      assay = input$exprsAssay
    )
    reactiveVals$exprsPlot
  })
  
  # ui for download button
  output$exprsPlotDownload <- renderUI({
    req(reactiveVals$exprsPlot)
    downloadButton("downloadPlotExprs", "Download Plot")
  })
  
  # function for downloading exprs plot
  output$downloadPlotExprs <- downloadHandler(
    filename = function(){
      paste0("Expr_Plot", ".pdf")
    },
    content = function(file){
      ggsave(file, plot = reactiveVals$exprsPlot, width=14, height=9)
    }
  )
  

  ## Exprs Heatmap
  
  # ui for exprs heatmap
  output$designExprsHeatmap <- renderUI({
    fluidRow(column(
      1,
      div(dropdownButton(
        tags$h3("Plot Options"),
        selectizeInput(
          "exprsHeatmapScale",
          "Scale:",
          c("never", "first", "last"),
          multiple = F
        ),
        selectizeInput(
          "exprsHeatmapAssay",
          "Raw or normalized counts:",
          possAssays,
          multiple = F
        ),
        selectizeInput("exprsHeatmapFeatures",
                       "Features:",
                       features,
                       multiple = F),
        circle = TRUE,
        status = "info",
        icon = icon("gear"),
        width = "400px",
        tooltip = tooltipOptions(title = "Click to see plot options")
      ),
      style = "position: relative; height: 500px;"
      )
    ),
    column(11, shinycssloaders::withSpinner(
      plotOutput("exprsHeatmapPlot", width = "100%", height = "500px")
    )),
    div(
      uiOutput("exprsHeatmapPlotDownload"),
      style = "position: absolute; bottom: 10px;right:10px;"
    ),)
  })
  
  # render exprs heatmap plot
  output$exprsHeatmapPlot <- renderPlot({
    feature <- input$exprsHeatmapFeatures
    if (feature == "all") {
      feature <- NULL
    }
    reactiveVals$exprsPlotHeatmap <- plotExprHeatmapCustom(
      sce,
      scale = input$exprsHeatmapScale,
      features = feature,
      assay = input$exprsHeatmapAssay
    )
    reactiveVals$exprsPlotHeatmap
  })
  
  # ui for download button
  output$exprsHeatmapPlotDownload <- renderUI({
    req(reactiveVals$exprsPlotHeatmap)
    library(ComplexHeatmap)
    downloadButton("downloadPlotExprsHeatmap", "Download Plot")
  })
  
  # function for downloading exprs heatmap
  output$downloadPlotExprsHeatmap <- downloadHandler(
    filename = "Expression_Heatmap.pdf", 
    content = function(file){
      pdf(file, width = 12, height = 8)
      draw(reactiveVals$exprsPlotHeatmap)
      dev.off()
    }
  )
  
}

# function anno_features from CATALYST with different color palette
.anno_factors <- function (x, ids, which, type = c("row", "column")) 
{
  type <- match.arg(type)
  cd <- colData(x)
  df <- data.frame(cd, check.names = FALSE)
  df <- select_if(df, ~!is.numeric(.))
  df <- mutate_all(df, ~droplevels(factor(.x)))
  m <- match(ids, df$sample_id)
  ns <- split(df, df$sample_id) %>% lapply(mutate_all, droplevels) %>% 
    lapply(summarize_all, nlevels) %>% do.call(what = "rbind")
  keep <- names(which(colMeans(ns) == 1))
  keep <- setdiff(keep, c("sample_id", "cluster_id"))
  if (is.character(which)) 
    keep <- intersect(keep, which)
  if (length(keep) == 0) 
    return(NULL)
  df <- df[m, keep, drop = FALSE]
  lvls <- lapply(as.list(df), levels)
  nlvls <- vapply(lvls, length, numeric(1))
  #pal <- brewer.pal(8, "Set2")
  pal <-  c("#999999","#009E73","#E69F00", "#56B4E9", "#F0E442","#D55E00","#0072B2", "#CC79A7")
  names(is) <- is <- colnames(df)
  cols <- lapply(is, function(i) {
    if (nlvls[i] > length(pal)) 
      pal_i <- colorRampPalette(pal)(max(nlvls))
    else pal_i <- pal
    u <- pal_i[seq_len(nlvls[i])]
    names(u) <- lvls[[i]]
    u
  })
  ComplexHeatmap::HeatmapAnnotation(which = type, df = df, col = cols, gp = gpar(col = "white"))
}

plotExprHeatmapCustom <- function (x, features = NULL, by = c("sample_id", "cluster_id", 
                                                        "both"), k = "meta20", m = NULL, assay = "exprs", fun = c("median", 
                                                                                                                  "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, 
                             row_anno = TRUE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
                             row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, 
                             bin_anno = FALSE, hm_pal = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), 
                             k_pal = CATALYST:::.cluster_cols, m_pal = k_pal, distance = c("euclidean", 
                                                                                           "maximum", "manhattan", "canberra", "binary", "minkowski"), 
                             linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
                                         "median", "centroid", "ward.D2")) 
{
  args <- as.list(environment())
  CATALYST:::.check_args_plotExprHeatmap(args)
  distance <- match.arg(distance)
  linkage <- match.arg(linkage)
  scale <- match.arg(scale)
  fun <- match.arg(fun)
  by <- match.arg(by)
  x <- x[unique(CATALYST:::.get_features(x, features)), ]
  if (by != "sample_id") {
    CATALYST:::.check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
  }
  if (by == "both") 
    by <- c("cluster_id", "sample_id")
  .do_agg <- function() {
    z <- CATALYST:::.agg(x, by, fun, assay)
    if (length(by) == 1) 
      return(z)
    magrittr::set_rownames(do.call("rbind", z), levels(x$cluster_id))
  }
  .do_scale <- function() {
    if (scale == "first") {
      z <- assay(x, assay)
      z <- CATALYST:::.scale_exprs(z, 1, q)
      assay(x, assay, FALSE) <- z
      return(x)
    }
    else CATALYST:::.scale_exprs(z, 1, q)
  }
  z <- switch(scale, first = {
    x <- .do_scale()
    .do_agg()
  }, last = {
    z <- .do_agg()
    .do_scale()
  }, never = {
    .do_agg()
  })
  if (length(by) == 1) 
    z <- t(z)
  if (scale != "never" && !(assay == "counts" && fun == "sum")) {
    qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  }
  else lgd_aes <- list()
  lgd_aes$title_gp <- grid::gpar(fontsize = 10, fontface = "bold", 
                                 lineheight = 0.8)
  if (!isFALSE(row_anno)) {
    left_anno <- switch(by[1], sample_id = .anno_factors(x, 
                                                         levels(x$sample_id), row_anno, "row"), CATALYST:::.anno_clusters(x, 
                                                                                                               k, m, k_pal, m_pal))
  }
  else left_anno <- NULL
  if (!isFALSE(col_anno) && length(by) == 2) {
    top_anno <- .anno_factors(x, levels(x$sample_id), col_anno, 
                              "colum")
  }
  else top_anno <- NULL
  if (bars) {
    right_anno <- CATALYST:::.anno_counts(x[[by[1]]], perc)
  }
  else right_anno <- NULL
  if (bin_anno) {
    cell_fun <- function(j, i, x, y, ...) grid.text(gp = grid::gpar(fontsize = 8), 
                                                    sprintf("%.2f", z[i, j]), x, y)
  }
  else cell_fun <- NULL
  a <- ifelse(assay == "exprs", "expression", assay)
  f <- switch(fun, median = "med", fun)
  hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
                                            "scaled", a), last = sprintf("%s %s\n%s", "scaled", 
                                                                         fun, a), never = paste(fun, a, sep = "\n"))
  if (length(by) == 2) {
    col_title <- features
  }
  else if (length(features) == 1 && features %in% c("type", 
                                                    "state")) {
    col_title <- paste0(features, "_markers")
  }
  else col_title <- ""
  ComplexHeatmap::Heatmap(matrix = z, name = hm_title, col = circlize::colorRamp2(seq(min(z), 
                                                                                      max(z), l = n <- 100), grDevices::colorRampPalette(hm_pal)(n)), 
                          column_title = col_title, column_title_side = ifelse(length(by) == 
                                                                                 2, "top", "bottom"), cell_fun = cell_fun, cluster_rows = row_clust, 
                          cluster_columns = col_clust, show_row_dend = row_dend, 
                          show_column_dend = col_dend, clustering_distance_rows = distance, 
                          clustering_method_rows = linkage, clustering_distance_columns = distance, 
                          clustering_method_columns = linkage, show_row_names = (is.null(left_anno) || 
                                                                                   isTRUE(by == "sample_id")) && !perc, row_names_side = ifelse(by[1] == 
                                                                                                                                                  "cluster_id" || isFALSE(row_anno) && !row_dend || 
                                                                                                                                                  isFALSE(row_clust), "left", "right"), top_annotation = top_anno, 
                          left_annotation = left_anno, right_annotation = right_anno, 
                          rect_gp = gpar(col = "white"), heatmap_legend_param = lgd_aes)
}