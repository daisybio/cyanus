reactiveVals$data <- list(upload = list(fcs=NULL, panel=NULL, md=NULL), example = list(fcs=NULL, panel=NULL, md=NULL) )

##### Observers #####

observeEvent(input$fcsFiles, {
  setup_logfile()
  fileTable <- input$fcsFiles
  fileTable <- fileTable[, c("name", "size")]
  fileTable$size <- sprintf("%.2f MB", fileTable$size / 1000000)
  reactiveVals$data$upload$fcs <- fileTable
})

observeEvent(input$metaFile, {
  if(endsWith(tolower(input$metaFile$datapath), ".csv")){
    reactiveVals$data$upload$md <- data.table::fread(input$metaFile$datapath)
    data.table::setDF(reactiveVals$data$upload$md)
  }else if(endsWith(tolower(input$metaFile$datapath), ".xls") | 
           endsWith(tolower(input$metaFile$datapath), ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files. If you can, please upload a .csv file", type = "warning")
    reactiveVals$data$upload$md <- read.xlsx2(input$metaFile$datapath, 1)
  } else {
    showNotification(HTML("<b>Unknown file extension.</b><br>Please upload a CSV (*.csv) or Excel (*.xls, *.xlsx) file.<br>Example: Check out the PBMC Example Data."), duration = NULL, type = "error")
  }
})

observeEvent(input$panelFile, {
  if(endsWith(tolower(input$panelFile$datapath), ".csv")){
    tmp_panel <-
      data.table::fread(input$panelFile$datapath)
    data.table::setDF(tmp_panel)
  }else if(endsWith(tolower(input$panelFile$datapath), ".xls") | 
           endsWith(tolower(input$panelFile$datapath), ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    tmp_panel <- read.xlsx2(input$panelFile$datapath, 1)
  } else {
    showNotification(HTML("<b>Unknown file extension.</b><br>Please upload a CSV (*.csv) or Excel (*.xls, *.xlsx) file.<br>Example: Check out the PBMC Example Data."), duration = NULL, type = "error")
    return()
  }
  if(any(!c("fcs_colname", "antigen") %in% names(tmp_panel))){
    # try setting header to TRUE
    tmp_panel <- data.table::fread(input$panelFile$datapath, header = TRUE)
    data.table::setDF(tmp_panel)
    if(any(!c("fcs_colname", "antigen") %in% names(tmp_panel))){
      showNotification(HTML("Error while reading the panel file:<br>A CSV or Excel file with headers describing the panel:<br>for each channel:<br>fcs_colname: its column name in the input data<br>antigen: targeted protein marker<br>marker_class: (optionally) class (type, state, or none)<br>i.e.:<br>fcs_colname,antigen[,marker_class]<br><b>Example: Check out the PBMC Example Data</b>"), duration = NULL, type = "error")
    }else{
      reactiveVals$data$upload$panel <- tmp_panel[, colnames(tmp_panel) %in% c("fcs_colname", "antigen", "marker_class")]
    }
  }else{
    reactiveVals$data$upload$panel <- tmp_panel[, colnames(tmp_panel) %in% c("fcs_colname", "antigen", "marker_class")]
  }
})

observeEvent(input$sceFile, {
  setup_logfile()
  library(CATALYST)
  tmp <- readRDS(file.path(input$sceFile$datapath))
  if (class(tmp) == "SingleCellExperiment"){
    reactiveVals$data$sce$rowdata <- rowData(tmp)
    reactiveVals$data$sce$coldata <- metadata(tmp)$experiment_info
  } else {
    showNotification("You have to upload a SingleCellExperiment object!", type="error")
    reset("sceFile")
  }
}, ignoreInit = TRUE)

observeEvent(input$exampleData, {
  setup_logfile()
  reactiveVals$data$example$fcs <- readRDS(file.path(input$exampleData, "fcs.rds"))
  reactiveVals$data$example$panel <- readRDS(file.path(input$exampleData, "panel.rds"))
  reactiveVals$data$example$md <- readRDS(file.path(input$exampleData, "md.rds"))
}, ignoreInit = TRUE)

observeEvent(input$loadData, {
  updateButton(session, "loadData", label = " Loading...", disabled = TRUE)
  resetPreprocessing()
  resetVisualization()
  resetClustering()
  resetDE()
  resetDEComparison()
  waiter_show(id = "app",html = tagList(spinner$logo, 
                                        HTML("<br>Loading Data...<br>Please be patient")), 
              color=spinner$color)
  
  library(CATALYST)
  reactiveVals$has_error <- FALSE
  if (input$chooseDataTab == "dataUpload") {
    reactiveVals$call_list <<- append(reactiveVals$call_list, quote(dn <- '~/path/to/fcs_data.fcs'))
    dn <- dirname(input$fcsFiles$datapath)[1]
    file.rename(input$fcsFiles$datapath, file.path(dn, "/", input$fcsFiles$name))
    
    conditions <- names(reactiveVals$data$upload$md)[!names(reactiveVals$data$upload$md) 
                                                                        %in% c("sample_id", "file_name")]
    md_cols <- list(file = "file_name", id = "sample_id", factors = conditions)
    
    tryCatch({
      withCallingHandlers({
      reactiveVals$sce <- evap(expression(reactiveVals$sce <- CATALYST::prepData(
        dn,
        panel = panel,
        md = md,
        transform = FALSE,
        md_cols = md_cols,
        FACS = FACS
      ))[[1]], params = list(panel = reactiveVals$data$upload$panel,
                             md = reactiveVals$data$upload$md,
                             md_cols = md_cols,
                             FACS = input$isFACSData))},
      message = function(m) {
        showNotification(HTML(sprintf("Loading the data produced with the following message:<br>
                                    <b>%s</b>", m$message)), duration = NULL, type = "message")
      },
      warning = function(w) {
        showNotification(HTML(sprintf("Loading the data produced with the following warning:<br>
                                    <b>%s</b>", w$message)), duration = NULL, type = "warning")
      }
    )},
    error = function(e){
      reactiveVals$has_error <- TRUE
      showNotification(HTML(sprintf("Loading the data failed with the following message:<br>
                                    <b>%s</b>", e$message)), duration = NULL, type = "error")
      start_tab <- which(tab_ids == "start")
      reactiveVals$continue[start_tab] <- FALSE
      if (isolate(reactiveVals$max_tab) > start_tab)
        reactiveVals$max_tab <- start_tab
    })
  } else if (input$chooseDataTab == "dataExample") {
    reactiveVals$call_list <<- append(reactiveVals$call_list, quote(reactiveVals$sce <- readRDS('~/path/to/example_data/sce.rds')))
    reactiveVals$sce <- readRDS(file.path(input$exampleData, "sce.rds"))
  } else if (input$chooseDataTab == "sceUpload") {
    reactiveVals$call_list <<- append(reactiveVals$call_list, quote(reactiveVals$sce <- readRDS('~/path/to/sce.rds')))
    reactiveVals$sce <- readRDS(file.path(input$sceFile$datapath))
    
  }else
    stop("Which tab is selected?")
  
  if (!is.null(reactiveVals$sce) & !reactiveVals$has_error) { # meaning the data loading worked
  
  #drop levels of markers
    SummarizedExperiment::rowData(reactiveVals$sce)$marker_class <- evap(expression(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class <- droplevels(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class))[[1]],
                                                                         params = NULL)
  #SummarizedExperiment::rowData(reactiveVals$sce)$marker_class <- droplevels(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
  # set negative values to zero
  if (length(which(assays(reactiveVals$sce)$counts < 0)) > 0){
    assays(reactiveVals$sce)$counts[assays(reactiveVals$sce)$counts < 0 ] <- evap(expression(assays(reactiveVals$sce)$counts[assays(reactiveVals$sce)$counts < 0 ] <- 0)[[1]])
    if ("exprs" %in% assayNames(reactiveVals$sce)){
      assays(reactiveVals$sce)$exprs[assays(reactiveVals$sce)$exprs < 0 ] <- evap(expression(assays(reactiveVals$sce)$exprs[assays(reactiveVals$sce)$exprs < 0 ] <- 0)[[1]])
    }
    showNotification(
      HTML(
        "<b>Negative values found.</b><br>
      Negative values were found in the data and now set to 0."
      ),
      duration = 10,
      type = "message"
    )
  }
  if(input$chooseDataTab == "dataUpload" & is.null(reactiveVals$data$upload$panel)){
    # create and show
    tmp_panel <- as.data.table(rowData(reactiveVals$sce))
    tmp_panel <- tmp_panel[, c('channel_name', 'marker_name', 'marker_class')]
    colnames(tmp_panel) <- c('fcs_colname', 'antigen', 'marker_class')
    reactiveVals$data$upload$panel <- tmp_panel
  }
  if(input$chooseDataTab == "dataUpload" & is.null(reactiveVals$data$upload$md)){
    # create and show
    file_names <- unlist(unname(input$fcsFiles['name']))
    sample_ids  <- tstrsplit(file_names, '.fcs', keep = 1)[[1]]
    reactiveVals$data$upload$md <- data.table('file_name' = file_names, 'sample_id' = sample_ids)
  }
  
  start_tab <- which(tab_ids == "start")
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
  reactiveVals$continue[start_tab] <- TRUE
  runjs("document.getElementById('nextTab').scrollIntoView();")
    
  }
  updateButton(session, "loadData", label = " Load Data", disabled = FALSE) 
  waiter_hide(id = "app")
  })

observeEvent(reactiveVals$data$upload$panel, {
  req(reactiveVals$data$upload$panel)
  if(!is.null(reactiveVals$data$upload$panel)){
  }
  if(any(duplicated(reactiveVals$data$upload$panel$fcs_colname))){
    showNotification(paste("Some fcs_colnames are duplicated:", 
                           reactiveVals$data$upload$panel$fcs_colname[anyDuplicated(reactiveVals$data$upload$panel$fcs_colname)], 
                           ". This can cause problems, please change it!"), 
                     type = "warning", 
                     duration = NULL)
  }
  if(any(duplicated(reactiveVals$data$upload$panel$antigen))){
    showNotification(paste("Some antigens are duplicated:", 
                           reactiveVals$data$upload$panel$antigen[anyDuplicated(reactiveVals$data$upload$panel$antigen)], 
                           ". This can cause problems, please change it!"), 
                     type = "warning", 
                     duration = NULL)
  }
})

observeEvent(input$panelDT_cell_edit, {
  start_tab <- which(tab_ids == "start")
  reactiveVals$continue[start_tab] <- FALSE
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
  reactiveVals$data$upload$panel <<- editData(reactiveVals$data$upload$panel, input$panelDT_cell_edit, "panelDT")
})


observeEvent(input$metadataDT_cell_edit, {
  start_tab <- which(tab_ids == "start")
  reactiveVals$continue[start_tab] <- FALSE
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
  reactiveVals$data$upload$md <<- editData(reactiveVals$data$upload$md, input$metadataDT_cell_edit, "metadataDT")
})


observeEvent(input$panel_add_row, {
  req(reactiveVals$data$upload$panel)
  new_row <- data.table(t(rep('*E D I T*', ncol(reactiveVals$data$upload$panel))))
  colnames(new_row) <- colnames(reactiveVals$data$upload$panel)
  new_row[, (colnames(new_row)) := lapply(.SD, as.factor), .SDcols = colnames(new_row)]
  reactiveVals$data$upload$panel[, (colnames(reactiveVals$data$upload$panel)) := lapply(.SD, as.factor), .SDcols = colnames(reactiveVals$data$upload$panel)]
  reactiveVals$data$upload$panel <- rbindlist(list(new_row,
                                                   reactiveVals$data$upload$panel))
  start_tab <- which(tab_ids == "start")
  reactiveVals$continue[start_tab] <- FALSE
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
}
)

observeEvent(input$panel_add_column, {
  req(reactiveVals$data$upload$panel)
  reactiveVals$data$upload$panel$marker_class <- 'none'
  start_tab <- which(tab_ids == "start")
  reactiveVals$continue[start_tab] <- FALSE
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
})

observeEvent(input$md_add_column, {
  req(reactiveVals$data$upload$md)
  feature_name <- paste0('condition_', ncol(reactiveVals$data$upload$md) - 1)
  reactiveVals$data$upload$md[, eval(feature_name)] <- '* E D I T *'
  start_tab <- which(tab_ids == "start")
  reactiveVals$continue[start_tab] <- FALSE
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
})


observeEvent(input$panel_delete_rows,{
  if (!is.null(input$panelDT_rows_selected)) {
    reactiveVals$data$upload$panel <- reactiveVals$data$upload$panel[-as.numeric(input$panelDT_rows_selected), ]
    start_tab <- which(tab_ids == "start")
    reactiveVals$continue[start_tab] <- FALSE
    if (isolate(reactiveVals$max_tab) > start_tab)
      reactiveVals$max_tab <- start_tab
  }
})

observeEvent(input$md_remove_column, {
  if (!is.null(input$metadataDT_columns_selected)) {
    to_delete <- c(as.numeric(input$metadataDT_columns_selected))
    if(!any(c('file_name', 'sample_id') %in% colnames(reactiveVals$data$upload$md)[to_delete])){
      reactiveVals$data$upload$md <-reactiveVals$data$upload$md[, colnames(reactiveVals$data$upload$md)[-to_delete], with = FALSE]
      if(ncol(reactiveVals$data$upload$md) > 2){
        colnames(reactiveVals$data$upload$md) <- c('file_name', 'sample_id', paste0('condition_', 
                                                                                    seq(ncol(reactiveVals$data$upload$md)-2)))
      }else if(ncol(reactiveVals$data$upload$md) == 2){
        colnames(reactiveVals$data$upload$md) <- c('file_name', 'sample_id')
      }
      start_tab <- which(tab_ids == "start")
      reactiveVals$continue[start_tab] <- FALSE
      if (isolate(reactiveVals$max_tab) > start_tab)
        reactiveVals$max_tab <- start_tab
    }else{
      showNotification('file_name and sample_id cannot be deleted!', type = 'error')
    }
  }
})

observeEvent(input$panel_reset, {
  req(reactiveVals$data$upload$panel)
  reactiveVals$data$upload$panel <- NULL
  start_tab <- which(tab_ids == "start")
  reactiveVals$continue[start_tab] <- FALSE
  if (isolate(reactiveVals$max_tab) > start_tab)
    reactiveVals$max_tab <- start_tab
})

observeEvent(input$chosenColorPalette, {
  if(input$chosenColorPalette == "Colorblind1"){
    reactiveVals$selected_palette <- c("#ff6db6", "#004949", "#db6d00",  "#B2DF8A", 
                                         "#FDB462", "#490092", "#009999", "#8f4e00", 
                                         "#ffdf4d", "#171723","#b66dff")
  }else if(input$chosenColorPalette == "Colorblind2: Wong"){
    reactiveVals$selected_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }else if(input$chosenColorPalette == "Colorblind3: Tol bright"){
    reactiveVals$selected_palette <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', 
                                         '#66CCEE', '#AA3377', '#BBBBBB')
  }else if(input$chosenColorPalette == "Colorblind4: Tol vibrant"){
    reactiveVals$selected_palette <- c('#EE7733', '#0077BB', '#33BBEE', '#EE3377', 
                                         '#CC3311', '#009988', '#BBBBBB')
  }else if(input$chosenColorPalette == "Colorblind5: Tol muted"){
    reactiveVals$selected_palette <- c('#CC6677', '#332288', '#DDCC77', '#117733', 
                                         '#88CCEE', '#882255', '#44AA99', '#999933', 
                                         '#AA4499', '#DDDDDD')
  }else if(input$chosenColorPalette  == "RColorBrewer: Set1"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(9, "Set1")
  }else if(input$chosenColorPalette  == "RColorBrewer: Set2"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(8, "Set2")
  }else if(input$chosenColorPalette  == "RColorBrewer: Set3"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(12, "Set3")
  }else if(input$chosenColorPalette  == "RColorBrewer: Pastel1"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(9, "Pastel1")
  }else if(input$chosenColorPalette  == "RColorBrewer: Pastel2"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(8, "Pastel2")
  }else if(input$chosenColorPalette  == "RColorBrewer: Paired"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(12, "Paired")
  }else if(input$chosenColorPalette  == "RColorBrewer: Dark2"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(8, "Dark2")
  }else if(input$chosenColorPalette  == "RColorBrewer: Accent"){
    reactiveVals$selected_palette <- RColorBrewer::brewer.pal(8, "Accent")
  }else if(input$chosenColorPalette  == "CATALYST colors"){
    reactiveVals$selected_palette <- CATALYST:::.cluster_cols
  }
})


observeEvent(input$addColor, {
  color <- input$colorInput
  # Check if input is a valid color code
  if (grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", color)) {
    reactiveVals$last_colors <- c(reactiveVals$last_colors, color)
  } else {
    # Try to convert color name to hexadecimal code
    hex <- tryCatch(grDevices::col2rgb(color, alpha = FALSE), error = function(e) NULL)
    if (!is.null(hex)) {
      reactiveVals$last_colors <- c(reactiveVals$last_colors, rgb(hex[1], hex[2], hex[3], maxColorValue = 255))
    } else {
      # If not a valid color, don't add it
      showNotification(HTML("Not a valid color"), duration = NULL, type = "error")
      return(NULL)
    }
  }
  # Clear the input field after adding the color
  updateTextInput(session, "colorInput", value = "")
})

observeEvent(input$clearColors, {
  reactiveVals$last_colors <- NULL
})

observeEvent(input$setColorPalette, {
  if(is.null(reactiveVals$last_colors)){
    showNotification(HTML("You did not set a color palette!"), type="error")
  }else{
    reactiveVals$selected_palette <- reactiveVals$last_colors
  }
})


#### Output elements ####

output$panelDT <- renderDT(
  checkNullTable(reactiveVals$data$upload$panel),
  editable = "cell"
)

output$metadataDT <- renderDT(
  checkNullTable(reactiveVals$data$upload$md),
  editable = "cell",
  selection = list(target = 'column')
  )

output$currentData <- renderInfoBox({
  if(input$chooseDataTab == "dataUpload"){
    fcs <- reactiveVals$data$upload$fcs
    panel <- reactiveVals$data$upload$panel
    md <- reactiveVals$data$upload$md
  } else if (input$chooseDataTab == "sceUpload"){
    rowdata <- reactiveVals$data$sce$rowdata
    coldata <- reactiveVals$data$sce$coldata
  }else {
    fcs <- reactiveVals$data$example$fcs
    panel <- reactiveVals$data$example$panel
    md <- reactiveVals$data$example$md
  }
  
  status <- "warning"

  # if current data tab is sce upload
  if (input$chooseDataTab == "sceUpload"){
    value <- list(
      renderTable(
        checkNullTable(rowdata),
        caption = "Rowdata",
        caption.placement = "top"
      ), 
      renderTable(
        checkNullTable(coldata),
        caption = "Coldata",
        caption.placement = "top"
      )
    )
  
  } else { 
  # if current data tab is upload data
  if(input$chooseDataTab == "dataUpload"){
    tablePanel <- fluidRow(
      column(12, h5(HTML("<p style=\"color:grey;\">Panel Data</p>"))),
      column(
      12, h5(HTML("<b style=\"color:grey;\">If you provide a panel, you can edit its cells by double-clicking. If you want to edit the automatically guessed panel, just upload your fcs files and click on 'Load Data'. You can then modify the guessed panel. </b>")), hr(), DTOutput("panelDT")),
      column(12, 
             bsButton("panel_add_column", "Add Marker Class Column", disabled = ifelse(is.null(panel), TRUE, ifelse('marker_class' %in% colnames(panel), TRUE, FALSE))),
             bsButton("panel_add_row", "Add Row", disabled = ifelse(is.null(panel), TRUE, FALSE)),
             bsButton("panel_delete_rows", "Delete Selected Rows", disabled = ifelse(is.null(panel), TRUE, FALSE)),
             bsButton('panel_reset', 'Reset Panel', disabled = ifelse(is.null(panel), TRUE, FALSE))
             )
      )
    tableMetadata <- fluidRow(
      column(12, h5(HTML("<p style=\"color:grey;\">Metadata</p>"))),
      column(
      12, h5(HTML("<b style=\"color:grey;\">If you provide metadata, you can alter its cells by double-clicking. If you want to create metadata from your fcs files, upload your files first and then click on 'Load Data'. Attention: If you write numbers in your column, the feature will be interpreted as continous. If you do not want that (e.g., in the case of patient ID), add some letters (e.g., P1). </b>")), hr(), DTOutput("metadataDT")),
      column(12, 
             bsButton("md_add_column", "Add Column", disabled = ifelse(is.null(md), TRUE, FALSE)),
             bsButton("md_remove_column", "Delete Selected Columns", disabled = ifelse(is.null(md), TRUE, FALSE))
      ))
  }else{
    tablePanel <- renderTable(
      checkNullTable(panel),
      caption = "Panel Data",
      caption.placement = "top"
    )
    tableMetadata <- renderTable(
      checkNullTable(md),
      caption = "Metadata",
      caption.placement = "top"
    )
  }
  value <-
    list(
      renderTable(
        checkNullTable(fcs),
        caption = "FCS Data",
        caption.placement = "top"
      ),
      column(tablePanel,
             div(
               downloadButton("downloadPanel", "Download Panel as CSV"),
               style = "float: right;"
             ),
             width = 12),
      column(tableMetadata,
             div(
               downloadButton("downloadMD", "Download Metadata as CSV"),
               style = "float: right;"
             ),
             width = 12)
    )
  }
  if (input$chooseDataTab == "dataUpload" &
      !is.null(input$fcsFiles)) {
    status <- "success"
  }else if (input$chooseDataTab == "dataExample" &
           input$exampleData != "") {
    status <- "success"
    info <- readRDS(file.path(input$exampleData, "help.rds"))
    value <- c(info,value)
  } else {
    if (input$chooseDataTab == "sceUpload" & !is.null(input$sceFile)){
      tmp <- readRDS(file.path(input$sceFile$datapath))
      if (class(tmp) == "SingleCellExperiment"){
        status <- "success"
      }
    }
  }
  
  if (status == "success") {
    value <- c(list(
      bsButton(
        "loadData",
        "Load Data",
        icon("database"),
        style = "success",
        block = TRUE
      )
    ),
    value)
  }
  
  shinydashboard::box(value, title = "Selected Data", status = status)
})


output$downloadPanel <- downloadHandler(
  filename = "panel.csv",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    if(input$chooseDataTab == "dataUpload"){
      panel <- reactiveVals$data$upload$panel
    }else if(input$chooseDataTab == "dataExample"){
      panel <- reactiveVals$data$example$panel
    }else{
      panel <- data.table()
    }
    write.csv(panel, file, row.names = FALSE, quote = FALSE)
    waiter_hide(id = "app")
  }
)

output$downloadMD <- downloadHandler(
  filename = "metadata.csv",
  content = function(file) {
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    if(input$chooseDataTab == "dataUpload"){
      md <- reactiveVals$data$upload$md
    }else if(input$chooseDataTab == "dataExample"){
      md <- reactiveVals$data$example$md
    }else{
      md <- data.table()
    }
    write.csv(md, file, row.names = FALSE, quote = FALSE)
    waiter_hide(id = "app")
  }
)

output$setColorPalette <- renderUI({
  div(
    radioButtons(
      "chosenColorPalette",
      label = "Available palettes",
      choices = c("Colorblind1", "Colorblind2: Wong", "Colorblind3: Tol bright", "Colorblind4: Tol vibrant", "Colorblind5: Tol muted", 
                  "RColorBrewer: Set1", "RColorBrewer: Set2", "RColorBrewer: Set3", "RColorBrewer: Pastel1", "RColorBrewer: Pastel2", "RColorBrewer: Paired", "RColorBrewer: Dark2", "RColorBrewer: Accent",     
                  "CATALYST colors", "Other"),
      selected = "Colorblind1",
      inline = TRUE
    ),
    uiOutput(
      "otherCustomPalette"
    ),
    uiOutput(
      "renderedPalette"
    )
  )
})

output$otherCustomPalette <- renderUI({
  req(input$chosenColorPalette == "Other")
    div(
      textInput("colorInput", "Enter color name or code:", ""),
      bsButton("addColor", "Add Color", style = "info"),
      bsButton("clearColors", "Clear All Colors", style = "warning"),
      bsButton("setColorPalette", "Set Colors as Color Palette", style = "success")
  )
})

output$showColorPalette <- renderPlot({
  req(reactiveVals$selected_palette)
  return(scales::show_col(reactiveVals$selected_palette))
})

output$showColorPaletteOther <- renderPlot({
  req(reactiveVals$last_colors)
  return(scales::show_col(reactiveVals$last_colors))
})

output$renderedPalette <- renderUI({
  if(input$chosenColorPalette == "Other"){
    div(
      column(width=6, 
             h4("New palette"),
             plotOutput("showColorPaletteOther",  width = "100%", height = "500px")),
      column(width=6,
             h4("Current palette"),
             plotOutput("showColorPalette",  width = "100%", height = "500px"))
    )
  }else{
    plotOutput("showColorPalette",  width = "100%", height = "500px")
  }
})





