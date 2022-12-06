reactiveVals$data <- list(upload = list(fcs=NULL, panel=NULL, md=NULL), example = list(fcs=NULL, panel=NULL, md=NULL) )

observeEvent(input$fcsFiles, {
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
  if(any(!names(tmp_panel) %in% c("fcs_colname", "antigen", "marker_class")))
    showNotification(HTML("Error while reading the panel file:<br>A CSV or Excel file with headers describing the panel:<br>for each channel:<br>fcs_colname: its column name in the input data<br>antigen: targeted protein marker<br>marker_class: (optionally) class (type, state, or none)<br>i.e.:<br>fcs_colname,antigen[,marker_class]<br><b>Example: Check out the PBMC Example Data</b>"), duration = NULL, type = "error")
  else
    reactiveVals$data$upload$panel <- tmp_panel
})

observeEvent(input$sceFile, {
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
  
  if (input$chooseDataTab == "dataUpload") {
    dn <- dirname(input$fcsFiles$datapath)[1]
    file.rename(input$fcsFiles$datapath, file.path(dn, "/", input$fcsFiles$name))
    
    conditions <- names(reactiveVals$data$upload$md)[!names(reactiveVals$data$upload$md) 
                                                     %in% c("sample_id", "file_name")]
    md_cols <- list(file = "file_name", id = "sample_id", factors = conditions)
    
    tryCatch({
      withCallingHandlers({
        reactiveVals$sce <- CATALYST::prepData(
          dn,
          panel = reactiveVals$data$upload$panel,
          md = reactiveVals$data$upload$md,
          transform = FALSE,
          md_cols = md_cols,
          FACS = input$isFACSData,
          emptyValue = input$isEmptyValue
        )},
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
        showNotification(HTML(sprintf("Loading the data failed with the following message:<br>
                                    <b>%s</b>", e$message)), duration = NULL, type = "error")
      })
  } else if (input$chooseDataTab == "dataExample") {
    reactiveVals$sce <- readRDS(file.path(input$exampleData, "sce.rds"))
  } else if (input$chooseDataTab == "sceUpload") {
    reactiveVals$sce <- readRDS(file.path(input$sceFile$datapath))
    
  }else
    stop("Which tab is selected?")
  
  
  if (!is.null(reactiveVals$sce)) { # meaning the data loading worked
    
    #drop levels of markers
    SummarizedExperiment::rowData(reactiveVals$sce)$marker_class <- droplevels(SummarizedExperiment::rowData(reactiveVals$sce)$marker_class)
    
    # set negative values to zero
    if (length(which(assays(reactiveVals$sce)$counts < 0)) > 0){
      assays(reactiveVals$sce)$counts[assays(reactiveVals$sce)$counts < 0 ] <- 0
      if ("exprs" %in% assayNames(reactiveVals$sce)){
        assays(reactiveVals$sce)$exprs[assays(reactiveVals$sce)$exprs < 0 ] <- 0
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
    
    start_tab <- which(tab_ids == "start")
    if (isolate(reactiveVals$max_tab) > start_tab)
      reactiveVals$max_tab <- start_tab
    reactiveVals$continue[start_tab] <- TRUE
    runjs("document.getElementById('nextTab').scrollIntoView();")
    
  }
  updateButton(session, "loadData", label = " Load Data", disabled = FALSE) 
  waiter_hide(id = "app")
})

output$panelDT <- renderDT(
  checkNullTable(reactiveVals$data$upload$panel),
  editable = "cell"
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
      tableObj <- fluidRow(column(
        12, h5(HTML("<b style=\"color:grey;\">If you provide a panel, you can alter its cells by double-clicking</b>")), hr(), DTOutput("panelDT")))
      
      # if current data tab is example data  
    }else{
      tableObj <- renderTable(
        checkNullTable(panel),
        caption = "Panel Data",
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
        tableObj,
        renderTable(
          checkNullTable(md),
          caption = "Metadata",
          caption.placement = "top"
        )
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

observeEvent(input$panelDT_cell_edit, {
  reactiveVals$data$upload$panel <<- editData(reactiveVals$data$upload$panel, input$panelDT_cell_edit, "panelDT")
})