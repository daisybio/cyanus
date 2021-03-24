checkNullTable <- function(toCheck) {
  if (is.null(toCheck))
    return(data.frame("Nothing" = ""))
  else
    return(toCheck)
}

reactiveVals$data <- list(upload = list(fcs=NULL, panel=NULL, md=NULL), example = list(fcs=NULL, panel=NULL, md=NULL) )

observeEvent(input$fcsFiles, {
  fileTable <- input$fcsFiles
  fileTable <- fileTable[, c("name", "size")]
  fileTable$size <- sprintf("%.2f MB", fileTable$size / 1000000)
  reactiveVals$data$upload$fcs <- fileTable
})

observeEvent(input$metaFile, {
  if(endsWith(input$metaFile$datapath, ".csv")){
    reactiveVals$data$upload$md <- read.table(input$metaFile$datapath, header = T, sep = ",")
  }else if(endsWith(input$metaFile$datapath, ".xls") | 
           endsWith(input$metaFile$datapath, ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    reactiveVals$data$upload$md <- read.xlsx2(input$metaFile$datapath, 1)
  }
})

observeEvent(input$panelFile, {
  if(endsWith(input$panelFile$datapath, ".csv")){
    reactiveVals$data$upload$panel <-
      read.table(input$panelFile$datapath, header = T, sep = ",")
  }else if(endsWith(input$panelFile$datapath, ".xls") | 
            endsWith(input$panelFile$datapath, ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    reactiveVals$data$upload$panel <- read.xlsx2(input$panelFile$datapath, 1)
  }
})

observeEvent(input$exampleData, {
  reactiveVals$data$example$fcs <- readRDS(file.path(input$exampleData, "fcs.rds"))
  reactiveVals$data$example$panel <- readRDS(file.path(input$exampleData, "panel.rds"))
  reactiveVals$data$example$md <- readRDS(file.path(input$exampleData, "md.rds"))
}, ignoreInit = TRUE)


observeEvent(input$loadData, {
  updateButton(session, "loadData", label = " Loading...", disabled = TRUE)
  library(CATALYST)
  if (input$chooseDataTab == "dataUpload") {
    dn <- dirname(input$fcsFiles$datapath)[1]
    file.rename(input$fcsFiles$datapath, file.path(dn, "/", input$fcsFiles$name))

    conditions <- names(reactiveVals$data$upload$md)[!names(reactiveVals$data$upload$md) 
                                                     %in% c("sample_id", "file_name")]
    md_cols <- list(file = "file_name", id = "sample_id", factors = conditions)
    
    reactiveVals$sce <- CATALYST::prepData(
      dn,
      panel = reactiveVals$data$upload$panel,
      md = reactiveVals$data$upload$md,
      transform = FALSE,
      md_cols = md_cols
      #TODO: check if we have other columns
      #panel_cols = names(reactiveVals$panel),
    )
  } else if (input$chooseDataTab == "dataExample") {
    reactiveVals$sce <- readRDS(file.path(input$exampleData, "sce.rds"))
  } else
    stop("Which tab is selected?")
  updateButton(session, "loadData", label = " Load Data", disabled = FALSE)
  reactiveVals$continue <- TRUE
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
  } else {
    fcs <- reactiveVals$data$example$fcs
    panel <- reactiveVals$data$example$panel
     md <- reactiveVals$data$example$md
  }
  
  status <- "warning"
  # if current data tab is upload data
  if(input$chooseDataTab == "dataUpload"){
    tableObj <- fluidRow(column(
      12, h5(HTML("<b style=\"color:grey;\">If you provide a panel, you can alter its cells by double-clicking</b>")), hr(), DTOutput("panelDT")))
  
  # if current data tab is example data  
  }else{
    tableObj <- renderTable(
      checkNullTable(panel),
      caption = "FCS Data",
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
  
  if (input$chooseDataTab == "dataUpload" &
      !is.null(input$fcsFiles)) {
    status <- "success"
  }
  else if (input$chooseDataTab == "dataExample" &
           input$exampleData != "") {
    status <- "success"
    info <- readRDS(file.path(input$exampleData, "help.rds"))
    value <- c(info,value)
    
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
