checkNullTable <- function(toCheck) {
  if (is.null(toCheck))
    return(data.frame("Nothing" = ""))
  else
    return(toCheck)
}

observeEvent(input$fcsFiles, {
  fileTable <- input$fcsFiles
  fileTable <- fileTable[, c("name", "size")]
  fileTable$size <- sprintf("%.2f MB", fileTable$size / 1000000)
  reactiveVals$fcs <- fileTable
})

observeEvent(input$metaFile, {
  if(endsWith(input$metaFile$datapath, ".csv")){
    reactiveVals$md <- read.table(input$metaFile$datapath, header = T, sep = ",")
  }else if(endsWith(input$metaFile$datapath, ".xls") | 
           endsWith(input$metaFile$datapath, ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    reactiveVals$md <- read.xlsx2(input$metaFile$datapath, 1)
  }
})

observeEvent(input$panelFile, {
  if(endsWith(input$panelFile$datapath, ".csv")){
    reactiveVals$panel <-
      read.table(input$panelFile$datapath, header = T, sep = ",")
  }else if(endsWith(input$panelFile$datapath, ".xls") | 
            endsWith(input$panelFile$datapath, ".xlsx")){
    library(xlsx)
    showNotification("There are often problems with reading Excel files in. If you can, please upload a .csv file", type = "warning")
    reactiveVals$panel <- read.xlsx2(input$panelFile$datapath, 1)
  }
})

observeEvent(input$exampleData, {
  reactiveVals$fcs <- readRDS(file.path(input$exampleData, "fcs.rds"))
  reactiveVals$panel <- readRDS(file.path(input$exampleData, "panel.rds"))
  reactiveVals$md <- readRDS(file.path(input$exampleData, "md.rds"))
}, ignoreInit = TRUE)

observeEvent(input$loadData, {
  updateButton(session, "loadData", label = " Loading...", disabled = TRUE)
  library(CATALYST)
  if (input$chooseDataTab == "dataUpload") {
    dn <- dirname(input$fcsFiles$datapath)[1]
    file.rename(input$fcsFiles$datapath, file.path(dn, "/", input$fcsFiles$name))
    reactiveVals$sce <- CATALYST::prepData(
      dn,
      reactiveVals$panel,
      reactiveVals$md,
      transform = FALSE
      #TODO: check if we have other columns
      #panel_cols = names(reactiveVals$panel),
      #md_cols = names(reactiveVals$md)
    )
  } else if (input$chooseDataTab == "dataExample") {
    reactiveVals$sce <- readRDS(file.path(input$exampleData, "sce.rds"))
  } else
    stop("Which tab is selected?")
  updateButton(session, "loadData", label = " Load Data", disabled = FALSE)
  updateButton(session, "continue", label = " Preprocessing")
  shinyjs::show("continue")
  runjs("document.getElementById('continue').scrollIntoView();")
})



output$panelDT <- renderDT(
  checkNullTable(reactiveVals$panel),
  editable = "cell"
)

output$currentData <- renderInfoBox({
  status <- "warning"
  if(input$chooseDataTab == "dataUpload"){
    tableObj <- fluidRow(column(
      12, h5(HTML("<b style=\"color:grey;\">If you provide a panel, you can alter its cells by double-clicking</b>")), hr(), DTOutput("panelDT")))
  }else{
    tableObj <- renderTable(
      checkNullTable(reactiveVals$panel),
      caption = "FCS Data",
      caption.placement = "top"
    )
  }
  value <-
    list(
      renderTable(
        checkNullTable(reactiveVals$fcs),
        caption = "FCS Data",
        caption.placement = "top"
      ),
      tableObj,
      renderTable(
        checkNullTable(reactiveVals$md),
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
    value <-
      c(list(div(
        sprintf(
          "Found %s: %s",
          input$exampleData,
          file.exists(input$exampleData)
        )
      ),
      div(
        sprintf("info about %s", input$exampleData)
      )),
      value)
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
  reactiveVals$panel <<- editData(reactiveVals$panel, input$panelDT_cell_edit, "panelDT")
})
