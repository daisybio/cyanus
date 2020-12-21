startBody <- function() {
  box_height <- "12em"
  exampleDataVector <- c("MouseData" = "data/mousedata.rds",
                         "Platelets" = "data/platelets.rds")

  
  dataUploadBox <- shinydashboard::box(
    fileInput(
      "fcsFiles",
      "Choose FCS File(s)",
      multiple = TRUE,
      accept = c(".fcs")
    ),
    title = "Upload FCS Data",
    height = box_height,
    width = 6
  )
  
  metaUploadBox <- shinydashboard::box(
    fileInput(
      "metaFile",
      "Choose Metadata File (optional)",
      accept = c(".csv")
    ),
    title = "Upload Metadata",
    height = box_height,
    width = 3
  )
  
  panelUploadBox <- shinydashboard::box(
    fileInput(
      "panelFile",
      "Choose Panel File (optional)",
      accept = c(".csv")
    ),
    title = "Upload Panel Data",
    height = box_height,
    width = 3
  )
  
  dataExampleBox <- shinydashboard::box(
    selectizeInput(
      "exampleData",
      "Example Datasets",
      exampleDataVector,
      options = list(
        placeholder = "Select Example Dataset...",
        onInitialize = I("function() { this.setValue(''); }")
      )
    ),
    title = "Choose Example Data",
    height = box_height,
    width = 12
  )
  
  startBody <- tabItem(
    tabName = "start",
    fluidRow(shinydashboard::box(
      div(
        "You can upload your own files in the FCS format or use an example dataset that we provide."
      ),
      title = h2("Get Started"),
      width = 12
    )),
    fluidRow(tabBox(
      tabPanel(fluidRow(dataUploadBox, metaUploadBox, panelUploadBox), value = "dataUpload", title = "Upload Data"),
      tabPanel(fluidRow(dataExampleBox), value = "dataExample", title = "Example Data"),
      id = "selectedData",
      title = "Choose Data",
      width = 12
    )),
    shinycssloaders::withSpinner(uiOutput("currentData"))
  )
  return(startBody)
}
