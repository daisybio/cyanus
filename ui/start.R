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
    title = span("Upload Metadata", icon("question-circle"), id = "metaQ"),
    height = box_height,
    width = 3
  )
  
  metaPopover <- 
    bsPopover(
      id = "metaQ",
      title = "A CSV file with headers describing the experiment",
      content = "e.g. 4 columns:<br>file_name,sample_id,patient_id,condition<br>file_name: the FCS file name<br>sample_id: a unique sample identifier<br>patient_id: the patient ID<br>condition: brief sample description (e.g. reference/stimulated, healthy/diseased)"
    )
  
  panelUploadBox <- shinydashboard::box(
    fileInput("panelFile",
              "Choose Panel File (optional)",
              accept = c(".csv")),
    title = span("Upload Panel Data", icon("question-circle"), id = "panelQ"),
    height = box_height,
    width = 3
  )
  
  panelPopover <- 
    bsPopover(
      id = "panelQ",
      title = "A CSV file with headers describing the panel",
      content = "for each channel, its column name in the input data, targeted protein marker, and (optionally) class (type, state, or none) i.e.:<br>channel,col_name,marker,class"
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
    fluidRow(
      tabBox(
        tabPanel(
          fluidRow(dataUploadBox, metaUploadBox, metaPopover, panelUploadBox, panelPopover),
          value = "dataUpload",
          title = "Upload Data"
        ),
        tabPanel(
          fluidRow(dataExampleBox),
          value = "dataExample",
          title = "Example Data"
        ),
        id = "selectedData",
        title = "Choose Data",
        width = 12
      )
    ),
    shinycssloaders::withSpinner(uiOutput("currentData"))
  )
  return(startBody)
}