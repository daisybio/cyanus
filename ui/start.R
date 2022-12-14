startBody <- function() {
  box_height <- "12em"
  exampleDataVector <- c("PBMC" = "data/pbmc",
                         "Semi-Simulated COVID-19 data" = "data/covid_spiked",
                         "Simulated CytoGLMM data" = "data/cytoGLMM_simulated"
  )
  availableExampleData <- exampleDataVector %in% list.files("data", full.names = T)
  exampleDataVector <- exampleDataVector[availableExampleData]
  
  dataUploadBox <- shinydashboard::box(
    column(fileInput(
      "fcsFiles",
      "Choose FCS File(s)",
      multiple = TRUE,
      accept = c(".fcs")
    ), 
    width=9),
    column(
      checkboxInput("isFACSData", HTML("<b>FACS Data</b>"), FALSE),
      width= 3
    ),
    title = "Upload FCS Data",
    height = box_height,
    width = 6
  )
  
  metaUploadBox <- shinydashboard::box(
    fileInput(
      "metaFile",
      "Choose Metadata File (optional)",
      accept = c(".csv", ".xlsx", ".xls")
    ),
    title = span("Upload Metadata", icon("circle-question"), id = "metaQ"),
    height = box_height,
    width = 3
  )
  
  metaPopover <- 
    bsPopover(
      id = "metaQ",
      title = "A CSV or Excel file with headers describing the experiment",
      content = "e.g. 4 columns:<br>file_name, sample_id, patient_id, condition<br>file_name: the FCS file name<br>sample_id: a unique sample identifier<br>patient_id: the patient ID<br>condition: brief sample description (e.g. reference/stimulated, healthy/diseased)<br><b>Example: Check out the PBMC Example Metadata</b>."
    )
  
  panelUploadBox <- shinydashboard::box(
    fileInput("panelFile",
              "Choose Panel File (optional)",
              accept = c(".csv", ".xlsx", ".xls")),
    title = span("Upload Panel Data", icon("circle-question"), id = "panelQ"),
    height = box_height,
    width = 3
  )
  
  panelPopover <- 
    bsPopover(
      id = "panelQ",
      title = "A CSV or Excel file with headers describing the panel",
      content = "for each channel:<br>fcs_colname: its column name in the input data<br>antigen: targeted protein marker<br>marker_class: (optionally) class (type, state, or none)<br>i.e.:<br>fcs_colname,antigen[,marker_class]<br><b>Example: Check out the PBMC Example Data</b>"
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
  
  sceUploadBox <- shinydashboard::box(
    fileInput(
      "sceFile",
      "Upload SCE object",
      accept = c(".rds")
    ),
    title = span("Upload SCE object", icon("circle-question"), id = "sceObj"),
    height = box_height,
    width = 12
  )
  
  scePopover <- 
    bsPopover(
      id = "sceObj",
      title = "An .rds file with the SCE object.",
      content = "If you have saved the SCE object from previous analysis, you can upload the object again and continue your analysis."
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
          fluidRow(dataUploadBox, panelUploadBox, panelPopover, metaUploadBox, metaPopover),
          value = "dataUpload",
          title = "Upload Data"
        ),
        tabPanel(
          fluidRow(sceUploadBox, scePopover),
          value = "sceUpload",
          title = "Upload SCE"
        ),
        tabPanel(
          fluidRow(dataExampleBox),
          value = "dataExample",
          title = "Example Data"
        ),
        id = "chooseDataTab",
        title = "Choose Data",
        width = 12
      )
    ),
    shinycssloaders::withSpinner(uiOutput("currentData"))
  )
  return(startBody)
}
