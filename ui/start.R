box_height = "12em"
exampleDataVector <- c("MouseData" = "data/mousedata.rds",
                       "Platelets" = "data/platelets.rds")


startBody <- tabItem(
  tabName = "start",
  box(
    div("first you can upload your data or use example data"),
    title = h2("Get Started"),
    width = 12
  ),
  fluidRow(
    box(
      fileInput(
        "fcsFiles",
        "Choose FCS File(s)",
        multiple = TRUE,
        accept = c(".fcs")
      ),
      title = "Upload Data",
      height = box_height
    ),
    box(
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
      height = box_height
    )
  ),
  uiOutput("currentData")
)