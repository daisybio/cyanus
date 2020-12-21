startBody <- function() {
  box_height <- "12em"
  exampleDataVector <- c("MouseData" = "data/mousedata.rds",
                         "Platelets" = "data/platelets.rds")
  startBody <- tabItem(
    tabName = "start",
    fluidRow(shinydashboard::box(
      div("You can upload your own files in the FCS format or use an example dataset that we provide."),
      title = h2("Get Started"),
      width = 12
    )),
    fluidRow(
      shinydashboard::box(
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
        width = 4
      ),
      shinydashboard::box(
        fileInput(
          "fcsFiles",
          "Choose FCS File(s)",
          multiple = TRUE,
          accept = c(".fcs")
        ),
        title = "Upload Data",
        height = box_height,
        width = 8
      )
    ),
    shinycssloaders::withSpinner(uiOutput("currentData"))
  )
  return(startBody)
}
