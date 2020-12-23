server <- function(input, output, session) {
  
  # make reactiveValues and server-wide variables
  tab_ids <- c("welcome", "start", "preprocessing", "visualization", "clustering", "de")
  reactiveVals <- reactiveValues()
  #reactiveVals$sce <- readRDS("data/plateletsSCE1-4.Rds")
  reactiveVals$current_tab <- 1
  
  # read all server files
  sapply(list.files("server", full.names = TRUE), source, environment())
  
  shinyjs::hide("visBox")
  shinyjs::hide("visPlotBox")
  ### menu ----
  
  
  ### start ----
  # observeEvent(input$exampleData, {
  #   if (input$exampleData != ""){
  #     reactiveVals$useExampleData <- T
  #     reactiveVals$useUploadedData <- F
  #   }
  # })
  # 
  # observeEvent(input$fcsFiles, {
  #   if (!is.null(input$fcsFiles)){
  #     reactiveVals$useExampleData <- F
  #     reactiveVals$useUploadedData <- T
  #   }
  # })
  
  #TODO: make tabs
  #TODO: make optional metadata input box
  #TODO: make panel information 
  
 
  
  
}
