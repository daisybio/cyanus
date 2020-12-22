visbody <- function(){
  box_height <- "12em"
  example_data <- c("Platelets (4)" = "data/plateletsSCE1-4.Rds", 
                    "Platelets (6)" = "data/plateletsSCE1-6.Rds")
  vis_methods <- c("UMAP", "T-SNE", "PCA", "Isomap")
  
  dataExampleBox <- shinydashboard::box(
    selectizeInput(
      "exampleDataVis",
      "Example Datasets Visualization",
      example_data,
      options = list(
        placeholder = "Select Example Dataset...",
        onInitialize = I("function() { this.setValue(''); }")
      )
    ),
    title = "Choose Example Data",
    height = box_height,
    width = 12
  )
  
  visualizationBox <- shinydashboard::box(
    selectizeInput(
      "selectedVisMethod",
      "Visualization method",
      vis_methods,
      options = list(
        placeholder = "Select how you want to visualize your data",
        onInitialize = I("function() { this.setValue(''); }")
      )
    ),
    title = "Choose your Visualization method",
    height = box_height,
    width = 12
  )
  
    
  visbody <- tabItem(
    tabName = "visualization",
    fluidRow(shinydashboard::box(
      div(
        "Data Visualization can be done here"
      ),
      title = h2("Data Visualization"),
      width = 12
    )),
    fluidRow(
      tabPanel(fluidRow(dataExampleBox), value = "dataExample", title = "Example Data"),
      id = "selectedDataVis",
      title = "Choose Data",
      width = 12
    ),
    fluidRow(
      tabPanel(fluidRow(visualizationBox), value = "visMethod", title = "Visualization Method"),
      id = "visMethod",
      title = "Choose Visualization",
      width = 12
    ),
    fluidRow(
      box(plotOutput("plot1", height = 250))
    ),
    shinycssloaders::withSpinner(uiOutput("currentDataVis"))
    
  )
  return(visbody)
}




