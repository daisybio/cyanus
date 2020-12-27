visbody <- function(){
  
  vis_methods <- c("UMAP", "T-SNE", "PCA", "Isomap")
  
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
    uiOutput("assayVis"),
    radioButtons(
      inputId = "scaleVis",
      label = "Scale ",
      choices = c("yes", "no"), 
      inline = TRUE
    ),
    numericInput(
      "nrCells",
      label = "From how many cells do you want to sample?",
      value = 100,
      min = 0,
      max = 100000,
      step = 100
      ),
    uiOutput("color_by"),
    uiOutput("facet_by"),
    div(
      bsButton(
        "startDimRed",
        "Start Dimensionality Reduction",
        icon = icon("palette"),
        style = "success", 
        disabled = TRUE
      ),
      style = "float: right;", 
      id = "divStartDimRed"
    ),
    id = "visBox", 
    title = "Choose your Visualization method",
    width = 6
  )
  
  
  plotBox <- shinydashboard::box(shinycssloaders::withSpinner(plotOutput("visPlot", width="100%", height="350px")),
                                 id = "visPlotBox",
                                 title = "Dimensionality Reduction",
                                 width = 12)
  
    
  visbody <- tabItem(
    tabName = "visualization",
    fluidRow(
      shinydashboard::box(
      div("Data Visualization can be done here."
      ),
      div("You can either apply the dimensionality reduction for a subset of markers or by marker class."),
      title = h2("Data Visualization"),
      width = 12)
      ),
    fluidRow(
      uiOutput("parametersVis"),
      visualizationBox
    ),
    fluidRow(
      plotBox
    )
  )
  return(visbody)
}




