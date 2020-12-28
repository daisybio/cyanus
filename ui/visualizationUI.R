markerSelectionBox <- shinydashboard::tabBox(
  tabPanel(uiOutput("markerClassVis"), 
           value = "classTab",
           title = "By Marker Class"),
  tabPanel(uiOutput("expressionVis"),
           value = "expressionTab",
           title = "By Marker"),
  id = "visTabs",
  title = "Run your dimensionality reduction", 
  width = 12
)

visbody <- function(){
  
  visualizationBox <- shinydashboard::box(
    uiOutput("methodsVis"),
    uiOutput("assayVis"),
    radioButtons(
      inputId = "scaleVis",
      label = "Scale ",
      choices = c("yes", "no"), 
      inline = TRUE
    ),
    uiOutput("color_by"),
    uiOutput("facet_by"),
    div(
      bsButton(
        "startDimRed",
        "Start Visualization",
        icon = icon("palette"),
        style = "success", 
        disabled = TRUE
      ),
      style = "float: right;", 
      id = "divStartDimRed"
    ),
    id = "visBox", 
    title = "Visualize your dimensionality reduction(s)", 
    width = 6
  )
  
  plotBox <- shinydashboard::box(
    fluidRow(
      box(shinycssloaders::withSpinner(plotlyOutput("visPlot", width = "90%")), 
          id = "imageBox",
          title = "Plot",
          width = 8),
      uiOutput("plotInfo")
    ),
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
      box(
        markerSelectionBox,
        uiOutput("runDRparBox"),
        width = 6
      ),
        visualizationBox,
    ),
    fluidRow(
      plotBox
    )
  )
  return(visbody)
}




