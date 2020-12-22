# Preprocessing Tab
preprocessingBody <- function() {
  box_height <- "20em"
  
  transformationBox <- shinydashboard::box(
    prettyRadioButtons(inputId="transformation", label="Possible transformations:", 
                       choices=c("no","log","arcsinh"),
                       selected = "no",
                       icon = icon("check"),
                       outline = TRUE),
    title = "Choose Transformation",
    height = "12em",
    width = 6
  )
  
  cofactorBox <- shinydashboard::box(
    conditionalPanel(condition = "input.transformation=='arcsinh'",
                     textInput("cofactor","Cofactor:", value="5")),
    title = "Choose Cofactor of Arcsinh transformation",
    height = "12em",
    width = 6
  )
  
  markersBox <- shinydashboard::box(
    selectizeInput("markerSelection", choices = NULL, label="Markers", multiple=TRUE),
    title = "Select Markers",
    height = "12em",
    width = 6
  )
  
  samplesBox <- shinydashboard::box(
    selectizeInput("sampleSelection", choices = NULL, label="Samples", multiple=TRUE),
    title = "Select Samples",
    height = "12em",
    width = 6
  )
  
  countsBox <- shinydashboard::box(
    shinycssloaders::withSpinner(plotOutput("countsPlot", width="100%", height="350px")),
    title = "Barplot showing the numbers of cells measured for each sample",
    width = 12,
    height = "30em"
  )
  
  mdsBox <- shinydashboard::box(
    shinycssloaders::withSpinner(plotOutput("mdsPlot", width="100%", height="350px")),
    title = "MDS plot",
    width = 12,
    height = "30em"
  )
  
  nrsBox <- shinydashboard::box(
    shinycssloaders::withSpinner(plotOutput("nrsPlot", width="100%", height="350px")),
    title = "NRS plot",
    width = 12,
    height = "30em"
  )
  
  exprsBox <- shinydashboard::box(
    shinycssloaders::withSpinner(plotOutput("exprsPlot", width="100%", height="350px")),
    title = "Plot with per-sample marker expression distributions",
    width = 12,
    height = "30em"
  )
  
  exprsHeatmapBox <- shinydashboard::box(
    shinycssloaders::withSpinner(plotOutput("exprsHeatmapPlot", width="100%", height="350px")),
    title = "Heatmap of marker expression across all cells measured for each sample",
    width = 12,
    height = "30em"
  )
  
  preprocessingBody <- tabItem(
    tabName = "preprocessing",
    fluidRow(shinydashboard::box(
      div("Preprocessing is essential in any mass cytometry analysis process. You have to choose a transformation to make the distributions more symmetric and to map them to a comparable range of expression."),
      title = h2("Data preprocessing"),
      width = 12
    )),
    #hide(id="plots", anim=TRUE)
   # hide(id="markers", anim=TRUE)
    #hide(id="samples", anim=TRUE)
    
  fluidRow(transformationBox, cofactorBox),
    
  actionButton("prepButton", "Start Preprocessing", icon("arrow-right"), class =
                   "btn-success btn-block"),
    
  fluidRow(div(id="plots", tabBox(
    tabPanel(fluidRow(countsBox), value = "plotCounts", title = "Counts"),
    tabPanel(fluidRow(mdsBox), value = "plotMDS", title = "MDS"),
    tabPanel(fluidRow(nrsBox), value = "plotNRS", title = "NRS"),
    tabPanel(fluidRow(exprsBox), value = "plotExpr", title = "Expr"),
    tabPanel(fluidRow(exprsHeatmapBox), value = "plotHeatmapExpr", title = "Heatmap"),
    id = "plots",
    title = "Simple Data Visualization",
    width = 12
  ))),
      
  fluidRow(div(id="markers", markersBox), div(id="samples",samplesBox))
    
  )
  
  return(preprocessingBody)
}
