# Preprocessing Tab

preprocessingBody <- function() {
  box_height <- "20em"
  plot_height <- "30em"
  
  # box with transformations: arcsinh, log or none
  transformationBox <- shinydashboard::box(
    prettyRadioButtons(
      inputId = "transformation",
      label = "Possible transformations:",
      choices = c("no", "log", "arcsinh"),
      selected = "no",
      icon = icon("check"),
      outline = TRUE
    ),
    div(
      bsButton(
        "prepButton",
        "Start Transformation",
        icon = icon("border-none"),
        style = "success"
      ),
      style = "float: right;"
    ),
    title = "Choose Transformation",
    height = "15em",
    width = 6
  )
  
  # box for specifying cofactor when selecting arcsinh transformation
  cofactorBox <- conditionalPanel(
    condition = "input.transformation=='arcsinh'",
    shinydashboard::box(
      textInput("cofactor", "Cofactor:", value =
                  "5"),
      title = "Choose Cofactor of Arcsinh transformation",
      height = "15em",
      width = 6
    )
  )
  
  # box with markers (all markers selected by default)
  markersBox <- shinydashboard::box(
    pickerInput(
      inputId = "markerSelection",
      label = "Markers",
      choices = NULL,
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    ),
      
    title = "Select Markers",
    height = "20em",
    width = 6
  )
  
  # box with samples (all samples selected by default)
  samplesBox <- shinydashboard::box(
    pickerInput(
      "sampleSelection",
      choices = NULL,
      label = "Samples",
      options = list(
        `actions-box` = TRUE,
        size = 4,
        `selected-text-format` = "count > 3",
        "dropup-auto" = FALSE
      ),
      multiple = TRUE
    ),
    title = "Select Samples",
    height = "20em",
    width = 6
  )
  
  # box for counts plots
  countsBox <- shinydashboard::box(
    uiOutput("designCounts"),
    shinydashboard::box(shinycssloaders::withSpinner(
      plotOutput("countsPlot", width = "100%", height = "350px")
    ), width = 8, heigth= plot_height),
    title = "Barplot showing the numbers of cells measured for each sample",
    width = 12,
    height = plot_height
  )
  
  mdsBox <- shinydashboard::box(
    uiOutput("designMds"),
    shinydashboard::box(
      plotOutput("mdsPlot", width = "100%", height = "350px"),
      width = 8,
    ),
    title = "MDS plot",
    width = 12
  )
  
  nrsBox <- shinydashboard::box(
    uiOutput("designNrs"),
    shinydashboard::box(
      plotOutput("nrsPlot", width = "100%", height = "350px"),
      width = 8,
    ),
    title = "NRS plot",
    width = 12
  )
  
  exprsBox <- shinydashboard::box(
    uiOutput("designExprs"),
    shinydashboard::box(
      plotOutput("exprsPlot", width = "100%", height = "350px"),
      width = 8,
    ),
    title = "Plot with per-sample marker expression distributions",
    width = 12
  )
  
  exprsHeatmapBox <- shinydashboard::box(
    uiOutput("designHeatmapExprs"),
    shinydashboard::box(
      plotOutput("exprsHeatmapPlot", width = "100%", height = "350px"),
      width = 8,
    ),
    title = "Heatmap of marker expression across all cells measured for each sample",
    width = 12
  )
  
  
  # Preprocessing body
  preprocessingBody <- tabItem(
    tabName = "preprocessing",
    fluidRow(shinydashboard::box(
      div(
        "Preprocessing is essential in any mass cytometry analysis process. You have to choose a transformation to make the distributions more symmetric and to map them to a comparable range of expression."
      ),
      title = h2("Data preprocessing"),
      width = 12
    )),
    
    # box for selecting transformations
    fluidRow(transformationBox, cofactorBox),
    
    # box selecting markers and box selecting samples
    fluidRow(markersBox, samplesBox),
    
    # tabBox with simple visualization plots
    fluidRow(id = "plots",
      tabBox(
        tabPanel(fluidRow(countsBox), value = "plotCounts", title = "Counts"),
        tabPanel(fluidRow(mdsBox), value = "plotMDS", title = "MDS"),
        tabPanel(fluidRow(nrsBox), value = "plotNRS", title = "NRS"),
        tabPanel(fluidRow(exprsBox), value = "plotExpr", title = "Expr"),
        tabPanel(
          fluidRow(exprsHeatmapBox),
          value = "plotHeatmapExpr",
          title = "Heatmap"
        ),
        id = "plots",
        title = "Simple Data Visualization",
        width = 12,
        height = "30em"
      )
    ),
  )
  return(preprocessingBody)
}
