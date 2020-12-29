# Preprocessing Tab

preprocessingBody <- function() {
  
  marker_sample_height <- "25em"
  panel_height <- "40em"
  plot_height <- "35em"
  
  # box with transformations: arcsinh, log or none
  transformationBox <- shinydashboard::box(
    textInput("cofactor", "Cofactor:", value ="5"),
    div(
      bsButton(
        "prepButton",
        "Start Transformation",
        icon = icon("tools"),
        style = "success"
      ),
      style = "float: right;"
    ),
    title = span("Choose Cofactor for Arcsinh Transformation", icon("question-circle"), id = "cofactor"),
    height = marker_sample_height,
    width = 6
  )
  
  cofactorPopover <- 
    bsPopover(
      id = "cofactor",
      title = "Cofactor of the inverse hyperbolic sine transformation",
      content = "Recommended values for the cofactor parameter are 5 for mass cytometry (CyTOF) or 150 for fluorescence flow cytometry."
    )
  
  # box with markers, samples and patients (all markers, patients, samples selected by default)
  selectingBox <- shinydashboard::box(
    uiOutput("markersBox"),
    uiOutput("samplesBox"),
    uiOutput("patientsBox"),
    div(
      bsButton(
        "prepVisButton",
        "Visualize Selection",
        icon = icon("palette"),
        style = "success"
      ),
      style = "float: right;"
    ),
    title = "Select Markers, Samples and Patients",
    height = marker_sample_height,
    width = 6
  )
  
  
  # box for counts plots
  countsBox <- shinydashboard::box(
    uiOutput("designCounts"),
    title = "Barplot showing the numbers of cells measured for each sample",
    width = 12,
    height = plot_height
  )
  
  # box for mds plots
  mdsBox <- shinydashboard::box(
    uiOutput("designMDS"),
    title = "MDS plot",
    width = 12,
    height = plot_height
  ) 
  
  # box for nrs plots
  nrsBox <- shinydashboard::box(
    uiOutput("designNRS"),
    title = "NRS plot",
    width = 12,
    height = plot_height
  ) 
  
  # box for exprs plots
  exprsBox <- shinydashboard::box(
    uiOutput("designExprs"),
    title = "Exprs",
    width = 12,
    height = plot_height
  ) 
  
  # box for exprs heatmpa plots
  exprsHeatmapBox <- shinydashboard::box(
    uiOutput("designExprsHeatmap"),
    title = "Heatmap",
    width = 12,
    height = plot_height
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
    
    # box for selecting transformation, markers, patients and samples
    fluidRow(transformationBox, cofactorPopover, selectingBox),
    

    # tabBox with simple visualization plots
    fluidRow(
      id = "plots",
      tabBox(
        tabPanel(fluidRow(countsBox), value = "plotCounts", title = "Counts"),
        tabPanel(fluidRow(mdsBox), value = "plotMDS", title = "MDS"),
        tabPanel(fluidRow(nrsBox), value = "plotNRS", title = "NRS"),
        tabPanel(fluidRow(exprsBox), value = "plotExpr", title = "Expr"),
        tabPanel(fluidRow(exprsHeatmapBox),value = "plotHeatmapExpr",title = "Heatmap"),
        id = "plots",
        title = "Simple Data Visualization",
        width = 12,
        height = panel_height
      )
    ),
  )
  return(preprocessingBody)
}
