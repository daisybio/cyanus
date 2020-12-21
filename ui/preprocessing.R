# Preprocessing Tab
preprocessingBody <- function() {
  box_height <- "20em"
  
  transformationBox <- shinydashboard::box(
    prettyRadioButtons(inputId="transformation", label="Possible transformations:", 
                       choices=c("no","log","arcsinh"),
                       selected = "no",
                       icon = icon("check"),
                       outline = TRUE),
    conditionalPanel(condition = "input.transformation=='arcsinh'",
                     textInput("5","Cofactor of arcsinh transformation")),
    title = "Choose Transformation",
    height = box_height,
    width = 12
  )

  plotCountsBox <- shinydashboard::box(
    plotOutput("plotCounts")
  )
  
  preprocessingBody <- tabItem(
    tabName = "preprocessing",
    fluidRow(shinydashboard::box(
      div("Preprocessing is essential in any mass cytometry analysis process. You have to choose a transformation to make the distributions more symmetric and to map them to a comparable range of expression."),
      title = h2("Data preprocessing"),
      width = 12
    )),
    fluidRow(transformationBox),
    
    actionButton("prepData", "Start Preprocessing", icon("arrow-right"), class =
                   "btn-success btn-block"),
    
  
    #shinycssloaders::withSpinner(uiOutput("currentTransformation"))
  )
  
  return(preprocessingBody)
}
