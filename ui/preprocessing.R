# Preprocessing Tab
preprocessingBody <- function() {
  box_height <- "20em"
  transformationVector <- c("arcsinh","log")
  
  preprocessingBody <- tabItem(
    tabName = "preprocessing",
    fluidRow(shinydashboard::box(
      div("Preprocessing is essential in any mass cytometry analysis process. You have to choose a transformation to make the distributions more symmetric and to map them to a comparable range of expression."),
      title = h2("Data preprocessing"),
      width = 12
    )),
    fluidRow(
      shinydashboard::box(
        prettyRadioButtons(inputId="transformation", label="Possible transformations:", 
                     choices=c("arcsinh","log"),
                     selected = "arcsinh",
                     icon = icon("check"),
                     outline = TRUE),
      
        textInput("5","Cofactor of arcsinh transformation"),
        title = "Choose Transformation",
        height = box_height,
        width = 12
      ))
  )
  return(preprocessingBody)
}
