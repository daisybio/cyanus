server <- function(input, output) {
  
  ### start ----
  output$currentData <- renderInfoBox({
    status <- "info"
    if (input$exampleData == "" & is.null(input$fcsFiles))
      value <- "You have currently no data selected"
    else if (!is.null(input$fcsFiles)){
      value <- "info about uploaded data"
    }
    else if (input$exampleData != ""){
      value <- "info about example data"
    } else stop("This is unexpected. Is Data selected or not?")
    box(value, title = "Selected Data", status = status)
  })
  
  
}