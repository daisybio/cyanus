# module goes here
library(shiny)

# Download UI element
downloadUI <- function(id) {  # nolint
  shiny::downloadButton(
    shiny::NS(id, "downloadData"),
    label = "Download",
    class = NULL,
    icon = shiny::icon("download"),
  )
}

# Downlaod server element
downloadServer <- function(id, harvestHere) {  # nolint
  # harvestHere should be name of the slot in which the plot is rendered to

  shiny::moduleServer(id, function(input, output, session) {

    output$downloadData <- shiny::downloadHandler(

      filename = function() {
        paste(
          as.character(harvestHere),
          "-plots-",
          Sys.Date(),
          ".zip",
          sep = ""
        )
      },
      content = function(file) {
        # zipping output[[harvestHere]] to file
        # We should run while loop to check if the file is ready
        # If not, we should wait for a while
        i <- 1
        while (i < 4) {
          Sys.sleep(3)
          i <- i + 1
        }
        zip(file, output[[harvestHere]])
      },
      contentType = "application/zip"
    )
  })
}
