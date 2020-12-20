library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)

reactiveVals <- reactiveValues()
tab_ids <- c("welcome", "start", "visualization", "clustering", "de")
tabs <- list(
  menuItem(
    "Welcome",
    tabName = tab_ids[1],
    icon = icon("home")
  ),
  menuItem(
    "Get Started",
    tabName = tab_ids[2],
    icon = icon("play-circle")
  ),
  menuItem(
    "Visualization",
    tabName = tab_ids[3],
    icon = icon("print")
  ),
  menuItem(
    "Clustering",
    tabName = tab_ids[4],
    icon = icon("border-none")
  ),
  menuItem(
    "DE analysis",
    tabName = tab_ids[5],
    icon = icon("chart-bar")
  )
)

reactiveVals$current_tab <- 1
reactiveVals$useExampleData <- F
reactiveVals$useUploadedData <- F


