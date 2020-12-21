library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)


# read all ui files
sapply(list.files("ui", full.names = TRUE), source, environment())

header <- dashboardHeader(title = "CyTOF Pipeline")

sidebar <- dashboardSidebar(uiOutput("sidebar"))

body <-
  dashboardBody(
    useShinyjs(),
    tabItems(welcomeBody,
             startBody()),
    actionButton("continue", "Start Analysis", icon("arrow-right"), class =
                   "btn-success btn-block")
  )


ui <- dashboardPage(header, sidebar, body)