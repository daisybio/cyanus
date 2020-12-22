library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(shinyBS)


# read all ui files

sapply(list.files("ui", full.names = TRUE), source, environment())

header <- dashboardHeader(title = "CyTOF Pipeline")

sidebar <- dashboardSidebar(uiOutput("sidebar"))

body <-
  dashboardBody(
    useShinyjs(),
    tabItems(welcomeBody,
             startBody(),
             visbody(),
             clusteringBody),
    bsButton("continue", "Start Analysis", icon("arrow-right"), style = "success", block = TRUE)
  )


ui <- dashboardPage(header, sidebar, body)