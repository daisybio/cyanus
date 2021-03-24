library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(plotly)
library(DT)

# read all ui files

sapply(list.files("ui", full.names = TRUE), source, environment())

header <-
  dashboardHeader(title = "CyTOF Pipeline", uiOutput("dashboard"))

sidebar <- dashboardSidebar(useShinyjs(), uiOutput("sidebar"))

body <-
  dashboardBody(
    fluidRow(column(6, bsButton(inputId = "previousTab", label = "Previous", icon = icon("spinner"), style = "warning", block = TRUE, disabled = TRUE)),
             column(6, bsButton(inputId = "nextTab", label = "Next", icon = icon("spinner"), style = "warning", block = TRUE, disabled = TRUE))),
    div(id = "loading", "This app is currently loading. Please be patient. This may take a minute."),
    tabItems(
      welcomeBody,
      startBody(),
      preprocessingBody(),
      visbody(),
      clusteringBody,
      deBody(),
      vennBody()
    )
  )


ui <- dashboardPage(header, sidebar, body)