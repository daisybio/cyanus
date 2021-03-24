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

sidebar <- dashboardSidebar(useShinyjs(),
                            tags$head(
                              tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }")
                            ),
                            uiOutput("sidebar"),
                            collapsed = TRUE)

body <-
  dashboardBody(
    div(
      id = "loading",
      "This app is currently loading. Please be patient. This may take a minute."
    ),
    tabItems(
      welcomeBody,
      startBody(),
      preprocessingBody(),
      visbody(),
      clusteringBody,
      deBody(),
      vennBody()
    ),
    fluidRow(column(
      6,
      bsButton(
        inputId = "previousTab",
        label = "App",
        icon = icon("spinner"),
        style = "warning",
        block = TRUE,
        disabled = TRUE
      )
    ),
    column(
      6,
      bsButton(
        inputId = "nextTab",
        label = "Loading",
        icon = icon("spinner"),
        style = "warning",
        block = TRUE,
        disabled = TRUE
      )
    ))
  )


ui <- dashboardPage(header, sidebar, body)