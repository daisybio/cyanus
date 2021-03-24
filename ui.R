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
    )),
     div(
       id = "loading",
       fluidRow(box(
         div ("CyTOF is an uprising method for discovering biomarkers of the immune system. 
                              Its advantage over Flow Cytometry is that it labels the antibodies with metal 
                              isotopes instead of fluorophores, allowing it to analyse more markers in a single 
                              run while needing fewer cells."),
         div ("Consequently, <b>CyTOF</b> experiments are becoming a powerful tool to find 
                              immune marker expression differences between different conditions. In order 
                              to facilitate the analysis of CyTOF data for biologists and physicians, 
                              a clear, understandable and user-friendly pipeline is needed."),
         title = h1("Welcome to the CyTOF Pipeline"),
         width = 12
       )),
     ),
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