library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(plotly)
library(DT)
library(waiter)

jscode <- "shinyjs.closewindow = function() { window.close(); }"

# read all ui files
sapply(list.files("ui", full.names = TRUE), source, environment())

header <-
  dashboardHeader(title = "CyTOF Pipeline", uiOutput("dashboard"))

sidebar <- dashboardSidebar(useShinyjs(),
                            extendShinyjs(text = jscode, functions = c("closewindow")),
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
    use_waiter(spinners = 5), # include dependencies
    tags$style( # manually make waiter overly for the whole page
      ".waiter-overlay-content{
      height: 100vh;
      position: absolute;
      top: 50px; /*30 pixels from the top*/
      right: 50px; /*48% from the right*/
    }"
    ),
    tags$head(tags$style(HTML('
        /* logo */
        .skin-blue .main-header .logo {
                              background-color: #6495ed;
                              }

        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #6495ed;
                              }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #92b4f2;
                              }        

        /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #92b4f2;
                              }

        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #3676e8;
                              }

        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #6495ed;
                              color: #fffff;
                              }

        /* other links in the sidebarmenu when hovered */
         .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #6495ed;
                              }
        /* toggle button when hovered  */                    
         .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #6495ed;
                              }
                              '))),
    waiter_show_on_load(html = tagList(spin_square_circle(), 
                                       HTML("<br>Loading App...")), 
                        color="rgb(0, 102, 204, .5)"),
    div(
      id = "loading",
      fluidRow(box(
        div ("CyTOF is an uprising method for discovering biomarkers of the immune system. 
                              Its advantage over Flow Cytometry is that it labels the antibodies with metal 
                              isotopes instead of fluorophores, allowing it to analyse more markers in a single 
                              run while needing fewer cells."),
        div ("Consequently, CyTOF experiments are becoming a powerful tool to find 
                              immune marker expression differences between different conditions. In order 
                              to facilitate the analysis of CyTOF data for biologists and physicians, 
                              a clear, understandable and user-friendly pipeline is needed."),
        title = h1("Welcome to the CyTOF Pipeline"),
        width = 12
      ))
    ),
    tabItems(
      welcomeBody,
      startBody(),
      preprocessingBody(),
      visbody(),
      clusteringBody,
      deBody(),
      vennBody(),
      goodbyeBody
    ),
    fluidRow(column(
      6,
      shinyjs::hidden(bsButton(
        inputId = "previousTab",
        label = "",
        icon = icon("spinner"),
        style = "warning",
        block = TRUE,
        disabled = TRUE
      ))
    ),
    column(
      6,
      bsButton(
        inputId = "nextTab",
        label = "App Loading",
        icon = icon("spinner"),
        style = "warning",
        block = TRUE,
        disabled = TRUE
      )
    ))
)

ui <- tags$div(id = "app", dashboardPage(header, sidebar, body))