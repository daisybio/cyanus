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
sapply(list.files("R/ui", full.names = TRUE), source, environment())

header <-
  dashboardHeader(title = span(img(src="logo_title.png", height="50px")), tags$li(uiOutput("downloadLog"), class='dropdown'), tags$li(uiOutput("dashboard"), class='dropdown'), tags$li(a(icon('envelope'), " Contact for problems/suggestions", href = "https://github.com/biomedbigdata/cyanus/issues", target = "_blank"), class='dropdown'))

sidebar <- dashboardSidebar(useShinyjs(),
                            extendShinyjs(text = jscode, functions = c("closewindow")),
                            tags$head(
                              tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }"),
                            ),
                            uiOutput("sidebar"),
                            collapsed = TRUE)

body <-
  dashboardBody(
    useWaiter(), # include dependencies
    titlePanel(windowTitle = "CYANUS",
               title = tags$head(tags$link(rel="icon",
                                           type="image/png",
                                           href="favicon.png"))),
    tags$style( # manually make waiter overly for the whole page
      ".waiter-overlay-content{
      color: black;
      font-weight: bold;
      height: 100vh;
      position: absolute;
      top: 50vh; /*30 pixels from the top*/
      right: 50%; /*48% from the right*/
    }"
    ),
    tags$head(
    tags$style(HTML('
        /* logo */
        .skin-blue .main-header .logo {
                              background-color: #ffffff;
                              }

        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #ffffff;
                              }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #8489c3;
                              }        

        /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #8489c3;
                              }

        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #b0b3e0;
                              }

        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #4a4d78;
                              color: #fffff;
                              }

        /* other links in the sidebarmenu when hovered */
         .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #4a4d78;
                              }
        /* toggle button when hovered  */                    
         .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #4a4d78;
         }
         /* remove focus from success buttons */
         .btn-success:focus{
                              background-color: #00a65a;
                              border-color: #008d4c;
          }
                              
                              '))),
    waiter_show_on_load(html = tagList(spin_loaders(id = 5, color = "black"), 
                                       HTML("<br>Loading App...")), 
                        color="rgb(146, 180, 242, .5)"),
    div(
      id = "loading",
      welcomeRow
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

ui <- tags$div(id = "app", dashboardPage(header, sidebar, body, title = "CYANUS"))