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
  dashboardHeader(title = span(img(src="logo_title.png", height="50px")), uiOutput("dashboard"))

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
    use_waiter(spinners = 5), # include dependencies
    tags$style( # manually make waiter overly for the whole page
      ".waiter-overlay-content{
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
         /* remove focus from success buttons */
         .btn-success:focus{
                              background-color: #00a65a;
                              border-color: #008d4c;
          }
                              
                              '))),
    waiter_show_on_load(html = tagList(spin_square_circle(), 
                                       HTML("<br>Loading App...")), 
                        color="rgb(146, 180, 242, .5)"),
    div(
      id = "loading",
      fluidRow(box(
        div ("Recently, high-dimensional time-of-flight mass cytometry (CyTOF) 
                               has emerged with the ability to identify more than 40 parameters simultaneously.
                               Traditional flow cytometry would require multiple tubes with different 
                               antibody panels to cover the same number of markers. Consequently, 
                               CyTOF experiments are becoming a powerful tool to unveil new cell subtypes, 
                               functions, and biomarkers in many fields, e.g. the discovery of disease-associated 
                               immunologic changes in cancer."),
        div ("In order to facilitate the analysis of CyTOF data for biologists and physicians, 
                               a clear, understandable and user-friendly pipeline is needed."),
        div("Here, we integrated the methods from the CATALYST package for preprocessing, 
                             visualization, and clustering. For differential abundance detection, we included the 
                             diffcyt methods diffcyt-DA-edgeR, diffcyt-DA-voom and diffcyt-DA-GLMM."),
        div("However, many experiments aim to detect differential states
                             within cell populations between samples in different conditions. 
                             For this, we integrated the published methods diffcyt-DS-limma, 
                             diffcyt-DS-LMM, CytoGLMM, and CytoGLM. Additionally, we performed a comprehensive analysis
                             of these existing methods and novel approaches published in …. 
                             Since the Wilcoxon rank-sum test and the t-test on sample medians, as well as
                             our novel method CyEMD performed well, we made them available in this interface as well."),
        div(img(src="cyanus_logo.png", height="150px", style="float:right; padding:20px;")),
        title = h1("Welcome to CYANUS: CYtometry ANalysis Using Shiny"),
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

ui <- tags$div(id = "app", dashboardPage(header, sidebar, body, title = "CYANUS"))