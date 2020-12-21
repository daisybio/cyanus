header <- dashboardHeader(title = "CyTOF Pipeline")

sidebar <- dashboardSidebar(uiOutput("sidebar"))

# read all ui files
sapply(list.files("ui", full.names = TRUE), source)


body <-
  dashboardBody(
    useShinyjs(),
    tabItems(welcomeBody,
             startBody(), preprocessingBody()),
    actionButton("continue", "Start Analysis", icon("arrow-right"), class =
                   "btn-success btn-block")
  )


ui <- dashboardPage(header, sidebar, body)