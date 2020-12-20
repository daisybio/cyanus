library(shiny)
library(shinydashboard)

header <- dashboardHeader(title = "CyTOF Pipeline")

sidebar <- dashboardSidebar(sidebarMenu(
  menuItem("Welcome", tabName = "welcome", icon = icon("home")),
  menuItem(
    "Get Started",
    tabName = "start",
    icon = icon("play-circle")
  ),
  menuItem(
    "Visualization",
    tabName = "visualization",
    icon = icon("print")
  ),
  menuItem(
    "Clustering",
    tabName = "clustering",
    icon = icon("border-none")
  ),
  menuItem(
    "DE analysis",
    tabName = "diff_expr",
    icon = icon("chart-bar")
  )
))

sapply(list.files("ui/", full.names = TRUE), source)


body <-
  dashboardBody(tabItems(# First tab content
    welcomeBody, startBody))


ui <- dashboardPage(header, sidebar, body)