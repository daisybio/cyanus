tabs <- list(
  menuItem(
    "Welcome",
    tabName = tab_ids[1],
    icon = icon("home")
  ),
  menuItem(
    "Get Started",
    tabName = tab_ids[2],
    icon = icon("play-circle")
  ),
  menuItem(
    "Preprocessing",
    tabName = tab_ids[3],
    icon = icon("tools")
  ),
  menuItem(
    "Visualization",
    tabName = tab_ids[4],
    icon = icon("palette")
  ),
  menuItem(
    "Clustering",
    tabName = tab_ids[5],
    icon = icon("border-none")
  ),
  menuItem(
    "DE analysis",
    tabName = tab_ids[6],
    icon = icon("chart-bar")
  )
)

observeEvent(input$continue, {
  reactiveVals$current_tab <- reactiveVals$current_tab + 1
  shinyjs::hide("continue")
})

# grow menu depending on current tab
output$sidebar <- renderUI({
  curr_menu <- sidebarMenu(id = "tabs",
                           tabs[1:reactiveVals$current_tab])
  updateTabItems(session, "tabs", tab_ids[reactiveVals$current_tab])
  shinyjs::runjs("window.scrollTo(0, 0)")
  return(curr_menu)
})