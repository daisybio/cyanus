tab_ids <-
  c("welcome",
    "start",
    "preprocessing",
    "visualization",
    "clustering",
    "de",
    "venn",
    "goodbye")
reactiveVals$current_tab <- 1
reactiveVals$max_tab <- 1
reactiveVals$continue <- TRUE

tabs <- list(
  menuItem("Welcome",
           tabName = tab_ids[1],
           icon = icon("home")),
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
    "DE Analysis",
    tabName = tab_ids[6],
    icon = icon("chart-bar")
  ),
  menuItem(
    "DE Method Comparison",
    tabName = tab_ids[7],
    icon = icon("object-group")
  ),
  menuItem(
    "Goodbye",
    tabName = tab_ids[8],
    icon = icon("stop")
  )
)

observeEvent({
  reactiveVals$current_tab
  reactiveVals$continue
}, {
  if (reactiveVals$current_tab == 1)
    shinyjs::hide("previousTab")
  else if (reactiveVals$current_tab > 1)
    shinyjs::show("previousTab")
  else
    stop("what is the current id?")
  if (reactiveVals$current_tab < length(tab_ids))
    shinyjs::show("nextTab")
  else if (reactiveVals$current_tab == length(tab_ids))
    #shinyjs::hide("nextTab")
    reactiveVals$continue <- TRUE #so that Exit App button is directly enabled
  else
    stop("what is the current id?")
  if (!reactiveVals$continue &&
      reactiveVals$current_tab == reactiveVals$max_tab)
    shinyjs::disable("nextTab")
  else if (reactiveVals$continue || reactiveVals$current_tab != reactiveVals$max_tab) {
    shinyjs::enable("nextTab")
  }
})

prevNames <- c(" Previous", " Welcome", " Get Started", " Preprocessing", " Visualization", " Clustering", " DE Analysis", " DE Method Comparison")
nextNames <- c(" Get Started", " Preprocessing", " Visualization", " Clustering", " DE Analysis", " DE Method Comparison", " Goodbye", " Exit App")

observeEvent(input$tabs, {
  reactiveVals$current_tab <- match(input$tabs, tab_ids)
  shinyBS::updateButton(session, inputId = "previousTab", label = prevNames[reactiveVals$current_tab])
  shinyBS::updateButton(session, inputId = "nextTab", label = nextNames[reactiveVals$current_tab])
})

observeEvent(input$previousTab, {
  shinyjs::enable("nextTab")
  reactiveVals$current_tab <- reactiveVals$current_tab - 1
  updateTabItems(session, "tabs", tab_ids[reactiveVals$current_tab]) # here we update in case the new current tab is not a new max tab
  shinyjs::runjs("window.scrollTo(0, 0)")
})

observeEvent(input$nextTab, {
  shinyjs::enable("previousTab")
  reactiveVals$current_tab <- reactiveVals$current_tab + 1
  if (reactiveVals$current_tab > reactiveVals$max_tab){
    reactiveVals$max_tab <- reactiveVals$current_tab
    if (reactiveVals$current_tab == length(tab_ids) + 1){
      stopApp(7)
    }
  } else
    updateTabItems(session, "tabs", tab_ids[reactiveVals$current_tab + 1]) # here we update in case the new current tab is not a new max tab
  shinyjs::runjs("window.scrollTo(0, 0)")
})

# grow menu depending on current tab
output$sidebar <- renderUI({
  curr_menu <- sidebarMenu(id = "tabs",
                           tabs[1:reactiveVals$max_tab])
  # if(reactiveVals$current_tab == 6){
  #   curr_menu <- sidebarMenu(id = "tabs",
  #                            tabs[1:7])
  # }
  # here we update in case the new current tab is also a new max tab
  updateTabItems(session, "tabs", tab_ids[reactiveVals$max_tab])
  if (reactiveVals$max_tab > 1) reactiveVals$continue <- FALSE
  shinyjs::runjs("window.scrollTo(0, 0)")
  return(curr_menu)
})
