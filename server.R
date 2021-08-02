server <- function(input, output, session) {
  
  # make reactiveValues and server-wide variables
  reactiveVals <- reactiveValues()
  
  reactiveVals$preprocessingShowed <- FALSE

  spinner <- list(logo = spin_loaders(id = 29, color = "black"), color="rgb(146, 180, 242, .5)")
  
  # read all server files
  sapply(list.files("server", full.names = TRUE), source, environment())
  sapply(list.files("functions", full.names = TRUE), source, environment())
  
  # additional functions
  #downloadPlotFunction <- function(name, ggplotObject, width = 7, height = 7){
  #  return(
  #    downloadHandler(
  #      filename = function(){
  #        paste0(name, ".pdf")
  #      },
  #      content = function(file){
  #        ggsave(file, plot = ggplotObject, width=width, height=height)
  #      }
  #    )
  #  )
  #}
  
  sortMarkerNames <- function(choices, classes, first = "type"){
    classes[classes == first] <- paste0("a", first)
    ret <- choices[order(factor(sprintf("%s %s", classes, choices)))]
    return(ret)
  }
  
  output$dashboard <- renderUI({
    req(reactiveVals$sce)
    downloadButton("dashboardButton", "Download current SCE object")
  })
  
  output$dashboardButton <- downloadHandler(
    filename = function(){
      paste("sce.rds")
    },
    content = function(file){
      waiter_show(id = "app", html = tagList(spinner$logo, 
                                            HTML("<br>Downloading...")), 
                  color=spinner$color)
      saveRDS(object = reactiveVals$sce, file = file)
      waiter_hide(id="app")
    }
  )
  
  output$licenseTable <- renderTable({
    deps_table <- sapply(unique(renv::dependencies()$Package), packageDescription, fields="License")
    data.frame(Package = names(deps_table), License = deps_table)
  })
  
  shinyBS::updateButton(session, inputId = "previousTab", icon = icon("arrow-left"), style = "success")
  shinyBS::updateButton(session, inputId = "nextTab", icon = icon("arrow-right"), style = "success", disabled = FALSE)
  shinyjs::hide("loading")
  waiter_hide()
}

