server <- function(input, output, session) {
  # make reactiveValues and server-wide variables
  reactiveVals <- reactiveValues()
  
  reactiveVals$preprocessingShowed <- FALSE
  
  # reactiveVals$enable_blacklist <- new.env()
  # my_enable <- function(id = NULL, selector = NULL, asis = FALSE) {
  #   rm(id, pos = reactiveVals$enable_blacklist)
  #   shinyjs::enable(id, selector, asis)
  # }
  # 
  # my_disable <- function(id = NULL, selector = NULL, asis = FALSE) {
  #   assign(id, TRUE, pos = reactiveVals$enable_blacklist)
  #   shinyjs::disable(id, selector, asis)
  # }
  # 
  # toggle_inputs <- function(enable_inputs = FALSE,
  #                           input_list = input)
  # {
  #   # Toggle elements 
  #   #also disables all downloadButtons automatically
  #   for (x in c(names(input_list), session$downloads$keys()))
  #     if (enable_inputs & !exists(x, reactiveVals$enable_blacklist)) {
  #       shinyjs::enable(x) # TODO: make reactive value with all ids that should not be enabled
  #     } else {
  #       shinyjs::disable(x)
  #     }
  # }
  
  # read all server files
  sapply(list.files("server", full.names = TRUE), source, environment())
  
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
      saveRDS(object = reactiveVals$sce, file = file)
    }
  )
  
  shinyBS::updateButton(session, inputId = "previousTab", label = " Previous", icon = icon("arrow-left"), style = "success")
  shinyBS::updateButton(session, inputId = "nextTab", label = " Next", icon = icon("arrow-right"), style = "success", disabled = FALSE)
  shinyjs::hide("loading")
}

