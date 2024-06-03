server <- function(input, output, session) {
  library(data.table)
  
  # make reactiveValues and server-wide variables
  reactiveVals <- reactiveValues()
  
  reactiveVals$call_list <- list() # saves all calls in a list
  reactiveVals$preprocessingShowed <- FALSE
  reactiveVals$selected_palette <- c("#ff6db6", "#004949", "#db6d00",  "#B2DF8A", "#FDB462", "#490092", "#009999", "#8f4e00", "#ffdf4d", "#171723","#b66dff")
  spinner <- list(logo = list(a(icon('envelope'), " Contact", href = "https://github.com/biomedbigdata/cyanus/issues", target = "_blank", style='color: black'), 
                              tags$br(), tags$br(), 
                              spin_loaders(id = 5, color = "black")), color="rgb(146, 180, 242, .5)")
  
  # read all server files
  sapply(list.files("R/server", full.names = TRUE), source, environment())
  sapply(list.files("R/functions", full.names = TRUE), source, environment())
  
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
    div(downloadButton("dashboardButton", "Download current SCE object"), style="height: 50px; display: flex; align-items: center;")
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
  
  output$downloadLog <- renderUI({
    req(reactiveVals$sce)
    req(reactiveVals$call_list)
    div(downloadButton("downloadLogButton", "Download R Script"), style="height: 50px; display: flex; align-items: center; margin-right: 10px")
  })
  
  output$downloadLogButton <- downloadHandler(
    filename = function(){
      "log.R"
    },
    content = function(file){
      waiter_show(id = "log", html = tagList(spinner$logo, 
                                             HTML("<br>Downloading...")), 
                  color=spinner$color)
      
      session_dir <- file.path(tempdir(), paste0("session_", session$token))
      dir.create(session_dir, showWarnings = FALSE)
      log_file <- file.path(session_dir, "log.R")
      log_content <- paste(reactiveVals$call_list, collapse = '\n')
      writeLines(log_content, con = log_file)
      file.copy(log_file, file)
      waiter_hide(id="log")
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
  onStop(function() {
    session_dir <- file.path(tempdir(), paste0("session_", session$token))
    if (file.exists(session_dir)) {
      unlink(session_dir, recursive = TRUE)
    }})
}

