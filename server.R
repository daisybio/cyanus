server <- function(input, output, session) {
  
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
  
  # make reactiveValues and server-wide variables
  tab_ids <- c("welcome", "start", "preprocessing", "visualization", "clustering", "de")
  reactiveVals <- reactiveValues()
  reactiveVals$current_tab <- 1
  
  # read all server files
  sapply(list.files("server", full.names = TRUE), source, environment())
}

