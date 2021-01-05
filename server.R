server <- function(input, output, session) {
  
  downloadPlotFunction <- function(name, ggplotObject, width = 7, height = 7){
    return(
      downloadHandler(
        filename = function(){
          paste0(name, ".pdf")
        },
        content = function(file){
          ggsave(file, plot = ggplotObject, width=width, height=height)
        }
      )
    )
  }
  
  # make reactiveValues and server-wide variables
  tab_ids <- c("welcome", "start", "preprocessing", "visualization", "clustering", "de")
  reactiveVals <- reactiveValues()
  reactiveVals$current_tab <- 1
  
  # read all server files
  sapply(list.files("server", full.names = TRUE), source, environment())
}


